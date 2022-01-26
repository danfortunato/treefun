function f = balance_v2(f, func)

L = leaves(f);
while ( ~isempty(L) )
    id = L(1);
    neighborIDs = neighbors(f, id);
    split = false;

    if ( ~split && ~isnan(neighborIDs(1)) ) % Left neighbor
        nbr = neighborIDs(1);
        if ( ~isLeaf(f, nbr) )
            childrenIDs = f.children([2 4], nbr); % SE, NE
            split = split | any( ~isLeaf(f, childrenIDs) );
        end
    end

    if ( ~split && ~isnan(neighborIDs(2)) ) % Right neighbor
        nbr = neighborIDs(2);
        if ( ~isLeaf(f, nbr) )
            childrenIDs = f.children([1 3], nbr); % SW, NW
            split = split | any( ~isLeaf(f, childrenIDs) );
        end
    end

    if ( ~split && ~isnan(neighborIDs(3)) ) % Down neighbor
        nbr = neighborIDs(3);
        if ( ~isLeaf(f, nbr) )
            childrenIDs = f.children([3 4], nbr); % NW, NE
            split = split | any( ~isLeaf(f, childrenIDs) );
        end
    end

    if ( ~split && ~isnan(neighborIDs(4)) ) % Up neighbor
        nbr = neighborIDs(4);
        if ( ~isLeaf(f, nbr) )
            childrenIDs = f.children([1 2], nbr); % SW, SE
            split = split | any( ~isLeaf(f, childrenIDs) );
        end
    end

    if ( ~split && ~isnan(neighborIDs(5)) ) % Left down corner
        nbr = neighborIDs(5);
        if ( ~isLeaf(f, nbr) )
            childrenIDs = f.children(4, nbr); % Right up
            split = split | any( ~isLeaf(f, childrenIDs) );
        end
    end

    if ( split )
        % Split into four child boxes
        f = refineBox(f, id);
        L = [L f.children(:,id).']; %#ok<AGROW>
        neighborIDs = neighborIDs( ~isnan(neighborIDs) );
        idx = f.level(neighborIDs) < f.level(id);
        L = [L neighborIDs(idx)]; %#ok<AGROW>
        L = unique(L, 'stable');
    end

    L(1) = [];
end

% Do a cumulative sum in reverse to correct the heights
for k = length(f.id):-1:1
    if ( ~isLeaf(f, k) )
        f.height(k) = 1 + max(f.height(f.children(:,k)));
    end
end

% Now fill in the leaf coefficients
L = leaves(f);
[xx0, yy0] = chebpts2(f.n, f.n, [0 1 0 1]);
for id = L(:).'
    if ( isempty(f.coeffs{id}) )
        dom = f.domain(:,id);
        sclx = diff(dom(1:2));
        scly = diff(dom(3:4));
        xx = sclx*xx0 + dom(1);
        yy = scly*yy0 + dom(3);
        vals = func(xx,yy);
        f.coeffs{id} = treefun2.vals2coeffs(vals);
    end
end

end
