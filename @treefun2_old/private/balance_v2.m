function f = balance_v2(f, func)

L = [leaves(f).id];
while ( ~isempty(L) )
    id = L(1);
    neighborIDs = neighbors(f, id);
    split = false;

    if ( ~split && ~isnan(neighborIDs(1)) ) % Left neighbor
        nbr = f.boxes(neighborIDs(1));
        if ( ~isLeaf(nbr) )
            childrenIDs = nbr.children([2 4]); % SE, NE
            split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
        end
    end

    if ( ~split && ~isnan(neighborIDs(2)) ) % Right neighbor
        nbr = f.boxes(neighborIDs(2));
        if ( ~isLeaf(nbr) )
            childrenIDs = nbr.children([1 3]); % SW, NW
            split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
        end
    end

    if ( ~split && ~isnan(neighborIDs(3)) ) % Down neighbor
        nbr = f.boxes(neighborIDs(3));
        if ( ~isLeaf(nbr) )
            childrenIDs = nbr.children([3 4]); % NW, NE
            split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
        end
    end

    if ( ~split && ~isnan(neighborIDs(4)) ) % Up neighbor
        nbr = f.boxes(neighborIDs(4));
        if ( ~isLeaf(nbr) )
            childrenIDs = nbr.children([1 2]); % SW, SE
            split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
        end
    end

    if ( ~split && ~isnan(neighborIDs(5)) ) % Left down corner
        nbr = f.boxes(neighborIDs(5));
        if ( ~isLeaf(nbr) )
            childrenIDs = nbr.children(4); % Right up
            split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
        end
    end

    if ( split )
        % Split into four child boxes
        f = refineBox(f, id);
        L = [L f.boxes(id).children]; %#ok<AGROW>
        neighborIDs = neighborIDs( ~isnan(neighborIDs) );
        idx = [f.boxes(neighborIDs).level] < f.boxes(id).level;
        L = [L neighborIDs(idx)]; %#ok<AGROW>
        L = unique(L, 'stable');
    end

    L(1) = [];
end

% Do a cumulative sum in reverse to correct the heights
for k = length(f.boxes):-1:1
    if ( ~isLeaf(f.boxes(k)) )
        f.boxes(k).height = 1 + max([f.boxes(f.boxes(k).children).height]);
    end
end

% Now fill in the leaf coefficients
L = [leaves(f).id];
[xx0, yy0] = chebpts2(f.n, f.n, [0 1 0 1]);
for k = 1:length(L)
    id = L(k);
    if ( isempty(f.boxes(id).coeffs) )
        dom = f.boxes(id).domain;
        sclx = diff(dom(1:2));
        scly = diff(dom(3:4));
        xx = sclx*xx0 + dom(1);
        yy = scly*yy0 + dom(3);
        vals = func(xx,yy);
        f.boxes(id).coeffs = treefun2.vals2coeffs(vals);
    end
end

end
