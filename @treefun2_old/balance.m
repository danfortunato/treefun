function f = balance(f)

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

    if ( ~split && ~isnan(neighborIDs(6)) ) % Right down corner
        nbr = f.boxes(neighborIDs(6));
        if ( ~isLeaf(nbr) )
            childrenIDs = nbr.children(3); % Left up
            split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
        end
    end

    if ( ~split && ~isnan(neighborIDs(7)) ) % Left up corner
        nbr = f.boxes(neighborIDs(7));
        if ( ~isLeaf(nbr) )
            childrenIDs = nbr.children(2); % Right down
            split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
        end
    end

    if ( ~split && ~isnan(neighborIDs(8)) ) % Right up corner
        nbr = f.boxes(neighborIDs(8));
        if ( ~isLeaf(nbr) )
            childrenIDs = nbr.children(1); % Left down
            split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
        end
    end

    if ( split )
        % This was a leaf, so we'll use its coeffs to evaluate
        % on the new children
        coeffs = f.boxes(id).coeffs;
        % Split into four child boxes
        f = refineBox(f, id);
        children = f.boxes(id).children;
        [LL, LR, UL, UR] = coeffs2children(coeffs);
        f.boxes(children(1)).coeffs = LL; % Lower left
        f.boxes(children(2)).coeffs = LR; % Lower right
        f.boxes(children(3)).coeffs = UL; % Upper left
        f.boxes(children(4)).coeffs = UR; % Upper right
        L = [L children]; %#ok<AGROW>
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

end