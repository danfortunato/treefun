function f = balance(f)

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

    if ( ~split && ~isnan(neighborIDs(6)) ) % Right down corner
        nbr = neighborIDs(6);
        if ( ~isLeaf(f, nbr) )
            childrenIDs = f.children(3, nbr); % Left up
            split = split | any( ~isLeaf(f, childrenIDs) );
        end
    end

    if ( ~split && ~isnan(neighborIDs(7)) ) % Left up corner
        nbr = neighborIDs(7);
        if ( ~isLeaf(f, nbr) )
            childrenIDs = f.children(2, nbr); % Right down
            split = split | any( ~isLeaf(f, childrenIDs) );
        end
    end

    if ( ~split && ~isnan(neighborIDs(8)) ) % Right up corner
        nbr = neighborIDs(8);
        if ( ~isLeaf(f, nbr) )
            childrenIDs = f.children(1, nbr); % Left down
            split = split | any( ~isLeaf(f, childrenIDs) );
        end
    end

    if ( split )
        % This was a leaf, so we'll use its coeffs to evaluate
        % on the new children
        coeffs = f.coeffs{id};
        % Split into four child boxes
        f = refineBox(f, id);
        children = f.children(:,id);
        [LL, LR, UL, UR] = coeffs2children(coeffs);
        f.coeffs{children(1)} = LL; % Lower left
        f.coeffs{children(2)} = LR; % Lower right
        f.coeffs{children(3)} = UL; % Upper left
        f.coeffs{children(4)} = UR; % Upper right
        L = [L children(:).']; %#ok<AGROW>
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

end
