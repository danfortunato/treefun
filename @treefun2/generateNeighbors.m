function [flatNeighbors, leafNeighbors] = generateNeighbors(f)

flevel    = f.level;
fparent   = f.parent;
fchildren = f.children;
fheight   = f.height;

nboxes = length(f.id);
flatNeighbors = zeros(8, nboxes);
leafNeighbors = cell(8, nboxes);
for id = 1:nboxes

    pid = fparent(id);

    % Left neighbor
    if ( flevel(id) == 0 )
        flatNeighbors(1, id) = NaN;
    elseif ( id == fchildren(2, pid) ) % SE
        flatNeighbors(1, id) = fchildren(1, pid); % SW
    elseif ( id == fchildren(4, pid) ) % NE
        flatNeighbors(1, id) = fchildren(3, pid); % NW
    else
        mu = flatNeighbors(1, pid);
        if ( isnan(mu) || fheight(mu) == 0 )
            flatNeighbors(1, id) = mu;
        elseif ( id == fchildren(1, pid) ) % SW
            flatNeighbors(1, id) = fchildren(2, mu); % SE
        else
            flatNeighbors(1, id) = fchildren(4, mu); % NE
        end
    end

    % Right neighbor
    if ( flevel(id) == 0 )
        flatNeighbors(2, id) = NaN;
    elseif ( id == fchildren(1, pid) ) % SW
        flatNeighbors(2, id) = fchildren(2, pid); % SE
    elseif ( id == fchildren(3, pid) ) % NW
        flatNeighbors(2, id) = fchildren(4, pid); % NE
    else
        mu = flatNeighbors(2, pid);
        if ( isnan(mu) || fheight(mu) == 0 )
            flatNeighbors(2, id) = mu;
        elseif ( id == fchildren(2, pid) ) % SE
            flatNeighbors(2, id) = fchildren(1, mu); % SW
        else
            flatNeighbors(2, id) = fchildren(3, mu); % NW
        end
    end

    % Down neighbor
    if ( flevel(id) == 0 )
        flatNeighbors(3, id) = NaN;
    elseif ( id == fchildren(3, pid) ) % NW
        flatNeighbors(3, id) = fchildren(1, pid); % SW
    elseif ( id == fchildren(4, pid) ) % NE
        flatNeighbors(3, id) = fchildren(2, pid); % SE
    else
        mu = flatNeighbors(3, pid);
        if ( isnan(mu) || fheight(mu) == 0 )
            flatNeighbors(3, id) = mu;
        elseif ( id == fchildren(1, pid) ) % SW
            flatNeighbors(3, id) = fchildren(3, mu); % NW
        else
            flatNeighbors(3, id) = fchildren(4, mu); % NE
        end
    end

    % Up neighbor
    if ( flevel(id) == 0 )
        flatNeighbors(4, id) = NaN;
    elseif ( id == fchildren(1, pid) ) % SW
        flatNeighbors(4, id) = fchildren(3, pid); % NW
    elseif ( id == fchildren(2, pid) ) % SE
        flatNeighbors(4, id) = fchildren(4, pid); % NE
    else
        mu = flatNeighbors(4, pid);
        if ( isnan(mu) || fheight(mu) == 0 )
            flatNeighbors(4, id) = mu;
        elseif ( id == fchildren(3, pid) ) % NW
            flatNeighbors(4, id) = fchildren(1, mu); % SW
        else
            flatNeighbors(4, id) = fchildren(2, mu); % SE
        end
    end

    % Left down neighbor
    if ( flevel(id) == 0 )
        flatNeighbors(5, id) = NaN;
    elseif ( id == fchildren(1, pid) ) % SW
        mu = flatNeighbors(5, pid);
        if ( isnan(mu) || fheight(mu) == 0 )
            flatNeighbors(5, id) = mu;
        else
            flatNeighbors(5, id) = fchildren(4, mu); % NE
        end
    elseif ( id == fchildren(2, pid) ) % SE
        mu = flatNeighbors(3, pid);
        if ( isnan(mu) || fheight(mu) == 0 )
            flatNeighbors(5, id) = NaN; % No corner neighbor
        else
            flatNeighbors(5, id) = fchildren(3, mu); % NW
        end
    elseif ( id == fchildren(3, pid) ) % NW
        mu = flatNeighbors(1, pid);
        if ( isnan(mu) || fheight(mu) == 0 )
            flatNeighbors(5, id) = NaN; % No corner neighbor
        else
            flatNeighbors(5, id) = fchildren(2, mu); % SE
        end
    elseif ( id == fchildren(4, pid) ) % NE
        flatNeighbors(5, id) = fchildren(1, pid); % SW
    end

    % Right down neighbor
    if ( flevel(id) == 0 )
        flatNeighbors(6, id) = NaN;
    elseif ( id == fchildren(2, pid) ) % SE
        mu = flatNeighbors(6, pid);
        if ( isnan(mu) || fheight(mu) == 0 )
            flatNeighbors(6, id) = mu;
        else
            flatNeighbors(6, id) = fchildren(3, mu); % NW
        end
    elseif ( id == fchildren(1, pid) ) % SW
        mu = flatNeighbors(3, pid);
        if ( isnan(mu) || fheight(mu) == 0 )
            flatNeighbors(6, id) = NaN; % No corner neighbor
        else
            flatNeighbors(6, id) = fchildren(4, mu); % NE
        end
    elseif ( id == fchildren(4, pid) ) % NE
        mu = flatNeighbors(2, pid);
        if ( isnan(mu) || fheight(mu) == 0 )
            flatNeighbors(6, id) = NaN; % No corner neighbor
        else
            flatNeighbors(6, id) = fchildren(1, mu); % SW
        end
    elseif ( id == fchildren(3, pid) ) % NW
        flatNeighbors(6, id) = fchildren(2, pid); % SE
    end

    % Left up neighbor
    if ( flevel(id) == 0 )
        flatNeighbors(7, id) = NaN;
    elseif ( id == fchildren(3, pid) ) % NW
        mu = flatNeighbors(7, pid);
        if ( isnan(mu) || fheight(mu) == 0 )
            flatNeighbors(7, id) = mu;
        else
            flatNeighbors(7, id) = fchildren(2, mu); % SE
        end
    elseif ( id == fchildren(1, pid) ) % SW
        mu = flatNeighbors(1, pid);
        if ( isnan(mu) || fheight(mu) == 0 )
            flatNeighbors(7, id) = NaN; % No corner neighbor
        else
            flatNeighbors(7, id) = fchildren(4, mu); % NE
        end
    elseif ( id == fchildren(4, pid) ) % NE
        mu = flatNeighbors(4, pid);
        if ( isnan(mu) || fheight(mu) == 0 )
            flatNeighbors(7, id) = NaN; % No corner neighbor
        else
            flatNeighbors(7, id) = fchildren(1, mu); % SW
        end
    elseif ( id == fchildren(2, pid) ) % SE
        flatNeighbors(7, id) = fchildren(3, pid); % NW
    end

    % Right up neighbor
    if ( flevel(id) == 0 )
        flatNeighbors(8, id) = NaN;
    elseif ( id == fchildren(4, pid) ) % NE
        mu = flatNeighbors(8, pid);
        if ( isnan(mu) || fheight(mu) == 0 )
            flatNeighbors(8, id) = mu;
        else
            flatNeighbors(8, id) = fchildren(1, mu); % SW
        end
    elseif ( id == fchildren(2, pid) ) % SE
        mu = flatNeighbors(2, pid);
        if ( isnan(mu) || fheight(mu) == 0 )
            flatNeighbors(8, id) = NaN; % No corner neighbor
        else
            flatNeighbors(8, id) = fchildren(3, mu); % NW
        end
    elseif ( id == fchildren(3, pid) ) % NW
        mu = flatNeighbors(4, pid);
        if ( isnan(mu) || fheight(mu) == 0 )
            flatNeighbors(8, id) = NaN; % No corner neighbor
        else
            flatNeighbors(8, id) = fchildren(2, mu); % SE
        end
    elseif ( id == fchildren(1, pid) ) % SW
        flatNeighbors(8, id) = fchildren(4, pid); % NE
    end

    if ( fheight(id) == 0 )

        nbr = flatNeighbors(1, id);
        if ( ~isnan(nbr) )
            if ( fheight(nbr) == 0 )
                leafNeighbors{1, id} = nbr;
            else
                leafNeighbors{1, id} = fchildren([2 4], nbr);
            end
        end

        nbr = flatNeighbors(2, id);
        if ( ~isnan(nbr) )
            if ( fheight(nbr) == 0 )
                leafNeighbors{2, id} = nbr;
            else
                leafNeighbors{2, id} = fchildren([1 3], nbr);
            end
        end

        nbr = flatNeighbors(3, id);
        if ( ~isnan(nbr) )
            if ( fheight(nbr) == 0 )
                leafNeighbors{3, id} = nbr;
            else
                leafNeighbors{3, id} = fchildren([3 4], nbr);
            end
        end

        nbr = flatNeighbors(4, id);
        if ( ~isnan(nbr) )
            if ( fheight(nbr) == 0 )
                leafNeighbors{4, id} = nbr;
            else
                leafNeighbors{4, id} = fchildren([1 2], nbr);
            end
        end

        nbr = flatNeighbors(5, id);
        if ( ~isnan(nbr) )
            if ( fheight(nbr) == 0 )
                leafNeighbors{5, id} = nbr;
            else
                leafNeighbors{5, id} = fchildren(4, nbr);
            end
        end

        nbr = flatNeighbors(6, id);
        if ( ~isnan(nbr) )
            if ( fheight(nbr) == 0 )
                leafNeighbors{6, id} = nbr;
            else
                leafNeighbors{6, id} = fchildren(3, nbr);
            end
        end

        nbr = flatNeighbors(7, id);
        if ( ~isnan(nbr) )
            if ( fheight(nbr) == 0 )
                leafNeighbors{7, id} = nbr;
            else
                leafNeighbors{7, id} = fchildren(2, nbr);
            end
        end

        nbr = flatNeighbors(8, id);
        if ( ~isnan(nbr) )
            if ( fheight(nbr) == 0 )
                leafNeighbors{8, id} = nbr;
            else
                leafNeighbors{8, id} = fchildren(1, nbr);
            end
        end

    end

end

end
