function [flatNeighbors, leafNeighbors] = generateNeighbors(f)

nboxes = length(f.id);
flatNeighbors = zeros(8, nboxes);
leafNeighbors = cell(8, nboxes);
for id = 1:nboxes

    flatNeighbors(1, id) =      leftNeighbor(f, id, flatNeighbors);
    flatNeighbors(2, id) =     rightNeighbor(f, id, flatNeighbors);
    flatNeighbors(3, id) =      downNeighbor(f, id, flatNeighbors);
    flatNeighbors(4, id) =        upNeighbor(f, id, flatNeighbors);
    flatNeighbors(5, id) =  leftDownNeighbor(f, id, flatNeighbors);
    flatNeighbors(6, id) = rightDownNeighbor(f, id, flatNeighbors);
    flatNeighbors(7, id) =    leftUpNeighbor(f, id, flatNeighbors);
    flatNeighbors(8, id) =   rightUpNeighbor(f, id, flatNeighbors);

    if ( isLeaf(f, id) )
        nbr = flatNeighbors(1, id);
        if ( ~isnan(nbr) )
            if ( isLeaf(f, nbr) )
                leafNeighbors(1, id) = {nbr};
            else
                leafNeighbors(1, id) = {f.children([2 4], nbr)};
            end
        end
        
        nbr = flatNeighbors(2, id);
        if ( ~isnan(nbr) )
            if ( isLeaf(f, nbr) )
                leafNeighbors(2, id) = {nbr};
            else
                leafNeighbors(2, id) = {f.children([1 3], nbr)};
            end
        end

        nbr = flatNeighbors(3, id);
        if ( ~isnan(nbr) )
            if ( isLeaf(f, nbr) )
                leafNeighbors(3, id) = {nbr};
            else
                leafNeighbors(3, id) = {f.children([3 4], nbr)};
            end
        end
        
        nbr = flatNeighbors(4, id);
        if ( ~isnan(nbr) )
            if ( isLeaf(f, nbr) )
                leafNeighbors(4, id) = {nbr};
            else
                leafNeighbors(4, id) = {f.children([1 2], nbr)};
            end
        end
        
        nbr = flatNeighbors(5, id);
        if ( ~isnan(nbr) )
            if ( isLeaf(f, nbr) )
                leafNeighbors(5, id) = {nbr};
            else
                leafNeighbors(5, id) = {f.children(4, nbr)};
            end
        end
        
        nbr = flatNeighbors(6, id);
        if ( ~isnan(nbr) )
            if ( isLeaf(f, nbr) )
                leafNeighbors(6, id) = {nbr};
            else
                leafNeighbors(6, id) = {f.children(3, nbr)};
            end
        end
        
        nbr = flatNeighbors(7, id);
        if ( ~isnan(nbr) )
            if ( isLeaf(f, nbr) )
                leafNeighbors(7, id) = {nbr};
            else
                leafNeighbors(7, id) = {f.children(2, nbr)};
            end
        end
        
        nbr = flatNeighbors(8, id);
        if ( ~isnan(nbr) )
            if ( isLeaf(f, nbr) )
                leafNeighbors(8, id) = {nbr};
            else
                leafNeighbors(8, id) = {f.children(1, nbr)};
            end
        end
    end

end

end

function out = leftNeighbor(f, id, flatNeighbors)

    if ( f.level(id) == 0 )
        out = NaN;
        return
    end

    pid = f.parent(id);
    if ( id == f.children(2, pid) ) % SE
        out = f.children(1, pid);   % SW
        return
    end
    if ( id == f.children(4, pid) ) % NE
        out = f.children(3, pid);   % NW
        return
    end

    mu = flatNeighbors(1, pid);

    if ( isnan(mu) || isLeaf(f, mu) )
        out = mu;
    else
        if ( id == f.children(1, pid) ) % SW
            out = f.children(2, mu);    % SE
        else
            out = f.children(4, mu);    % NE
        end
    end

end


function out = rightNeighbor(f, id, flatNeighbors)

    if ( f.level(id) == 0 )
        out = NaN;
        return
    end

    pid = f.parent(id);
    if ( id == f.children(1, pid) ) % SW
        out = f.children(2, pid);   % SE
        return
    end
    if ( id == f.children(3, pid) ) % NW
        out = f.children(4, pid);   % NE
        return
    end

    mu = flatNeighbors(2, pid);

    if ( isnan(mu) || isLeaf(f, mu) )
        out = mu;
    else
        if ( id == f.children(2, pid) ) % SE
            out = f.children(1, mu);    % SW
        else
            out = f.children(3, mu);    % NW
        end
    end

end

function out = downNeighbor(f, id, flatNeighbors)

    if ( f.level(id) == 0 )
        out = NaN;
        return
    end

    pid = f.parent(id);
    if ( id == f.children(3, pid) ) % NW
        out = f.children(1, pid);   % SW
        return
    end
    if ( id == f.children(4, pid) ) % NE
        out = f.children(2, pid);   % SE
        return
    end

    mu = flatNeighbors(3, pid);

    if ( isnan(mu) || isLeaf(f, mu) )
        out = mu;
    else
        if ( id == f.children(1, pid) ) % SW
            out = f.children(3, mu);    % NW
        else
            out = f.children(4, mu);    % NE
        end
    end

end

function out = upNeighbor(f, id, flatNeighbors)

    if ( f.level(id) == 0 )
        out = NaN;
        return
    end

    pid = f.parent(id);
    if ( id == f.children(1, pid) ) % SW
        out = f.children(3, pid);   % NW
        return
    end
    if ( id == f.children(2, pid) ) % SE
        out = f.children(4, pid);   % NE
        return
    end

    mu = flatNeighbors(4, pid);

    if ( isnan(mu) || isLeaf(f, mu) )
        out = mu;
    else
        if ( id == f.children(3, pid) ) % NW
            out = f.children(1, mu);    % SW
        else
            out = f.children(2, mu);    % SE
        end
    end

end

function out = leftDownNeighbor(f, id, flatNeighbors)

    if ( f.level(id) == 0 )
        out = NaN;
        return
    end

    pid = f.parent(id);
    switch id
        case f.children(1, pid) % SW
            % mu = leftDownNeighbor(f, pid);
            mu = flatNeighbors(5, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = mu;
            else
                out = f.children(4, mu); % NE
            end
        case f.children(2, pid) % SE
            % mu = downNeighbor(f, pid);
            mu = flatNeighbors(3, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = NaN; % No corner neighbor
            else
                out = f.children(3, mu); % NW
            end
        case f.children(3, pid) % NW
            % mu = leftNeighbor(f, pid);
            mu = flatNeighbors(1, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = NaN; % No corner neighbor
            else
                out = f.children(2, mu); % SE
            end
        case f.children(4, pid) % NE
            out = f.children(1, pid); % SW
    end

end

function out = rightDownNeighbor(f, id, flatNeighbors)

    if ( f.level(id) == 0 )
        out = NaN;
        return
    end

    pid = f.parent(id);
    switch id
        case f.children(2, pid) % SE
            % mu = rightDownNeighbor(f, pid);
            mu = flatNeighbors(6, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = mu;
            else
                out = f.children(3, mu); % NW
            end
        case f.children(1, pid) % SW
            % mu = downNeighbor(f, pid);
            mu = flatNeighbors(3, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = NaN; % No corner neighbor
            else
                out = f.children(4, mu); % NE
            end
        case f.children(4, pid) % NE
            % mu = rightNeighbor(f, pid);
            mu = flatNeighbors(2, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = NaN; % No corner neighbor
            else
                out = f.children(1, mu); % SW
            end
        case f.children(3, pid) % NW
            out = f.children(2, pid); % SE
    end

end

function out = leftUpNeighbor(f, id, flatNeighbors)

    if ( f.level(id) == 0 )
        out = NaN;
        return
    end

    pid = f.parent(id);
    switch id
        case f.children(3, pid) % NW
            % mu = leftUpNeighbor(f, pid);
            mu = flatNeighbors(7, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = mu;
            else
                out = f.children(2, mu); % SE
            end
        case f.children(1, pid) % SW
            % mu = leftNeighbor(f, pid);
            mu = flatNeighbors(1, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = NaN; % No corner neighbor
            else
                out = f.children(4, mu); % NE
            end
        case f.children(4, pid) % NE
            % mu = upNeighbor(f, pid);
            mu = flatNeighbors(4, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = NaN; % No corner neighbor
            else
                out = f.children(1, mu); % SW
            end
        case f.children(2, pid) % SE
            out = f.children(3, pid); % NW
    end

end

function out = rightUpNeighbor(f, id, flatNeighbors)

    if ( f.level(id) == 0 )
        out = NaN;
        return
    end

    pid = f.parent(id);
    switch id
        case f.children(4, pid) % NE
            % mu = rightUpNeighbor(f, pid);
            mu = flatNeighbors(8, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = mu;
            else
                out = f.children(1, mu); % SW
            end
        case f.children(2, pid) % SE
            % mu = rightNeighbor(f, pid);
            mu = flatNeighbors(2, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = NaN; % No corner neighbor
            else
                out = f.children(3, mu); % NW
            end
        case f.children(3, pid) % NW
            % mu = upNeighbor(f, pid);
            mu = flatNeighbors(4, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = NaN; % No corner neighbor
            else
                out = f.children(2, mu); % SE
            end
        case f.children(1, pid) % SW
            out = f.children(4, pid); % NE
    end

end
