function [flatNeighbors, leafNeighbors] = generateNeighbors(f)

nboxes = length(f.boxes);
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

    if ( isLeaf(f.boxes(id)) )
        nbr = flatNeighbors(1, id);
        if ( ~isnan(nbr) )
            if ( isLeaf(f.boxes(nbr)) )
                leafNeighbors(1, id) = {nbr};
            else
                leafNeighbors(1, id) = {f.boxes(nbr).children([2 4])};
            end
        end
        
        nbr = flatNeighbors(2, id);
        if ( ~isnan(nbr) )
            if ( isLeaf(f.boxes(nbr)) )
                leafNeighbors(2, id) = {nbr};
            else
                leafNeighbors(2, id) = {f.boxes(nbr).children([1 3])};
            end
        end

        nbr = flatNeighbors(3, id);
        if ( ~isnan(nbr) )
            if ( isLeaf(f.boxes(nbr)) )
                leafNeighbors(3, id) = {nbr};
            else
                leafNeighbors(3, id) = {f.boxes(nbr).children([3 4])};
            end
        end
        
        nbr = flatNeighbors(4, id);
        if ( ~isnan(nbr) )
            if ( isLeaf(f.boxes(nbr)) )
                leafNeighbors(4, id) = {nbr};
            else
                leafNeighbors(4, id) = {f.boxes(nbr).children([1 2])};
            end
        end
        
        nbr = flatNeighbors(5, id);
        if ( ~isnan(nbr) )
            if ( isLeaf(f.boxes(nbr)) )
                leafNeighbors(5, id) = {nbr};
            else
                leafNeighbors(5, id) = {f.boxes(nbr).children(4)};
            end
        end
        
        nbr = flatNeighbors(6, id);
        if ( ~isnan(nbr) )
            if ( isLeaf(f.boxes(nbr)) )
                leafNeighbors(6, id) = {nbr};
            else
                leafNeighbors(6, id) = {f.boxes(nbr).children(3)};
            end
        end
        
        nbr = flatNeighbors(7, id);
        if ( ~isnan(nbr) )
            if ( isLeaf(f.boxes(nbr)) )
                leafNeighbors(7, id) = {nbr};
            else
                leafNeighbors(7, id) = {f.boxes(nbr).children(2)};
            end
        end
        
        nbr = flatNeighbors(8, id);
        if ( ~isnan(nbr) )
            if ( isLeaf(f.boxes(nbr)) )
                leafNeighbors(8, id) = {nbr};
            else
                leafNeighbors(8, id) = {f.boxes(nbr).children(1)};
            end
        end
    end

end

end

function out = leftNeighbor(f, id, flatNeighbors)

    if ( f.boxes(id).level == 0 )
        out = NaN;
        return
    end

    pid = f.boxes(id).parent;
    if ( id == f.boxes(pid).children(2) ) % SE
        out = f.boxes(pid).children(1);   % SW
        return
    end
    if ( id == f.boxes(pid).children(4) ) % NE
        out = f.boxes(pid).children(3);   % NW
        return
    end

    mu = flatNeighbors(1, pid);

    if ( isnan(mu) || isLeaf(f.boxes(mu)) )
        out = mu;
    else
        if ( id == f.boxes(pid).children(1) ) % SW
            out = f.boxes(mu).children(2);    % SE
        else
            out = f.boxes(mu).children(4);    % NE
        end
    end

end


function out = rightNeighbor(f, id, flatNeighbors)

    if ( f.boxes(id).level == 0 )
        out = NaN;
        return
    end

    pid = f.boxes(id).parent;
    if ( id == f.boxes(pid).children(1) ) % SW
        out = f.boxes(pid).children(2);   % SE
        return
    end
    if ( id == f.boxes(pid).children(3) ) % NW
        out = f.boxes(pid).children(4);   % NE
        return
    end

    mu = flatNeighbors(2, pid);

    if ( isnan(mu) || isLeaf(f.boxes(mu)) )
        out = mu;
    else
        if ( id == f.boxes(pid).children(2) ) % SE
            out = f.boxes(mu).children(1);    % SW
        else
            out = f.boxes(mu).children(3);    % NW
        end
    end

end

function out = downNeighbor(f, id, flatNeighbors)

    if ( f.boxes(id).level == 0 )
        out = NaN;
        return
    end

    pid = f.boxes(id).parent;
    if ( id == f.boxes(pid).children(3) ) % NW
        out = f.boxes(pid).children(1);   % SW
        return
    end
    if ( id == f.boxes(pid).children(4) ) % NE
        out = f.boxes(pid).children(2);   % SE
        return
    end

    mu = flatNeighbors(3, pid);

    if ( isnan(mu) || isLeaf(f.boxes(mu)) )
        out = mu;
    else
        if ( id == f.boxes(pid).children(1) ) % SW
            out = f.boxes(mu).children(3);    % NW
        else
            out = f.boxes(mu).children(4);    % NE
        end
    end

end

function out = upNeighbor(f, id, flatNeighbors)

    if ( f.boxes(id).level == 0 )
        out = NaN;
        return
    end

    pid = f.boxes(id).parent;
    if ( id == f.boxes(pid).children(1) ) % SW
        out = f.boxes(pid).children(3);   % NW
        return
    end
    if ( id == f.boxes(pid).children(2) ) % SE
        out = f.boxes(pid).children(4);   % NE
        return
    end

    mu = flatNeighbors(4, pid);

    if ( isnan(mu) || isLeaf(f.boxes(mu)) )
        out = mu;
    else
        if ( id == f.boxes(pid).children(3) ) % NW
            out = f.boxes(mu).children(1);    % SW
        else
            out = f.boxes(mu).children(2);    % SE
        end
    end

end

function out = leftDownNeighbor(f, id, flatNeighbors)

    if ( f.boxes(id).level == 0 )
        out = NaN;
        return
    end

    pid = f.boxes(id).parent;
    switch id
        case f.boxes(pid).children(1) % SW
            % mu = leftDownNeighbor(f, pid);
            mu = flatNeighbors(5, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = mu;
            else
                out = f.boxes(mu).children(4); % NE
            end
        case f.boxes(pid).children(2) % SE
            % mu = downNeighbor(f, pid);
            mu = flatNeighbors(3, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = NaN; % No corner neighbor
            else
                out = f.boxes(mu).children(3); % NW
            end
        case f.boxes(pid).children(3) % NW
            % mu = leftNeighbor(f, pid);
            mu = flatNeighbors(1, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = NaN; % No corner neighbor
            else
                out = f.boxes(mu).children(2); % SE
            end
        case f.boxes(pid).children(4) % NE
            out = f.boxes(pid).children(1); % SW
    end

end

function out = rightDownNeighbor(f, id, flatNeighbors)

    if ( f.boxes(id).level == 0 )
        out = NaN;
        return
    end

    pid = f.boxes(id).parent;
    switch id
        case f.boxes(pid).children(2) % SE
            % mu = rightDownNeighbor(f, pid);
            mu = flatNeighbors(6, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = mu;
            else
                out = f.boxes(mu).children(3); % NW
            end
        case f.boxes(pid).children(1) % SW
            % mu = downNeighbor(f, pid);
            mu = flatNeighbors(3, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = NaN; % No corner neighbor
            else
                out = f.boxes(mu).children(4); % NE
            end
        case f.boxes(pid).children(4) % NE
            % mu = rightNeighbor(f, pid);
            mu = flatNeighbors(2, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = NaN; % No corner neighbor
            else
                out = f.boxes(mu).children(1); % SW
            end
        case f.boxes(pid).children(3) % NW
            out = f.boxes(pid).children(2); % SE
    end

end

function out = leftUpNeighbor(f, id, flatNeighbors)

    if ( f.boxes(id).level == 0 )
        out = NaN;
        return
    end

    pid = f.boxes(id).parent;
    switch id
        case f.boxes(pid).children(3) % NW
            % mu = leftUpNeighbor(f, pid);
            mu = flatNeighbors(7, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = mu;
            else
                out = f.boxes(mu).children(2); % SE
            end
        case f.boxes(pid).children(1) % SW
            % mu = leftNeighbor(f, pid);
            mu = flatNeighbors(1, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = NaN; % No corner neighbor
            else
                out = f.boxes(mu).children(4); % NE
            end
        case f.boxes(pid).children(4) % NE
            % mu = upNeighbor(f, pid);
            mu = flatNeighbors(4, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = NaN; % No corner neighbor
            else
                out = f.boxes(mu).children(1); % SW
            end
        case f.boxes(pid).children(2) % SE
            out = f.boxes(pid).children(3); % NW
    end

end

function out = rightUpNeighbor(f, id, flatNeighbors)

    if ( f.boxes(id).level == 0 )
        out = NaN;
        return
    end

    pid = f.boxes(id).parent;
    switch id
        case f.boxes(pid).children(4) % NE
            % mu = rightUpNeighbor(f, pid);
            mu = flatNeighbors(8, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = mu;
            else
                out = f.boxes(mu).children(1); % SW
            end
        case f.boxes(pid).children(2) % SE
            % mu = rightNeighbor(f, pid);
            mu = flatNeighbors(2, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = NaN; % No corner neighbor
            else
                out = f.boxes(mu).children(3); % NW
            end
        case f.boxes(pid).children(3) % NW
            % mu = upNeighbor(f, pid);
            mu = flatNeighbors(4, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = NaN; % No corner neighbor
            else
                out = f.boxes(mu).children(2); % SE
            end
        case f.boxes(pid).children(1) % SW
            out = f.boxes(pid).children(4); % NE
    end

end
