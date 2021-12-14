function neighbors = neighbors(f, id, type)

narginchk(2, 3);
if ( nargin < 3 )
    type = 'flat';
end

if ( strcmpi(type, 'leaf') )
    
    neighbors = cell(8, 1);

    v = [];
    search(f, @(node) cond(node, 'left'),  @accum, 1);
    neighbors{1} = v;

    v = [];
    search(f, @(node) cond(node, 'right'), @accum, 1);
    neighbors{2} = v;

    v = [];
    search(f, @(node) cond(node, 'down'),  @accum, 1);
    neighbors{3} = v;

    v = [];
    search(f, @(node) cond(node, 'up'),    @accum, 1);
    neighbors{4} = v;

    start = leftDownNeighbor(f, id);
    if ( ~isnan(start) )
        v = [];
        search(f, @(node) cond(node, 'leftdown'),  @accum, start);
        neighbors{5} = v;
    end

    start = rightDownNeighbor(f, id);
    if ( ~isnan(start) )
        v = [];
        search(f, @(node) cond(node, 'rightdown'), @accum, start);
        neighbors{6} = v;
    end

    start = leftUpNeighbor(f, id);
    if ( ~isnan(start) )
        v = [];
        search(f, @(node) cond(node, 'leftup'),    @accum, start);
        neighbors{7} = v;
    end

    start = rightUpNeighbor(f, id);
    if ( ~isnan(start) )
        v = [];
        search(f, @(node) cond(node, 'rightup'),   @accum, start);
        neighbors{8} = v;
    end

else

    neighbors = zeros(1, 8);
    neighbors(1) =      leftNeighbor(f, id);
    neighbors(2) =     rightNeighbor(f, id);
    neighbors(3) =      downNeighbor(f, id);
    neighbors(4) =        upNeighbor(f, id);
    neighbors(5) =  leftDownNeighbor(f, id);
    neighbors(6) = rightDownNeighbor(f, id);
    neighbors(7) =    leftUpNeighbor(f, id);
    neighbors(8) =   rightUpNeighbor(f, id);

end

    function out = cond(node, side)
        out = isNeighbor(f, id, node, side) || ...
              isParent(f, id, node);
    end

    function accum(node)
        v(end+1) = node;
    end

end

function search(f, cond, accum, id)

    if ( cond(id) )
        if ( ~isLeaf(f, id) )
            search(f, cond, accum, f.children(1, id));
            search(f, cond, accum, f.children(2, id));
            search(f, cond, accum, f.children(3, id));
            search(f, cond, accum, f.children(4, id));
        else
            accum(id);
        end
    end

end

function out = leftNeighbor(f, id)

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

    mu = leftNeighbor(f, pid);

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

function out = rightNeighbor(f, id)

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

    mu = rightNeighbor(f, pid);

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

function out = downNeighbor(f, id)

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

    mu = downNeighbor(f, pid);

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

function out = upNeighbor(f, id)

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

    mu = upNeighbor(f, pid);

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


function out = leftDownNeighbor(f, id)

    if ( f.level(id) == 0 )
        out = NaN;
        return
    end

    pid = f.parent(id);
    switch id
        case f.children(1, pid) % SW
            mu = leftDownNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = mu;
            else
                out = f.children(4, mu); % NE
            end
        case f.children(2, pid) % SE
            mu = downNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = NaN; % No corner neighbor
            else
                out = f.children(3, mu); % NW
            end
        case f.children(3, pid) % NW
            mu = leftNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = NaN; % No corner neighbor
            else
                out = f.children(2, mu); % SE
            end
        case f.children(4, pid) % NE
            out = f.children(1, pid); % SW
    end

end

function out = rightDownNeighbor(f, id)

    if ( f.level(id) == 0 )
        out = NaN;
        return
    end

    pid = f.parent(id);
    switch id
        case f.children(2, pid) % SE
            mu = rightDownNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = mu;
            else
                out = f.children(3, mu); % NW
            end
        case f.children(1, pid) % SW
            mu = downNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = NaN; % No corner neighbor
            else
                out = f.children(4, mu); % NE
            end
        case f.children(4, pid) % NE
            mu = rightNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = NaN; % No corner neighbor
            else
                out = f.children(1, mu); % SW
            end
        case f.children(3, pid) % NW
            out = f.children(2, pid); % SE
    end

end

function out = leftUpNeighbor(f, id)

    if ( f.level(id) == 0 )
        out = NaN;
        return
    end

    pid = f.parent(id);
    switch id
        case f.children(3, pid) % NW
            mu = leftUpNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = mu;
            else
                out = f.children(2, mu); % SE
            end
        case f.children(1, pid) % SW
            mu = leftNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = NaN; % No corner neighbor
            else
                out = f.children(4, mu); % NE
            end
        case f.children(4, pid) % NE
            mu = upNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = NaN; % No corner neighbor
            else
                out = f.children(1, mu); % SW
            end
        case f.children(2, pid) % SE
            out = f.children(3, pid); % NW
    end

end

function out = rightUpNeighbor(f, id)

    if ( f.level(id) == 0 )
        out = NaN;
        return
    end

    pid = f.parent(id);
    switch id
        case f.children(4, pid) % NE
            mu = rightUpNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = mu;
            else
                out = f.children(1, mu); % SW
            end
        case f.children(2, pid) % SE
            mu = rightNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = NaN; % No corner neighbor
            else
                out = f.children(3, mu); % NW
            end
        case f.children(3, pid) % NW
            mu = upNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f, mu) )
                out = NaN; % No corner neighbor
            else
                out = f.children(2, mu); % SE
            end
        case f.children(1, pid) % SW
            out = f.children(4, pid); % NE
    end

end
