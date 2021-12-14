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
        out = isNeighbor(f.boxes(id), node, side) || ...
              isParent(f.boxes(id), node);
    end

    function accum(node)
        v(end+1) = node.id;
    end

end

function search(f, cond, accum, id)

    if ( cond(f.boxes(id)) )
        if ( ~isLeaf(f.boxes(id)) )
            search(f, cond, accum, f.boxes(id).children(1));
            search(f, cond, accum, f.boxes(id).children(2));
            search(f, cond, accum, f.boxes(id).children(3));
            search(f, cond, accum, f.boxes(id).children(4));
        else
            accum(f.boxes(id));
        end
    end

end

function out = leftNeighbor(f, id)

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

    mu = leftNeighbor(f, pid);

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

function out = rightNeighbor(f, id)

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

    mu = rightNeighbor(f, pid);

    %if ( isnan(mu) || isLeaf(f.boxes(mu)) )
    ia = isnan(mu);
    if ( ia )
        out = mu;
    else
        mid = f.boxes(mu);
        ib = isLeaf(mid);
        if ( ib )
            out = mu;
        else
            if ( id == f.boxes(pid).children(2) ) % SE
                out = f.boxes(mu).children(1);    % SW
            else
                out = f.boxes(mu).children(3);    % NW
            end
        end
    end
end

function out = downNeighbor(f, id)

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

    mu = downNeighbor(f, pid);

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

function out = upNeighbor(f, id)

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

    mu = upNeighbor(f, pid);

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

function out = leftDownNeighbor(f, id)

    if ( f.boxes(id).level == 0 )
        out = NaN;
        return
    end

    pid = f.boxes(id).parent;
    switch id
        case f.boxes(pid).children(1) % SW
            mu = leftDownNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = mu;
            else
                out = f.boxes(mu).children(4); % NE
            end
        case f.boxes(pid).children(2) % SE
            mu = downNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = NaN; % No corner neighbor
            else
                out = f.boxes(mu).children(3); % NW
            end
        case f.boxes(pid).children(3) % NW
            mu = leftNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = NaN; % No corner neighbor
            else
                out = f.boxes(mu).children(2); % SE
            end
        case f.boxes(pid).children(4) % NE
            out = f.boxes(pid).children(1); % SW
    end

end

function out = rightDownNeighbor(f, id)

    if ( f.boxes(id).level == 0 )
        out = NaN;
        return
    end

    pid = f.boxes(id).parent;
    switch id
        case f.boxes(pid).children(2) % SE
            mu = rightDownNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = mu;
            else
                out = f.boxes(mu).children(3); % NW
            end
        case f.boxes(pid).children(1) % SW
            mu = downNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = NaN; % No corner neighbor
            else
                out = f.boxes(mu).children(4); % NE
            end
        case f.boxes(pid).children(4) % NE
            mu = rightNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = NaN; % No corner neighbor
            else
                out = f.boxes(mu).children(1); % SW
            end
        case f.boxes(pid).children(3) % NW
            out = f.boxes(pid).children(2); % SE
    end

end

function out = leftUpNeighbor(f, id)

    if ( f.boxes(id).level == 0 )
        out = NaN;
        return
    end

    pid = f.boxes(id).parent;
    switch id
        case f.boxes(pid).children(3) % NW
            mu = leftUpNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = mu;
            else
                out = f.boxes(mu).children(2); % SE
            end
        case f.boxes(pid).children(1) % SW
            mu = leftNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = NaN; % No corner neighbor
            else
                out = f.boxes(mu).children(4); % NE
            end
        case f.boxes(pid).children(4) % NE
            mu = upNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = NaN; % No corner neighbor
            else
                out = f.boxes(mu).children(1); % SW
            end
        case f.boxes(pid).children(2) % SE
            out = f.boxes(pid).children(3); % NW
    end

end

function out = rightUpNeighbor(f, id)

    if ( f.boxes(id).level == 0 )
        out = NaN;
        return
    end

    pid = f.boxes(id).parent;
    switch id
        case f.boxes(pid).children(4) % NE
            mu = rightUpNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = mu;
            else
                out = f.boxes(mu).children(1); % SW
            end
        case f.boxes(pid).children(2) % SE
            mu = rightNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = NaN; % No corner neighbor
            else
                out = f.boxes(mu).children(3); % NW
            end
        case f.boxes(pid).children(3) % NW
            mu = upNeighbor(f, pid);
            if ( isnan(mu) || isLeaf(f.boxes(mu)) )
                out = NaN; % No corner neighbor
            else
                out = f.boxes(mu).children(2); % SE
            end
        case f.boxes(pid).children(1) % SW
            out = f.boxes(pid).children(4); % NE
    end

end
