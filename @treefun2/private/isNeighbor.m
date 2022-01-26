function out = isNeighbor(f, a, b, location)
%ISNEIGHBOR   Is node B a neighbor of node A?

switch lower(location)
    case {'left', 'right', 'down', 'up'}
        out = isNeighborSide(f, a, b, location);
    case {'leftdown', 'rightdown', 'leftup', 'rightup'}
        out = isNeighborCorner(f, a, b, location);
end

end

function out = isNeighborSide(f, a, b, side)

kLeft = 0;
kRight = 1;

% In the same plane?
switch lower(side)
    case 'left'
        dir = kLeft;
        dim = 1;
    case 'right'
        dir = kRight;
        dim = 1;
    case 'down'
        dir = kLeft;
        dim = 3;
    case 'up'
        dir = kRight;
        dim = 3;
end

if ( dir == kLeft )
    planeA = f.domain(dim,   a);
    planeB = f.domain(dim+1, b);
else
    planeA = f.domain(dim+1, a);
    planeB = f.domain(dim,   b);
end

if (planeA ~= planeB)
    out = false;
    return
end

if ( dim == 1 )
    i = 3;
else
    i = 1;
end

if ( ~(f.domain(i+1, b) > f.domain(i, a) && f.domain(i, b) < f.domain(i+1, a)) )
    out = false;
else
    out = true;
end

end

function out = isNeighborCorner(f, a, b, corner)

switch lower(corner)
    case 'leftdown'
        out = all(f.domain([1 3], a) == f.domain([2 4], b));
    case 'rightdown'
        out = all(f.domain([2 3], a) == f.domain([1 4], b));
    case 'leftup'
        out = all(f.domain([1 4], a) == f.domain([2 3], b));
    case 'rightup'
        out = all(f.domain([2 4], a) == f.domain([1 3], b));
end

end
