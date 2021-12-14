function out = isNeighbor(a, b, location)
%ISNEIGHBOR   Is node B a neighbor of node A?

switch lower(location)
    case {'left', 'right', 'down', 'up'}
        out = isNeighborSide(a, b, location);
    case {'leftdown', 'rightdown', 'leftup', 'rightup'}
        out = isNeighborCorner(a, b, location);
end

end

function out = isNeighborSide(a, b, side)

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
    planeA = a.domain(dim);
    planeB = b.domain(dim+1);
else
    planeA = a.domain(dim+1);
    planeB = b.domain(dim);
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

if ( ~(b.domain(i+1) > a.domain(i) && b.domain(i) < a.domain(i+1)) )
    out = false;
else
    out = true;
end

end

function out = isNeighborCorner(a, b, corner)

switch lower(corner)
    case 'leftdown'
        out = all(a.domain([1 3]) == b.domain([2 4]));
    case 'rightdown'
        out = all(a.domain([2 3]) == b.domain([1 4]));
    case 'leftup'
        out = all(a.domain([1 4]) == b.domain([2 3]));
    case 'rightup'
        out = all(a.domain([2 4]) == b.domain([1 3]));
end

end
