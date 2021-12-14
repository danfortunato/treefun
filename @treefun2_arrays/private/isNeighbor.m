function out = isNeighbor(f, a, b, side)
%ISNEIGHBOR   Is node B a neighbor of node A?

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
    planeA = f.domain(a, dim);
    planeB = f.domain(b, dim+1);
else
    planeA = f.domain(a, dim+1);
    planeB = f.domain(b, dim);
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

if ( ~(f.domain(b, i+1) > f.domain(a, i) && f.domain(b, i) < f.domain(a, i+1)) )
    out = false;
else
    out = true;
end

end
