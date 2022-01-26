function h = plus(f, g)
%+   Plus for TREEFUN2.
%   F + G adds the TREEFUN2 objects F and G. F and G must have the same
%   domains and discretization sizes. F and G may also be scalars.
%
%   See also MINUS.

if ( isnumeric( f ) )
    h = plus(g, f);
    return
elseif ( isnumeric( g ) )
    h = f;
    boxes = leaves(h);
    for k = 1:length(boxes)
        id = boxes(k).id;
        h.boxes(id).coeffs(1,1) = h.boxes(id).coeffs(1,1) + g;
    end
elseif ( isa(f, 'treefun2') && isa(g, 'treefun2') )
    % TODO: Assume on the same grid for now.
    h = f;
    boxes = leaves(h);
    for k = 1:length(boxes)
        id = boxes(k).id;
        h.boxes(id).coeffs = f.boxes(id).coeffs + g.boxes(id).coeffs;
    end
end

end
