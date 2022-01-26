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
    ids = leaves(h);
    for id = ids(:).'
        h.coeffs{id}(1,1) = h.coeffs{id}(1,1) + g;
    end
elseif ( isa(f, 'treefun2') && isa(g, 'treefun2') )
    % TODO: Assume on the same grid for now.
    h = f;
    ids = leaves(h);
    for id = ids(:).'
        h.coeffs{id} = f.coeffs{id} + g.coeffs{id};
    end
end

end
