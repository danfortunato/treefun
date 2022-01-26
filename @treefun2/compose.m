function f = compose(op, f, g)
%COMPOSE   Compose a function handle with a TREEFUN2.
%   COMPOSE(OP, F) returns a TREEFUN2 representing OP(F), where OP is a
%   function handle and F is a TREEFUN2. The function handle OP is applied
%   to F in value space.
%
%   COMPOSE(OP, F, G) returns a TREEFUN2 representing OP(F, G), where OP
%   is a function handle and at least one of F and G is a TREEFUN2. The
%   function handle OP is applied to F and G in value space.

if ( nargin == 2 )
    ids = leaves(f);
    for id = ids(:).'
        ff = treefun2.coeffs2vals(f.coeffs{id});
        ff(abs(ff) < 1e-20) = 1e-20;
        f.coeffs{id} = treefun2.vals2coeffs(op(ff));
    end
elseif ( isa(f, 'treefun2') && isa(g, 'treefun2') )
    % F and G are both TREEFUN2s:
    ids_f = leaves(f);
    ids_g = leaves(g);
    for k = 1:length(ids_f)
        ff = treefun2.coeffs2vals(f.coeffs{ids_f(k)});
        gg = treefun2.coeffs2vals(g.coeffs{ids_g(k)});
        f.coeffs{ids_f(k)} = treefun2.vals2coeffs(op(ff,gg));
    end
elseif ( isa(f, 'treefun2') )
    % G is not a TREEFUN2:
    ids = leaves(f);
    for id = ids(:).'
        ff = treefun2.coeffs2vals(f.coeffs{id});
        f.coeffs{id} = treefun2.vals2coeffs(op(ff,g));
    end
elseif ( isa(g, 'treefun2') )
    % F is not a TREEFUN2:
    ids = leaves(g);
    for id = ids(:).'
        gg = treefun2.coeffs2vals(g.coeffs{id});
        g.coeffs{id} = treefun2.vals2coeffs(op(f,gg));
    end
    f = g;
end

end
