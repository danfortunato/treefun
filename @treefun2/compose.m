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
    leaf = leaves(f);
    for k = 1:length(leaf)
        ff = coeffs2vals(leaf(k).coeffs);
        ff(abs(ff) < 1e-20) = 1e-20;
        f.boxes(leaf(k).id).coeffs = vals2coeffs(op(ff));
    end
elseif ( isa(f, 'treefun2') && isa(g, 'treefun2') )
    % F and G are both TREEFUN2s:
    leaf_f = leaves(f);
    leaf_g = leaves(g);
    for k = 1:length(leaff)
        ff = coeffs2vals(leaf_f(k).coeffs);
        gg = coeffs2vals(leaf_g(k).coeffs);
        f.boxes(leaf_f(k).id).coeffs = vals2coeffs(op(ff,gg));
    end
elseif ( isa(f, 'treefun2') )
    % G is not a TREEFUN2:
    leaf = leaves(f);
    for k = 1:length(leaf)
        ff = coeffs2vals(leaf(k).coeffs);
        f.boxes(leaf(k).id).coeffs = vals2coeffs(op(ff,g));
    end
elseif ( isa(g, 'treefun2') )
    % F is not a TREEFUN2:
    leaf = leaves(g);
    for k = 1:length(leaf)
        gg = coeffs2vals(leaf(k).coeffs);
        g.boxes(leaf(k).id).coeffs = vals2coeffs(op(f,gg));
    end
    f = g;
end

end
