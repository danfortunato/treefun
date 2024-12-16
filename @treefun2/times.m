function f = times(f, g)
%.*   Pointwise multiplication for TREEFUN2.
%   F.*G multiplies F and G pointwise, where F and G may be TREEFUN2 objects or
%   scalars. F and G must have the same tree structure.
%
%   See also MTIMES, COMPOSE.

if ( ~isa(f, 'treefun2') )
    % Ensure F is the TREEFUN2:
    f = times(g, f);
    return
elseif ( isa(g, 'treefun2' ) )
    % Multiply two TREEFUN2s:
    f = compose(@times, f, g);
elseif ( isnumeric(g) && isscalar(g) )
    % Multiply TREEFUN2 F by scalar G:
    f = compose(@(x) g*x, f);
else
    error('TREEFUN2:times:invalid', ...
        'F and G must be scalars or treefun2 objects.')
end

end
