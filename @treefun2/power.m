function f = power(f, g)
%.^   Pointwise power of a TREEFUN2.
%   F.^G returns a TREEFUN2 F to the scalar power G, a scalar F to the TREEFUN2
%   power G, or a TREEFUN2 F to the TREEFUN2 power G.
%
%   See also COMPOSE.

if ( isempty(f) )
    return
elseif ( isempty(g) )
    f = g;
    return
else
    f = compose(@power, f, g);
end

end
