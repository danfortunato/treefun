function f = ldivide(f, g)
%.\   Pointwise left divide for TREEFUN2.
%   F.\G divides G by F, where F and G may be TREEFUN2 objects or scalars.
%
%   See also RDIVIDE, COMPOSE.

f = rdivide(g, f);

end
