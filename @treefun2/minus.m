function h = minus(f, g)
%-   Minus for TREEFUN2.
%   F - G subtracts G from F, where F and G are TREEFUN2 objects. F and G
%   must have the same domains and discretization sizes. F and G may also
%   be scalars.
%
%   See also PLUS.

h = plus( f, (-g) );

end
