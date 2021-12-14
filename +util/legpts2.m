function [xx, yy] = legpts2(nx, ny, dom)
%LEGPTS2   Legendre tensor product grid.
%   [XX, YY] = LEGPTS2(N) constructs an N x N grid of tensor-product
%   Legendre points on [-1,1]^2.
%
%   [XX, YY] = LEGPTS2(NX, NY) constructs an NX x NY grid of
%   tensor-product Legendre points on [-1,1]^2.
%
%   [XX, YY] = LEGPTS2(NX, NY, DOM) constructs an NX x NY grid of
%   tensor-product Legendre points on the rectangle [a,b] x [c,d], where
%   DOM = [a b c d].
%
%   See also LEGPTS.

if ( nargin == 1 )
    ny = nx;
end

% Default domain:
if ( nargin < 3 )
    dom = [-1 1 -1 1];
end

x = legpts(nx, dom(1:2));
y = legpts(ny, dom(3:4));
[xx, yy] = meshgrid(x, y);

end
