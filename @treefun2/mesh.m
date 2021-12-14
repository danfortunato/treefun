function varargout = mesh(f, varargin)
%MESH   Plot the mesh of an TREEFUN2.
%   MESH(F) plots the tensor product Chebyshev grids for the patches of F
%   colored according to the values of F.
%
%   See also PLOT, SURF.

holdState = ishold();

hold on
[x, y] = leafpts(f);
vals = leafvals(f);
for k = 1:length(f)
    u = vals{k};
    if ( ~isreal(u) )
        u = abs(u);
    end
    hk = mesh(x{k}, y{k}, u, varargin{:});
    if ( nargout > 0 )
        h(k) = hk; %#ok<AGROW>
    end
    hold on
end
view(2)
axis equal
xlim(f.domain(1:2))
ylim(f.domain(3:4))

if ( ~holdState )
    hold off
end

if ( nargout > 0 )
    varargout = {h};
end

end
