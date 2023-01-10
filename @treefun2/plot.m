function varargout = plot(f, varargin)
%PLOT   Plot a TREEFUN2.
%   PLOT(F) gives a 2D color plot of the TREEFUN2 F and shows the tree on
%   which F is defined.
%
%   See also SURF, MESH.

holdState = ishold();
nplotpts = 800;

% Plot the function
hold on
ids = leaves(f);
% for k = 1:length(ids)
%     id = ids(k);
%     [x, y] = meshgrid(linspace(f.domain(1,id), f.domain(2,id), nplotpts), ...
%                       linspace(f.domain(3,id), f.domain(4,id), nplotpts));
%     u = coeffs2plotvals(f.coeffs{id});
%     hk = surface(x, y, 0*u, u, 'EdgeAlpha', 1, varargin{:});
%     if ( nargout > 0 )
%         h(k) = hk; %#ok<AGROW>
%     end
% end

[x, y] = meshgrid(linspace(f.domain(1,1), f.domain(2,1), nplotpts), ...
                  linspace(f.domain(3,1), f.domain(4,1), nplotpts));
u = feval(f, x, y);
h = surface(x, y, 0*u, u, 'EdgeAlpha', 1, varargin{:});
shading interp
view(2)

% Plot the boxes
xdata = [f.domain([1 2 2 1 1], ids) ; nan(1, length(ids))];
ydata = [f.domain([3 3 4 4 3], ids) ; nan(1, length(ids))];
line('XData', xdata(:), 'YData', ydata(:), 'LineWidth', 1)

axis equal
xlim(f.domain(1:2, 1))
ylim(f.domain(3:4, 1))

if ( ~holdState )
    hold off
end

if ( nargout > 0 )
    varargout = {h};
end

end
