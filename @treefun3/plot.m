function varargout = plot(f, varargin)
%PLOT   Plot a TREEFUN3.
%   PLOT(F) gives a 3D plot of the TREEFUN3 F and shows the tree on
%   which F is defined.
%
%   See also SURF, MESH.

% func = varargin{1};

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

[x, y, z] = meshgrid(linspace(f.domain(1), f.domain(2), nplotpts), ...
                     linspace(f.domain(3), f.domain(4), nplotpts), ...
                     linspace(f.domain(5), f.domain(6), nplotpts));
% v = func( x, y, z);

% Plot the boxes
xdata = [f.domain([1 2 2 1 1 1 2 2 1 1], ids) ; nan(1, length(ids)); ... % bottom & top
         f.domain([2 2 2 2 1 1], ids) ; nan(1, length(ids));]; % the rest to complete the box
ydata = [f.domain([3 3 4 4 3 3 3 4 4 3], ids) ; nan(1, length(ids)); ...
         f.domain([3 3 4 4 4 4], ids) ; nan(1, length(ids));];
zdata = [f.domain([5 5 5 5 5 6 6 6 6 6], ids) ; nan(1, length(ids)); ...
         f.domain([5 6 6 5 5 6], ids) ; nan(1, length(ids));];
line('XData', xdata(:), 'YData', ydata(:), 'ZData', zdata(:), 'LineWidth', 1)

axis equal
xlim(f.domain(1:2, 1))
ylim(f.domain(3:4, 1))
zlim(f.domain(5:6, 1))

if ( ~holdState )
    hold off
end

if ( nargout > 0 )
    varargout = {h};
end

end
