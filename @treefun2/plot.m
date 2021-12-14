function varargout = plot(f, varargin)
%PLOT   Plot a TREEFUN2.
%   PLOT(F) gives a 2D color plot of the TREEFUN2 F and shows the tree on
%   which F is defined.
%
%   See also SURF, MESH.

holdState = ishold();
nplotpts = 400;

% Plot the function
hold on
boxes = leaves(f);
% for k = 1:length(boxes)
%     box = boxes(k);
%     [x, y] = meshgrid(linspace(box.domain(1), box.domain(2), nplotpts), ...
%                       linspace(box.domain(3), box.domain(4), nplotpts));
%     u = coeffs2plotvals(box.coeffs);
%     hk = surface(x, y, 0*u, u, 'EdgeAlpha', 1, varargin{:});
%     if ( nargout > 0 )
%         h(k) = hk; %#ok<AGROW>
%     end
% end

[x, y] = meshgrid(linspace(f.boxes(1).domain(1), f.boxes(1).domain(2), nplotpts), ...
                  linspace(f.boxes(1).domain(3), f.boxes(1).domain(4), nplotpts));
u = feval(f, x, y);
h = surface(x, y, 0*u, u, 'EdgeAlpha', 1, varargin{:});
shading interp
view(2)

% Plot the boxes
for k = 1:length(boxes)
    vertices = boxes(k).domain([1 3; 2 3; 2 4; 1 4]);
    line('XData', vertices([1:end 1], 1), ...
         'YData', vertices([1:end 1], 2), ...
         'LineWidth', 1)
end

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
