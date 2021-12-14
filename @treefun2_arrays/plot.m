function varargout = plot(f, varargin)
%PLOT   Plot a TREEFUN2.
%   PLOT(F) gives a 2D color plot of the TREEFUN2 F and shows the tree on
%   which F is defined.
%
%   See also SURF, MESH.

holdState = ishold();
nplotpts = 100;

% Plot the function
hold on
boxes = leaves(f);
for k = 1:length(boxes)
    id = boxes(k);
    [x, y] = meshgrid(linspace(f.domain(1, id), f.domain(2, id), nplotpts), ...
                      linspace(f.domain(3, id), f.domain(4, id), nplotpts));
    u = coeffs2plotvals(f.coeffs{id});
    hk = surface(x, y, 0*u, u, 'EdgeAlpha', 1, varargin{:});
    if ( nargout > 0 )
        h(k) = hk; %#ok<AGROW>
    end
end
shading interp
view(0, 90)

% Plot the boxes
for k = 1:length(boxes)
    id = boxes(k);
    vertices = f.domain([1 3; 2 3; 2 4; 1 4], id);
    line('XData', vertices([1:end 1], 1), ...
         'YData', vertices([1:end 1], 2), ...
         'LineWidth', 1)
end

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
