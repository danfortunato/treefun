function varargout = surf(f, varargin)
%SURF   Surface plot of a TREEFUN2.
%   SURF(F) plots a colored parametric surface whose height and color is
%   defined by the values of F.
%
%   SURF(..., 'PropertyName', PropertyValue, ...) sets the value of the
%   specified surface property. Multiple property values can be set with a
%   single statement.
%
%   H = SURF(...) returns a handle to a surface plot object.
%
%   See also PLOT, MESH.

holdState = ishold();
nplotpts = 100;

% Plot the function
hold on
boxes = leaves(f);
for k = 1:length(boxes)
    box = boxes(k);
    [x, y] = meshgrid(linspace(box.domain(1), box.domain(2), nplotpts), ...
                      linspace(box.domain(3), box.domain(4), nplotpts));
    u = coeffs2plotvals(box.coeffs);
    hk = surf(x, y, u, varargin{:});
    if ( nargout > 0 )
        h(k) = hk; %#ok<AGROW>
    end
end
shading interp
view(3)

% Plot the boxes
for k = 1:length(boxes)
    vertices = boxes(k).domain([1 3; 2 3; 2 4; 1 4]);
    line('XData', vertices([1:end 1], 1), ...
         'YData', vertices([1:end 1], 2), ...
         'LineWidth', 1)
end

xlim(f.domain(1:2))
ylim(f.domain(3:4))

if ( ~holdState )
    hold off
end

if ( nargout > 0 )
    varargout = {h};
end

end
