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
ids = leaves(f);
for k = 1:length(ids)
    id = ids(k);
    [x, y] = meshgrid(linspace(f.domain(1,id), f.domain(2,id), nplotpts), ...
                      linspace(f.domain(3,id), f.domain(4,id), nplotpts));
    u = coeffs2plotvals(f.coeffs{id});
    hk = surf(x, y, u, varargin{:});
    if ( nargout > 0 )
        h(k) = hk; %#ok<AGROW>
    end
end
shading interp
view(3)

% Plot the boxes
xdata = [f.domain([1 2 2 1 1], ids) ; nan(1, length(ids))];
ydata = [f.domain([3 3 4 4 3], ids) ; nan(1, length(ids))];
line('XData', xdata(:), 'YData', ydata(:), 'LineWidth', 1)

xlim(f.domain(1:2))
ylim(f.domain(3:4))

if ( ~holdState )
    hold off
end

if ( nargout > 0 )
    varargout = {h};
end

end
