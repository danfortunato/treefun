function varargout = pcolor(f, varargin)
%PCOLOR   Pseudocolor (checkerboard) plot of a TREEFUN2.
%   PCOLOR(F) is a pseudocolor or "checkerboard" plot of the TREEFUN2 F.
%
%   H = PCOLOR(...) returns a handle to a SURFACE object.
%
%   See also PLOT, SURF.

holdState = ishold();
nplotpts = 100;

% Plot the function
hold on
boxes = leaves(f);
for k = 1:length(boxes)
    dom = boxes(k).domain;
    [x, y] = meshgrid(linspace(dom(1), dom(2), nplotpts), ...
                      linspace(dom(3), dom(4), nplotpts));
    u = coeffs2plotvals(boxes(k).coeffs);
    hk = surface(x, y, 0*u, u, 'EdgeAlpha', 1, varargin{:});
    if ( nargout > 0 )
        h(k) = hk; %#ok<AGROW>
    end
end
shading interp
view(0, 90)

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
