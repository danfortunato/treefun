function varargout = pcolor(f, varargin)
%PCOLOR   Pseudocolor (checkerboard) plot of an ULTRASEM.SOL.
%   PCOLOR(SOL) is a pseudocolor or "checkerboard" plot of the solution
%   SOL.
%
%   H = PCOLOR(...) returns a handle to a SURFACE object.
%
%   See also PLOT, SURF.

holdState = ishold();

% Plot the function
hold on
boxes = leaves(f);
for k = 1:length(boxes)
    dom = boxes(k).domain;
    [x, y] = meshgrid(linspace(dom(1), dom(2), 200), ...
                      linspace(dom(3), dom(4), 200));
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
