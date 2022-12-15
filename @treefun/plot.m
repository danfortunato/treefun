function varargout = plot(f, varargin)
%PLOT   Plot a TREEFUN.
%   PLOT(F) plots the TREEFUN F and shows the tree on which F is defined.

holdState = ishold();
nplotpts = 100;
ms = 10;
lw = 1;

% Plot the function
hold on
boxes = leaves(f);
for k = 1:length(boxes)
    box = boxes(k);
    x = linspace(box.domain(1), box.domain(2), nplotpts);
    y = coeffs2plotvals(box.coeffs);
    if ( isreal(y) )
        hk = plot(x, y, 'LineWidth', lw, varargin{:});
    else
        hk = plot(y, 'LineWidth', lw, varargin{:});
    end
    if ( nargout > 0 )
        h(k) = hk; %#ok<AGROW>
    end
    if ( isreal(y) )
        plot(boxes(k).domain(1), y(1), 'k.', 'MarkerSize', ms)
    else
        plot(y(1), 'k.', 'MarkerSize', ms)
    end
end
if ( isreal(y) )
    plot(boxes(end).domain(2), y(end), 'k.', 'MarkerSize', ms)
else
    plot(y(end), 'k.', 'MarkerSize', ms)
end

if ( isreal(y) )
    xlim(f.domain)
end

if ( ~holdState )
    hold off
end

if ( nargout > 0 )
    varargout = {h};
end

end
