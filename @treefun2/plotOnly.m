function varargout = plotOnly(f, ids, varargin)
%PLOT   Plot a TREEFUN2.
%   PLOT(F) gives a 2D color plot of the TREEFUN2 F and shows the tree on
%   which F is defined.
%
%   See also SURF, MESH.

holdState = ishold();
nplotpts = 100;
doBoxes  = true;
doValues = true;
filter = [];
for k = 1:length(varargin)
    if ( isstring(varargin{k}) && lower(varargin{k}) == "nplotpts" )
        nplotpts = varargin{k+1};
        filter = [filter k:k+1];
        continue
    end
    if ( isstring(varargin{k}) && lower(varargin{k}) == "boxes" )
        doBoxes = varargin{k+1};
        filter = [filter k:k+1];
        continue
    end
    if ( isstring(varargin{k}) && lower(varargin{k}) == "values" )
        doValues = varargin{k+1};
        filter = [filter k:k+1];
        continue
    end
end
varargin(filter) = [];

hold on

% Plot the function
if ( doValues )
    for k = 1:length(ids)
        id = ids(k);
        [x, y] = meshgrid(linspace(f.domain(1,id), f.domain(2,id), nplotpts), ...
                          linspace(f.domain(3,id), f.domain(4,id), nplotpts));
        u = coeffs2plotvals(f.coeffs{id});
        %u = feval(f, x, y);
        hk = surface(x, y, 0*u, u, 'EdgeAlpha', 0, varargin{:});
        if ( nargout > 0 )
            h(k) = hk; %#ok<AGROW>
        end
    end
    shading interp
    view(2)
end

% Plot the boxes
if ( doBoxes )
    xdata = [f.domain([1 2 2 1 1], ids) ; nan(1, length(ids))];
    ydata = [f.domain([3 3 4 4 3], ids) ; nan(1, length(ids))];
    line('XData', xdata(:), 'YData', ydata(:), 'LineWidth', 1, varargin{:})
end

axis equal
xlim(f.domain(1:2, f.root))
ylim(f.domain(3:4, f.root))

if ( ~holdState )
    hold off
end

if ( nargout > 0 )
    varargout = {h};
end

end
