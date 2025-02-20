function varargout = plot(f, varargin)
%PLOT   Plot a TREEFUN2.
%   PLOT(F) gives a 2D color plot of the TREEFUN2 F and shows the tree on
%   which F is defined.
%
%   See also SURF, MESH.

if ( isempty(f) )
    return
end

holdState = ishold();
root = f.root;
dom = f.domain(:,root);
nplotpts = 400;
doLabel  = false;
doBoxes  = true;
doValues = true;
filter = [];
for k = 1:length(varargin)
    if ( isstring(varargin{k}) && lower(varargin{k}) == "axis" )
        dom = varargin{k+1};
        filter = [filter k:k+1];
        continue
    end
    if ( isstring(varargin{k}) && lower(varargin{k}) == "nplotpts" )
        nplotpts = varargin{k+1};
        filter = [filter k:k+1];
        continue
    end
    if ( isstring(varargin{k}) && lower(varargin{k}) == "label" )
        doLabel = varargin{k+1};
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
ids = leaves(f);

% Plot the function
if ( doValues )
    % nplotpts = 200;
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
    [x, y] = meshgrid(linspace(dom(1), dom(2), nplotpts), ...
                      linspace(dom(3), dom(4), nplotpts));
    u = feval(f, x, y);
    h = surface(x, y, 0*u, u, 'EdgeAlpha', 1, varargin{:});
    shading interp
    view(2)
end

% Plot the boxes
if ( doBoxes )
    xdata = [f.domain([1 2 2 1 1], ids) ; nan(1, length(ids))];
    ydata = [f.domain([3 3 4 4 3], ids) ; nan(1, length(ids))];
    line('XData', xdata(:), 'YData', ydata(:), 'LineWidth', 1)
end

axis equal
xlim(dom(1:2))
ylim(dom(3:4))

if ( doLabel )
    cx = sum(f.domain(1:2,ids)) / 2;
    cy = sum(f.domain(3:4,ids)) / 2;
    text(cx, cy, int2str(ids(:)), 'HorizontalAlignment', 'center')
end

if ( ~holdState )
    hold off
end

if ( nargout > 0 )
    varargout = {h};
end

end
