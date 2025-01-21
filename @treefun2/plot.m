function varargout = plot(f, dom, varargin)
%PLOT   Plot a TREEFUN2.
%   PLOT(F) gives a 2D color plot of the TREEFUN2 F and shows the tree on
%   which F is defined.
%
%   See also SURF, MESH.

if ( isempty(f) )
    return
end

doLabel = false;
for k = 1:length(varargin)
    if ( isstring(varargin{k}) && lower(varargin{k}) == "label" )
        doLabel = varargin{k+1};
        varargin(k:k+1) = [];
        break
    end
end

holdState = ishold();
nplotpts = 800;
root = f.root;

if ( nargin < 2 )
    dom = f.domain(:,root);
end

% Plot the function
hold on
ids = leaves(f);
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

% Plot the boxes
xdata = [f.domain([1 2 2 1 1], ids) ; nan(1, length(ids))];
ydata = [f.domain([3 3 4 4 3], ids) ; nan(1, length(ids))];
line('XData', xdata(:), 'YData', ydata(:), 'LineWidth', 1)

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
