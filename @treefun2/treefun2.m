classdef treefun2  %#ok<*PROP,*PROPLC>

    properties

        n = 16
        domain
        level
        height
        id
        parent
        children
        coeffs
        col
        row
        flatNeighbors
        leafNeighbors

    end

    methods

        function f = treefun2(varargin)

            if ( nargin < 1 )
                return
            end

            dom = [-1 1 -1 1];
            opts = struct();
            opts.balance = true;
            opts.neighbors = true;
            
            if ( nargin == 2 )
                if ( isa(varargin{2}, 'treefun2') ) % TREEFUN2(F, TF)
                    % We were given the tree structure
                    f = varargin{2};
                    func = varargin{1};
                    if ( isnumeric(func) && isscalar(func) )
                        func = @(x,y) func + 0*x;
                    elseif ( isa(func, 'chebfun2') )
                        func = @(x,y) feval(func, x, y);
                    end
                    % We just need to fill in the leaf coefficients
                    L = leaves(f);
                    [xx0, yy0] = chebpts2(f.n, f.n, [0 1 0 1]);
                    for id = L(:).'
                        dom = f.domain(:,id);
                        sclx = diff(dom(1:2));
                        scly = diff(dom(3:4));
                        xx = sclx*xx0 + dom(1);
                        yy = scly*yy0 + dom(3);
                        vals = func(xx,yy);
                        f.coeffs{id} = treefun2.vals2coeffs(vals);
                    end
                    return
                elseif ( isscalar(varargin{2}) ) % TREEFUN2(F, N)
                    f.n = varargin{2};
                else
                    dom = varargin{2};           % TREEFUN2(F, [A B C D])
                end
            elseif ( nargin == 3 )
                dom = varargin{2};               % TREEFUN2(F, [A B C D], N)
                f.n = varargin{3};
            elseif ( nargin == 4 )
                dom = varargin{2};               % TREEFUN2(F, [A B C D], N, OPTS)
                f.n = varargin{3};
                opts1 = varargin{4};
                if ( isfield(opts1, 'balance') )
                    opts.balance = opts1.balance;
                end
                if ( isfield(opts, 'neighbors') )
                    opts.neighbors = opts1.neighbors;
                end
            elseif ( nargin == 9 )
                % TREEFUN2(DOMAIN, LEVEL, HEIGHT, ID, PARENT, CHILDREN,
                %   COEFFS, COL, ROW)
                [f.domain, f.level, f.height, f.id, f.parent, ...
                    f.children, f.coeffs, f.col, f.row] = deal(varargin{:});
                f.n = size(f.coeffs{end}, 1);
                f = balance(f);
                [f.flatNeighbors, f.leafNeighbors] = generateNeighbors(f);
                return
            end

            func = varargin{1};
            if ( isnumeric(func) && isscalar(func) )
                func = @(x,y) func + 0*x;
            elseif ( isa(func, 'chebfun2') )
                func = @(x,y) feval(func, x, y);
            end

            f.domain(:,1)   = dom(:);
            f.level(1)      = 0;
            f.height(1)     = 0;
            f.id(1)         = 1;
            f.parent(1)     = 0;
            f.children(:,1) = zeros(4, 1);
            f.coeffs{1}     = [];
            f.col(1)        = 1;
            f.row(1)        = 1;
            
            % f = buildDepthFirst(f, func);
            f = buildBreadthFirst(f, func);

            % Now do level restriction
            if ( opts.balance )
                f = balance(f);
            else
                % Do a cumulative sum in reverse to correct the heights
                for k = length(f.id):-1:1
                    if ( ~isLeaf(f, k) )
                        f.height(k) = 1 + max(f.height(f.children(:,k)));
                    end
                end
            end

            if ( opts.neighbors )
                [f.flatNeighbors, f.leafNeighbors] = generateNeighbors(f);
            end

        end

    end

    methods ( Access = private )

        f = balance(f);
        f = refineBox(f, id);

        function f = buildBreadthFirst(f, func)

            % Note: the length changes at each iteration here
            id = 1;
            while ( id <= length(f.id) )
                [resolved, coeffs] = isResolved(func, f.domain(:,id), f.n);
                if ( resolved )
                    f.coeffs{id} = coeffs;
                    f.height(id) = 0;
                else
                    % Split into four child boxes
                    f = refineBox(f, id);
                    f.height(id) = 1;
                end
                id = id + 1;
            end

            % Do a cumulative sum in reverse to correct the heights
            for k = length(f.id):-1:1
                if ( ~isLeaf(f, k) )
                    %f.height(k) = f.height(k) + max(f.height(f.children(:,k)));
                    f.height(k) = 1 + max(f.height(f.children(:,k)));
                end
            end

        end
        
        function f = buildDepthFirst(f, func, id, level)

            if ( nargin == 2 )
                id = 1;
                level = 0;
            end

            f.level(id) = level;
            f.height(id) = 0;

            [resolved, coeffs] = isResolved(func, f.domain(:,id), f.n);

            if ( resolved )
                f.coeffs{id} = coeffs;
            else
                % Split into four child boxes
                f = refineBox(f, id);

                % Recurse
                f = buildDepthFirst(f, func, f.children(1,id), level+1);
                f = buildDepthFirst(f, func, f.children(2,id), level+1);
                f = buildDepthFirst(f, func, f.children(3,id), level+1);
                f = buildDepthFirst(f, func, f.children(4,id), level+1);

                % Set height
                f.height(id) = 1 + max(f.height(f.children(:,id)));
            end

        end

    end

    methods ( Static )

        u = poisson(f);
        coeffs = vals2coeffs(vals);
        vals = coeffs2vals(coeffs);
        vals = coeffs2refvals(coeffs);
        refvals = chebvals2refvals(chebvals);

    end

end

function [resolved, coeffs] = isResolved(f, dom, n)

persistent xx0 yy0 xxx0 yyy0 nstored

%tol = 1e-9;
tol = 1e-12;
nrefpts = 2*n; % Sample at equispaced points to test error

if ( isempty(xx0) || isempty(xxx0) || n ~= nstored )
    nstored = n;
    [xx0, yy0] = chebpts2(n, n, [0 1 0 1]);
    [xxx0, yyy0] = meshgrid(linspace(0, 1, nrefpts));
end
sclx = diff(dom(1:2));
scly = diff(dom(3:4));
xx = sclx*xx0 + dom(1); xxx = sclx*xxx0 + dom(1);
yy = scly*yy0 + dom(3); yyy = scly*yyy0 + dom(3);

vals = f(xx,yy);
coeffs = treefun2.vals2coeffs(vals);
Ex = sum(abs(coeffs(end-1:end,:)), 'all') / (2*n);
Ey = sum(abs(coeffs(:,end-1:end)), 'all') / (2*n);
err_cfs = (Ex + Ey) / 2;

% F = f(xxx,yyy);
% G = coeffs2refvals(coeffs);
% err_vals = norm(F(:) - G(:), inf);

err = err_cfs;

vmax = max(abs(vals(:)));
resolved = ( err < tol * max(vmax, 1) );

end
