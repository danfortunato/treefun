classdef treefun3  %#ok<*PROP,*PROPLC>
%TREEFUN3   Piecewise polynomial on an adaptive quadtree.
%   TREEFUN3(F) constructs a TREEFUN3 object representing the function F on
%   the domain [-1, 1] x [-1, 1] x [-1, 1]. F may be a function handle.
%   A TREEFUN3 is constructed by recursively subdividing the domain until
%   each piece is well approximated by a trivariate polynomial of degree
%   (N-1) x (N-1). The default is N = 16. 
%
%   TREEFUN3(F, N) uses piecewise polynomials of degree (N-1) x (N-1).
%
%   TREEFUN3(F, [A B C D E F]) specifies a domain [A, B] x [C, D] x [E, F]
%   on which the TREEFUN3 is defined.
%
%   TREEFUN3(F, [A B C D E F], N) specifies both a degree and a domain.
%
%   TREEFUN3(F, TF) creates a TREEFUN3 approximation to F using the tree
%   structure from a previously-defined TREEFUN3. No adaptive construction
%   takes place; the function F is simply initialized on the grid inherited
%   from TF.
%
%   not all preperties are used ...

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
        morton
        flatNeighbors
        leafNeighbors

    end

    methods

        function f = treefun3(varargin)

            if ( nargin < 1 )
                return
            end

            dom = [-1 1 -1 1 -1 1];
            opts = struct();
            opts.balance = false; % not yet
            opts.neighbors = false;
            
            if ( nargin == 2 )
                if ( isa(varargin{2}, 'treefun3') ) % TREEFUN3(F, TF)
                    % We were given the tree structure
                    f = varargin{2};
                    func = varargin{1};
                    if ( isnumeric(func) && isscalar(func) )
                        func = @(x,y) func + 0*x;
                    elseif ( isa(func, 'chebfun3') )
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
                        f.coeffs{id} = treefun3.vals2coeffs(vals);
                    end
                    return
                elseif ( isscalar(varargin{2}) ) % TREEFUN3(F, N)
                    f.n = varargin{2};
                else
                    dom = varargin{2};           % TREEFUN3(F, [A B C D E F])
                end
            elseif ( nargin == 3 )
                dom = varargin{2};               % TREEFUN3(F, [A B C D E F], N)
                f.n = varargin{3};
            elseif ( nargin == 4 )
                dom = varargin{2};               % TREEFUN3(F, [A B C D E F], N, OPTS)
                f.n = varargin{3};
                opts1 = varargin{4};
                if ( isfield(opts1, 'balance') )
                    opts.balance = opts1.balance;
                end
                if ( isfield(opts1, 'neighbors') )
                    opts.neighbors = opts1.neighbors;
                end
            elseif ( nargin == 9 )
                % TREEFUN3(DOMAIN, LEVEL, HEIGHT, ID, PARENT, CHILDREN,
                %   COEFFS, COL, ROW)
                [f.domain, f.level, f.height, f.id, f.parent, ...
                    f.children, f.coeffs, f.col, f.row] = deal(varargin{:});
                f.morton = cartesian2morton(f.col, f.row);
                f.n = size(f.coeffs{end}, 1);
                f = balance(f);
                [f.flatNeighbors, f.leafNeighbors] = generateNeighbors(f);
                return
            end

            func = varargin{1};
            if ( isnumeric(func) && isscalar(func) )
                func = @(x,y) func + 0*x;
            elseif ( isa(func, 'chebfun3') )
                func = @(x,y) feval(func, x, y);
            end

            f.domain(:,1)   = dom(:);
            f.level(1)      = 0;
            f.height(1)     = 0;
            f.id(1)         = 1;
            f.parent(1)     = 0;
            f.children(:,1) = zeros(8, 1);
            f.coeffs{1}     = [];
            f.col           = uint64(0);
            f.row           = uint64(0);
            
            % f = buildDepthFirst(f, func);
            f = buildBreadthFirst(f, func);
            % f.morton = cartesian2morton(f.col, f.row);

            % Now do level restriction
            opts.balance = false;
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

    end

    methods ( Static )

        % u = poisson(f, isource);
        coeffs = vals2coeffs(vals);
        vals = coeffs2vals(coeffs);
        % vals = coeffs2refvals(coeffs);
        % refvals = chebvals2refvals(chebvals);

    end
    

end

function [resolved, coeffs] = isResolved(f, dom, n)

persistent xx0 yy0 zz0 xxx0 yyy0 zzz0 nstored

tol = 1e-12;
nalias = n;
nrefpts = 2*n; % Sample at equispaced points to test error

if ( isempty(xx0) || isempty(xxx0) || n ~= nstored )
    nstored = n;
    x0 = (1-cos(pi*(2*(1:nalias)'-1)/(2*nalias)))/2;
    [xx0, yy0, zz0] = ndgrid(x0);
    [xxx0, yyy0, zzz0] = ndgrid(linspace(0, 1, nrefpts));
end
sclx = diff(dom(1:2));
scly = diff(dom(3:4));
sclz = diff(dom(5:6));
xx = sclx*xx0 + dom(1); xxx = sclx*xxx0 + dom(1);
yy = scly*yy0 + dom(3); yyy = scly*yyy0 + dom(3);
zz = sclz*zz0 + dom(5); zzz = sclz*zzz0 + dom(5);

vals = f(xx,yy,zz);
coeffs = treefun3.vals2coeffs(vals);
coeffs = coeffs(1:n,1:n);
Ex = sum(abs(coeffs(end-1:end,:)), 'all') / (2*n);
Ey = sum(abs(coeffs(:,end-1:end)), 'all') / (2*n);
err_cfs = (Ex + Ey) / 2;

% F = f(xxx,yyy,zzz);
% G = coeffs2refvals(coeffs);
% err_vals = max(abs(F(:) - G(:)));

err = err_cfs;
%err = min(err_cfs, err_vals);
h = sclx;
eta = 0;

vmax = max(abs(vals(:)));
resolved = ( err * h^eta < tol * max(vmax, 1) );

end
