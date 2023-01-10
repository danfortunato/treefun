classdef treefun
%TREEFUN   Piecewise polynomial on an adaptive binary tree.
%   TREEFUN(F) constructs a TREEFUN object representing the function F on
%   the interval [-1, 1]. F may be a function handle, scalar, or chebfun
%   object. A TREEFUN is constructed by recursively subdividing the
%   domain until each piece is well approximated by a polynomial of degree
%   N-1. The default is N = 16.
%
%   TREEFUN(F, N) uses piecewise polynomials of degree N-1.
%
%   TREEFUN(F, [A B]) specifies a domain [A, B] on which the TREEFUN is
%   defined.
%
%   TREEFUN(F, [A B], N) specifies both a degree and a domain.
%
%   TREEFUN(F, N, TOL), TREEFUN(F, [A B], TOL), or TREEFUN(F, [A B], N,
%   TOL) resolves each piece to the specified tolerance TOL. The default is
%   TOL = 1e-12.

    properties

        boxes
        n = 16
        domain = [-1 1]
        tol = 1e-12

    end

    methods

        function f = treefun(varargin)

            if ( nargin < 1 )
                return
            end

            func = varargin{1};
            if ( isnumeric(func) && isscalar(func) )
                func = @(x) func + 0*x;
            elseif ( isa(func, 'chebfun') )
                func = @(x) feval(func, x);
            end

            if ( nargin == 2 )
                if ( isa(varargin{2}, 'treefun') ) % TREEFUN(F, TF)
                    % We were given the tree structure
                    f = varargin{2};
                    % We just need to fill in the leaf coefficients
                    L = [leaves(f).id];
                    x0 = chebpts(f.n, [0 1]);
                    for k = 1:length(L)
                        id = L(k);
                        dom = f.boxes(id).domain;
                        sclx = diff(dom);
                        x = sclx*x0 + dom(1);
                        vals = func(x);
                        f.boxes(id).coeffs = treefun.vals2coeffs(vals);
                    end
                    return
                elseif ( isscalar(varargin{2}) )   % TREEFUN(F, N)
                    f.n = varargin{2};
                else                               % TREEFUN(F, [A B])
                    f.domain = varargin{2};
                end
            elseif ( nargin == 3 )
                if ( length(varargin{2}) == 1 )    % TREEFUN(F, N, TOL)
                    f.n = varargin{2};
                    f.tol = varargin{3};
                elseif ( varargin{3} > 1 )         % TREEFUN(F, [A B], N)
                    f.domain = varargin{2};
                    f.n = varargin{3};
                else                               % TREEFUN(F, [A B], TOL)
                    f.domain = varargin{2};
                    f.tol = varargin{3};
                end
            elseif ( nargin == 4 )                 % TREEFUN(F, [A B], N, TOL)
                f.domain = varargin{2};
                f.n = varargin{3};
                f.tol = varargin{4};
            end

            f.boxes = struct();
            f.boxes(1).domain   = f.domain;
            f.boxes(1).level    = 0;
            f.boxes(1).height   = 0;
            f.boxes(1).id       = 1;
            f.boxes(1).parent   = 0;
            f.boxes(1).children = [];
            f.boxes(1).coeffs   = [];
            f.boxes(1).col      = 1;

            f = buildBreadthFirst(f, func, @isResolvedValues);

        end

    end

    methods ( Access = private )

        function f = refineBox(f, id)

            % Split into two child boxes
            dom = f.boxes(id).domain;
            xmid = mean(dom);
            parent = f.boxes(id);

            child1 = struct();
            child1.domain   = [dom(1) xmid];
            child1.id       = length(f.boxes)+1;
            child1.parent   = id;
            child1.children = [];
            child1.level    = parent.level+1;
            child1.height   = 0;
            child1.coeffs   = [];
            child1.col      = 2*parent.col;
            f.boxes(end+1) = child1;

            child2 = struct();
            child2.domain   = [xmid dom(2)];
            child2.id       = length(f.boxes)+1;
            child2.parent   = id;
            child2.children = [];
            child2.level    = parent.level+1;
            child2.height   = 0;
            child2.coeffs   = [];
            child2.col      = 2*parent.col + 1;
            f.boxes(end+1) = child2;

            f.boxes(id).children = [child1.id, child2.id];
            f.boxes(id).height = 1;
            f.boxes(id).coeffs = [];

        end

        function f = buildBreadthFirst(f, func, isResolved)

            % Note: the length changes at each iteration here
            id = 1;
            while ( id <= length(f.boxes) )
                [resolved, coeffs] = isResolved(func, f.boxes(id).domain, f.n, f.tol);
                if ( resolved )
                    f.boxes(id).coeffs = coeffs;
                    f.boxes(id).height = 0;
                else
                    % Split into two child boxes
                    f = refineBox(f, id);
                    f.boxes(id).height = 1;
                end
                id = id + 1;
            end

            % Do a cumulative sum in reverse to correct the heights
            for k = length(f.boxes):-1:1
               if ( f.boxes(k).height ~= 0 )
                   %f.boxes(k).height = f.boxes(k).height + max([f.boxes(f.boxes(k).children).height]);
                   f.boxes(k).height = 1 + max([f.boxes(f.boxes(k).children).height]);
               end
            end

        end

    end

    methods ( Static )

        coeffs = vals2coeffs(vals);
        vals = coeffs2vals(coeffs);

    end

end

function [resolved, coeffs] = isResolvedValues(f, dom, n, tol)

nrefpts = 2*n; % Sample at equispaced points to test error

persistent xx0 xxx0 nstored
if ( isempty(xx0) || isempty(xxx0) || n ~= nstored )
    nstored = n;
    xx0 = chebpts(n, [0 1]);
    xxx0 = linspace(0, 1, nrefpts).';
end
sclx = diff(dom);
xx  = sclx*xx0  + dom(1);
xxx = sclx*xxx0 + dom(1);

vals = f(xx);
coeffs = treefun.vals2coeffs(vals);
F = f(xxx);
G = coeffs2refvals(coeffs);
err = norm(F(:) - G(:), inf);
vmax = max(abs(vals(:)));
resolved = ( err < tol * max(vmax, 1) );

end

function [resolved, coeffs] = isResolvedCoeffs(f, dom, n, tol)

persistent xx0 nstored
if ( isempty(xx0) || n ~= nstored )
    nstored = n;
    xx0 = chebpts(n, [0 1 0 1]);
end
sclx = diff(dom);
xx = sclx*xx0  + dom(1);

vals = f(xx);
coeffs = treefun.vals2coeffs(vals);
err = sum(abs(coeffs(end-1:end))) / 2;
vmax = max(abs(vals(:)));
resolved = ( err < tol * max(vmax, 1) );

end
