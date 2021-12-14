classdef treefun

    properties

        boxes
        n = 16
        domain = [-1 1]

    end

    methods

        function f = treefun(varargin)

            if ( nargin < 1 )
                return
            end

            if ( nargin == 2 )
                if ( isa(varargin{2}, 'treefun') ) % TREEFUN(F, TF)
                    % We were given the tree structure
                    f = varargin{2};
                    func = varargin{1};
                    if ( isnumeric(func) && isscalar(func) )
                        func = @(x) func + 0*x;
                    elseif ( isa(func, 'chebfun') )
                        func = @(x) feval(func, x);
                    end
                    % We just need to fill in the leaf coefficients
                    L = [leaves(f).id];
                    x0 = chebpts(f.n, [0 1]);
                    for k = 1:length(L)
                        id = L(k);
                        dom = f.boxes(id).domain;
                        sclx = diff(dom);
                        x = sclx*x0 + dom(1);
                        vals = func(x);
                        f.boxes(id).coeffs = vals2coeffs(vals);
                    end
                    return
                elseif ( isscalar(varargin{2}) ) % TREEFUN(F, N)
                    f.n = varargin{2};
                else
                    f.domain = varargin{2};  % TREEFUN(F, [A B])
                end
            elseif ( nargin == 3 )
                f.domain = varargin{2};      % TREEFUN(F, [A B], N)
                f.n = varargin{3};
            end

            func = varargin{1};
            if ( isnumeric(func) && isscalar(func) )
                func = @(x) func + 0*x;
            elseif ( isa(func, 'chebfun') )
                func = @(x) feval(func, x);
            elseif ( isstruct(func) )
                f.boxes = func;
                f.domain = f.boxes(1).domain;
                f.n = size(f.boxes(end).coeffs, 1);
                %f = balance2(f);
                %[f.flatNeighbors, f.leafNeighbors] = generateNeighbors(f);
                return
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

            % f = buildDepthFirst(f, func);
            f = buildBreadthFirst(f, func);

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
            child1.col      = 2*(parent.col-1) + 1;
            f.boxes(end+1) = child1;

            child2 = struct();
            child2.domain   = [xmid dom(2)];
            child2.id       = length(f.boxes)+1;
            child2.parent   = id;
            child2.children = [];
            child2.level    = parent.level+1;
            child2.height   = 0;
            child2.coeffs   = [];
            child2.col      = 2*(parent.col-1) + 2;
            f.boxes(end+1) = child2;

            f.boxes(id).children = [child1.id, child2.id];
            f.boxes(id).height = 1;
            f.boxes(id).coeffs = [];

        end

        function f = buildBreadthFirst(f, func)

            % Note: the length changes at each iteration here
            id = 1;
            while ( id <= length(f.boxes) )
                [resolved, coeffs] = isResolved(func, f.boxes(id).domain, f.n);
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

end

function [resolved, coeffs] = isResolved(f, dom, n)

persistent xx0 xxx0 nstored

tol = 1e-10;
nrefpts = 2*n; % Sample at equispaced points to test error

if ( isempty(xx0) || isempty(xxx0) || n ~= nstored )
    nstored = n;
    xx0 = chebpts(n, [0 1 0 1]);
    xxx0 = linspace(0, 1, nrefpts).';
end
sclx = diff(dom);
xx  = sclx*xx0  + dom(1);
xxx = sclx*xxx0 + dom(1);

vals = f(xx);
coeffs = vals2coeffs(vals);
err_cfs = sum(abs(coeffs(end-1:end))) / 2;

F = f(xxx);
G = coeffs2refvals(coeffs);
err_vals = norm(F(:) - G(:), inf);

% err = err_cfs;
err = err_vals;

% tech = chebtech2;
% tech.coeffs = coeffs;
% pref = tech.techPref();
% pref.chebfuneps = 1e-1;
% data = struct();
% data = chebtech.parseDataInputs(data, pref);
% data.vscale = max(data.vscale, vscale(tech));
% ishappy = standardCheck(tech, vals, data, pref);
% ishappy

vmax = max(abs(vals(:)));
resolved = ( err < tol * max(vmax, 1) );

end