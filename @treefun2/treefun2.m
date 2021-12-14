classdef treefun2

    properties

        boxes
        n = 16
        domain = [-1 1 -1 1]
        flatNeighbors
        leafNeighbors

    end

    methods

        function f = treefun2(varargin)

            if ( nargin < 1 )
                return
            end

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
                    L = [leaves(f).id];
                    [xx0, yy0] = chebpts2(f.n, f.n, [0 1 0 1]);
                    for k = 1:length(L)
                        id = L(k);
                        dom = f.boxes(id).domain;
                        sclx = diff(dom(1:2));
                        scly = diff(dom(3:4));
                        xx = sclx*xx0 + dom(1);
                        yy = scly*yy0 + dom(3);
                        vals = func(xx,yy);
                        f.boxes(id).coeffs = vals2coeffs(vals);
                    end
                    return
                elseif ( isscalar(varargin{2}) ) % TREEFUN2(F, N)
                    f.n = varargin{2};
                else
                    f.domain = varargin{2};  % TREEFUN2(F, [A B C D])
                end
            elseif ( nargin == 3 )
                f.domain = varargin{2};      % TREEFUN2(F, [A B C D], N)
                f.n = varargin{3};
            elseif ( nargin == 4 )
                f.domain = varargin{2};      % TREEFUN2(F, [A B C D], N, OPTS)
                f.n = varargin{3};
                opts1 = varargin{4};
                if ( isfield(opts1, 'balance') )
                    opts.balance = opts1.balance;
                end
                if ( isfield(opts, 'neighbors') )
                    opts.neighbors = opts1.neighbors;
                end
            end

            func = varargin{1};
            if ( isnumeric(func) && isscalar(func) )
                func = @(x,y) func + 0*x;
            elseif ( isa(func, 'chebfun2') )
                func = @(x,y) feval(func, x, y);
            elseif ( isstruct(func) )
                f.boxes = func;
                f.domain = f.boxes(1).domain;
                f.n = size(f.boxes(end).coeffs, 1);
                f = balance2(f);
                [f.flatNeighbors, f.leafNeighbors] = generateNeighbors(f);
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
            f.boxes(1).row      = 1;

            % f = buildDepthFirst(f, func);
            f = buildBreadthFirst(f, func);

            % Now do level restriction
            if ( opts.balance )
                %f = balance(f, func);
                f = balance2(f);
            else
                % Do a cumulative sum in reverse to correct the heights
                for k = length(f.boxes):-1:1
                   if ( ~isLeaf(f.boxes(k)) )
                       f.boxes(k).height = 1 + max([f.boxes(f.boxes(k).children).height]);
                   end
                end
            end

            if ( opts.neighbors )
                [f.flatNeighbors, f.leafNeighbors] = generateNeighbors(f);
            end

        end

    end

    methods ( Static )

        u = poisson(f);

    end

    methods ( Access = private )

        function f = refineBox(f, id)

            % Split into four child boxes
            dom = f.boxes(id).domain;
            xmid = mean(dom(1:2));
            ymid = mean(dom(3:4));
            parent = f.boxes(id);

            child1 = struct();
            child1.domain   = [dom(1) xmid dom(3) ymid];
            child1.id       = length(f.boxes)+1;
            child1.parent   = id;
            child1.children = [];
            child1.level    = parent.level+1;
            child1.height   = 0;
            child1.coeffs   = [];
            child1.col      = 2*(parent.col-1) + 1;
            child1.row      = 2*(parent.row-1) + 1;
            f.boxes(end+1) = child1;

            child2 = struct();
            child2.domain   = [xmid dom(2) dom(3) ymid];
            child2.id       = length(f.boxes)+1;
            child2.parent   = id;
            child2.children = [];
            child2.level    = parent.level+1;
            child2.height   = 0;
            child2.coeffs   = [];
            child2.col      = 2*(parent.col-1) + 2;
            child2.row      = 2*(parent.row-1) + 1;
            f.boxes(end+1) = child2;

            child3 = struct();
            child3.domain   = [dom(1) xmid ymid dom(4)];
            child3.id       = length(f.boxes)+1;
            child3.parent   = id;
            child3.children = [];
            child3.level    = parent.level+1;
            child3.height   = 0;
            child3.coeffs   = [];
            child3.col      = 2*(parent.col-1) + 1;
            child3.row      = 2*(parent.row-1) + 2;
            f.boxes(end+1) = child3;

            child4 = struct();
            child4.domain   = [xmid dom(2) ymid dom(4)];
            child4.id       = length(f.boxes)+1;
            child4.parent   = id;
            child4.children = [];
            child4.level    = parent.level+1;
            child4.height   = 0;
            child4.coeffs   = [];
            child4.col      = 2*(parent.col-1) + 2;
            child4.row      = 2*(parent.row-1) + 2;
            f.boxes(end+1) = child4;

            f.boxes(id).children = [child1.id, child2.id, ...
                                    child3.id, child4.id];

            %f.boxes(id).children = [child3.id, child4.id, ...
            %                        child2.id, child1.id];
            
            f.boxes(id).height = 1;
            f.boxes(id).coeffs = [];

        end

        function f = balance(f, func)
            L = [leaves(f).id];
            while ( ~isempty(L) )
                id = L(1);
                neighborIDs = neighbors(f, id);
                split = false;

                if ( ~split && ~isnan(neighborIDs(1)) ) % Left neighbor
                    nbr = f.boxes(neighborIDs(1));
                    if ( ~isLeaf(nbr) )
                        childrenIDs = nbr.children([2 4]); % SE, NE
                        split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
                    end
                end

                if ( ~split && ~isnan(neighborIDs(2)) ) % Right neighbor
                    nbr = f.boxes(neighborIDs(2));
                    if ( ~isLeaf(nbr) )
                        childrenIDs = nbr.children([1 3]); % SW, NW
                        split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
                    end
                end

                if ( ~split && ~isnan(neighborIDs(3)) ) % Down neighbor
                    nbr = f.boxes(neighborIDs(3));
                    if ( ~isLeaf(nbr) )
                        childrenIDs = nbr.children([3 4]); % NW, NE
                        split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
                    end
                end

                if ( ~split && ~isnan(neighborIDs(4)) ) % Up neighbor
                    nbr = f.boxes(neighborIDs(4));
                    if ( ~isLeaf(nbr) )
                        childrenIDs = nbr.children([1 2]); % SW, SE
                        split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
                    end
                end
                
                if ( ~split && ~isnan(neighborIDs(5)) ) % Left down corner
                    nbr = f.boxes(neighborIDs(5));
                    if ( ~isLeaf(nbr) )
                        childrenIDs = nbr.children(4); % Right up
                        split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
                    end
                end

                if ( split )
                    % Split into four child boxes
                    f = refineBox(f, id);
                    L = [L f.boxes(id).children]; %#ok<AGROW>
                    neighborIDs = neighborIDs( ~isnan(neighborIDs) );
                    idx = [f.boxes(neighborIDs).level] < f.boxes(id).level;
                    L = [L neighborIDs(idx)]; %#ok<AGROW>
                    L = unique(L, 'stable');
                end

                L(1) = [];
            end

            % Do a cumulative sum in reverse to correct the heights
            for k = length(f.boxes):-1:1
               if ( ~isLeaf(f.boxes(k)) )
                   f.boxes(k).height = 1 + max([f.boxes(f.boxes(k).children).height]);
               end
            end

            % Now fill in the leaf coefficients
            L = [leaves(f).id];
            [xx0, yy0] = chebpts2(f.n, f.n, [0 1 0 1]);
            for k = 1:length(L)
                id = L(k);
                if ( isempty(f.boxes(id).coeffs) )
                    dom = f.boxes(id).domain;
                    sclx = diff(dom(1:2));
                    scly = diff(dom(3:4));
                    xx = sclx*xx0 + dom(1);
                    yy = scly*yy0 + dom(3);
                    vals = func(xx,yy);
                    f.boxes(id).coeffs = vals2coeffs(vals);
                end
            end
        end
        
        function f = balance2(f)
            L = [leaves(f).id];
            while ( ~isempty(L) )
                id = L(1);
                neighborIDs = neighbors(f, id);
                split = false;

                if ( ~split && ~isnan(neighborIDs(1)) ) % Left neighbor
                    nbr = f.boxes(neighborIDs(1));
                    if ( ~isLeaf(nbr) )
                        childrenIDs = nbr.children([2 4]); % SE, NE
                        split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
                    end
                end

                if ( ~split && ~isnan(neighborIDs(2)) ) % Right neighbor
                    nbr = f.boxes(neighborIDs(2));
                    if ( ~isLeaf(nbr) )
                        childrenIDs = nbr.children([1 3]); % SW, NW
                        split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
                    end
                end

                if ( ~split && ~isnan(neighborIDs(3)) ) % Down neighbor
                    nbr = f.boxes(neighborIDs(3));
                    if ( ~isLeaf(nbr) )
                        childrenIDs = nbr.children([3 4]); % NW, NE
                        split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
                    end
                end

                if ( ~split && ~isnan(neighborIDs(4)) ) % Up neighbor
                    nbr = f.boxes(neighborIDs(4));
                    if ( ~isLeaf(nbr) )
                        childrenIDs = nbr.children([1 2]); % SW, SE
                        split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
                    end
                end

                if ( ~split && ~isnan(neighborIDs(5)) ) % Left down corner
                    nbr = f.boxes(neighborIDs(5));
                    if ( ~isLeaf(nbr) )
                        childrenIDs = nbr.children(4); % Right up
                        split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
                    end
                end
                
                if ( ~split && ~isnan(neighborIDs(6)) ) % Right down corner
                    nbr = f.boxes(neighborIDs(6));
                    if ( ~isLeaf(nbr) )
                        childrenIDs = nbr.children(3); % Left up
                        split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
                    end
                end
                
                if ( ~split && ~isnan(neighborIDs(7)) ) % Left up corner
                    nbr = f.boxes(neighborIDs(7));
                    if ( ~isLeaf(nbr) )
                        childrenIDs = nbr.children(2); % Right down
                        split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
                    end
                end

                if ( ~split && ~isnan(neighborIDs(8)) ) % Right up corner
                    nbr = f.boxes(neighborIDs(8));
                    if ( ~isLeaf(nbr) )
                        childrenIDs = nbr.children(1); % Left down
                        split = split | any( ~isLeaf(f.boxes(childrenIDs)) );
                    end
                end

                if ( split )
                    % This was a leaf, so we'll use its coeffs to evaluate
                    % on the new children
                    coeffs = f.boxes(id).coeffs;
                    % Split into four child boxes
                    f = refineBox(f, id);
                    children = f.boxes(id).children;
                    [LL, LR, UL, UR] = coeffs2children(coeffs);
                    f.boxes(children(1)).coeffs = LL; % Lower left
                    f.boxes(children(2)).coeffs = LR; % Lower right
                    f.boxes(children(3)).coeffs = UL; % Upper left
                    f.boxes(children(4)).coeffs = UR; % Upper right
                    L = [L children]; %#ok<AGROW>
                    neighborIDs = neighborIDs( ~isnan(neighborIDs) );
                    idx = [f.boxes(neighborIDs).level] < f.boxes(id).level;
                    L = [L neighborIDs(idx)]; %#ok<AGROW>
                    L = unique(L, 'stable');
                end

                L(1) = [];
            end

            % Do a cumulative sum in reverse to correct the heights
            for k = length(f.boxes):-1:1
               if ( ~isLeaf(f.boxes(k)) )
                   f.boxes(k).height = 1 + max([f.boxes(f.boxes(k).children).height]);
               end
            end
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
                    % Split into four child boxes
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

        function f = buildDepthFirst(f, func, id, level)

            if ( nargin == 2 )
                id = 1;
                level = 0;
            end

            f.boxes(id).level = level;
            f.boxes(id).height = 0;

            [resolved, coeffs] = isResolved(func, f.boxes(id).domain, f.n);

            if ( resolved )
                f.boxes(id).coeffs = coeffs;
            else
                % Split into four child boxes
                f = refineBox(f, id);

                % Recurse
                f = buildDepthFirst(f, func, f.boxes(id).children(1), level+1);
                f = buildDepthFirst(f, func, f.boxes(id).children(2), level+1);
                f = buildDepthFirst(f, func, f.boxes(id).children(3), level+1);
                f = buildDepthFirst(f, func, f.boxes(id).children(4), level+1);

                % Set height
                f.boxes(id).height = 1 + max([f.boxes(f.boxes(id).children).height]);
            end

        end

    end
    
    methods ( Static )

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
coeffs = vals2coeffs(vals);
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
