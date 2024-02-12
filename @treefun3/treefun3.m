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
        rint

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
            opts.tol = 1e-12;
            opts.checkpts = [];
            
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
                if ( isfield(opts1, 'tol') )
                    opts.tol = opts1.tol;
                end
                if ( isfield(opts1, 'checkpts') )
                    opts.checkpts = opts1.checkpts;
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

            %%% rint %%%
            nalias = f.n;
            x0 = (1-cos(pi*(2*(1:nalias)'-1)/(2*nalias)))/2;
            [xx0, yy0, zz0] = ndgrid(x0);
            l = floor(nalias/2)+1;
            v = [2*exp(1i*pi*(0:nalias-l)/nalias)./(1-4*(0:nalias-l).^2)  zeros(1,l)];
            w0 = real(ifft(v(1:nalias) + conj(v(nalias+1:-1:2))))'/2;
            [wx0, wy0, wz0] = ndgrid(w0); 
            sclx = diff(dom(1:2));
            scly = diff(dom(3:4));
            sclz = diff(dom(5:6));
            xx = sclx*xx0 + dom(1); 
            yy = scly*yy0 + dom(3); 
            zz = sclz*zz0 + dom(5); 
            ww0 = wx0.*wy0.*wz0;
            vals = func(xx,yy,zz);
            f.rint = max(sqrt((sclx*scly*sclz)*sum(vals.^2.*ww0, 'all')),1e-16); % initialze l2
            

            % f = buildDepthFirst(f, func);
            f = buildBreadthFirst(f, func, opts.tol, opts.checkpts);
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

        function f = buildBreadthFirst(f, func, tol, checkpts)

            persistent xx0 yy0 zz0 ww0 nstored
            nalias = f.n;
            if ( isempty(xx0) || f.n ~= nstored )
                nstored = f.n;
                x0 = (1-cos(pi*(2*(1:nalias)'-1)/(2*nalias)))/2;
                [xx0, yy0, zz0] = ndgrid(x0);
                l = floor(nalias/2)+1;
                v = [2*exp(1i*pi*(0:nalias-l)/nalias)./(1-4*(0:nalias-l).^2)  zeros(1,l)];
                w0 = real(ifft(v(1:nalias) + conj(v(nalias+1:-1:2))))'/2;
                [wx0, wy0, wz0] = ndgrid(w0); 
                ww0 = wx0.*wy0.*wz0;
            end

            % Note: the length changes at each iteration here
            id = 1;
            while ( id <= length(f.id) )
                [resolved, coeffs, rintvals] = isResolved(func, f.domain(:,id), f.n, tol, checkpts, f.rint);
                if ( resolved )
                    f.coeffs{id} = coeffs;
                    f.height(id) = 0;
                else
                    % Split into eight child boxes
                    f = refineBox(f, id);
                    f.height(id) = 1;
                    % whenever split, update rint to be more accurate ...
                    % (f.rint)^2 - (rintvals)^2 + eight child (rintvals)^2
                    r8intvals = 0;
                    for k = 1:8
                      dom = f.domain(:,f.children(k,id));
                      sclx = diff(dom(1:2));
                      scly = diff(dom(3:4));
                      sclz = diff(dom(5:6));
                      xx = sclx*xx0 + dom(1); 
                      yy = scly*yy0 + dom(3); 
                      zz = sclz*zz0 + dom(5); 
                      r8intvals = r8intvals + (sclx*scly*sclz)*sum(func(xx,yy,zz).^2.*ww0, 'all');
                    end
                    f.rint = sqrt(f.rint^2 - rintvals^2 + r8intvals);
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
        checkvals = coeffs2checkvals(coeffs,x,y,z);

    end
    

end

function [resolved, coeffs, rintvals] = isResolved(f, dom, n, tol, checkpts, rint)

persistent xx0 yy0 zz0 ww0 xxx0 yyy0 zzz0 nstored

nalias = n;
nrefpts = 2*n; % Sample at equispaced points to test error

if ( isempty(xx0) || isempty(xxx0) || n ~= nstored )
    nstored = n;
    x0 = (1-cos(pi*(2*(1:nalias)'-1)/(2*nalias)))/2;
    [xx0, yy0, zz0] = ndgrid(x0);
    [xxx0, yyy0, zzz0] = ndgrid(linspace(0, 1, nrefpts));
    l = floor(nalias/2)+1;
    v = [2*exp(1i*pi*(0:nalias-l)/nalias)./(1-4*(0:nalias-l).^2)  zeros(1,l)];
    w0 = real(ifft(v(1:nalias) + conj(v(nalias+1:-1:2))))'/2;
    [wx0, wy0, wz0] = ndgrid(w0); 
    ww0 = wx0.*wy0.*wz0;
end
sclx = diff(dom(1:2));
scly = diff(dom(3:4));
sclz = diff(dom(5:6));
xx = sclx*xx0 + dom(1); xxx = sclx*xxx0 + dom(1);
yy = scly*yy0 + dom(3); yyy = scly*yyy0 + dom(3);
zz = sclz*zz0 + dom(5); zzz = sclz*zzz0 + dom(5);

vals = f(xx,yy,zz);
rintvals = sqrt((sclx*scly*sclz)*sum(vals.^2.*ww0, 'all'));
coeffs = treefun3.vals2coeffs(vals);
coeffs = coeffs(1:n,1:n,1:n);
h = sclx;
eta = 0;

vmax = max(abs(vals(:)));

if 0
  Ex = sum(abs(coeffs(end-1:end,:,:)), 'all') / (3*n^2);
  Ey = sum(abs(coeffs(:,end-1:end,:)), 'all') / (3*n^2);
  Ez = sum(abs(coeffs(:,:,end-1:end)), 'all') / (3*n^2);
  err_cfs = (Ex + Ey + Ez) / 3;
  
  % F = f(xxx,yyy,zzz);
  % G = coeffs2refvals(coeffs);
  % err_vals = max(abs(F(:) - G(:)));
  
  err = err_cfs;
  %err = min(err_cfs, err_vals);
  
  resolved = ( err * h^eta < tol * max(vmax, 1) );
else
  erra = sqrt(sum(coeffs.^2, 'all') - sum(coeffs(1:end-2,1:end-2,1:end-2).^2, 'all')) / (n^3 - (n-2)^3);
  resolved = ( erra < tol* sqrt(1/(sclx*scly*sclz)) * rint );
end

if ( ~isempty(checkpts) ) % check if func values @ checkpts agree
  % func = @(x,y,z) exp(-(x.^2+y.^2+z.^2)*5000000) + exp(-((x-1/2).^2+(y-1/3).^2+(z-3/5).^2)*1000000) + exp(-((x+1/2).^2+(y+1/3).^2+(z+3/5).^2)*2000000);
  % f = treefun3(func,[-2 2 -2 2 -2 2],10,struct('checkpts',[0 1/2 -1/2;0 1/3 -1/3;0 3/5 -3/5])); vs f = treefun3(func,[-2 2 -2 2 -2 2],10);
  % plot(f,func)
  % for later, reduce checkpts if resolved
  xxx = 2 * ((checkpts(1,:)' - dom(1))/sclx) - 1;
  yyy = 2 * ((checkpts(2,:)' - dom(3))/sclx) - 1;
  zzz = 2 * ((checkpts(3,:)' - dom(5))/sclx) - 1;
  in = ( xxx>=-1 & xxx<=1 & ...
         yyy>=-1 & yyy<=1 & ...
         zzz>=-1 & zzz<=1);
  if ( any(in) )
    F = f(checkpts(1,in),checkpts(2,in),checkpts(3,in));
    G = treefun3.coeffs2checkvals(coeffs,xxx(in),yyy(in),zzz(in));
    err_checkvals = max(abs(F(:) - G(:)));
    resolved = resolved && ( err_checkvals * h^eta < tol * max(vmax, 1) );
  end
  
end

end
