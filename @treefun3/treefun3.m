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
        dep
        morton
        flatNeighbors
        leafNeighbors
        rint0
        rint
        vmax

    end

    methods

        function f = treefun3(varargin)

            if ( nargin < 1 )
                return
            end

            dom = [-1 1 -1 1 -1 1];
            opts = struct();
            opts.balance = false; % fortran way exists, morton way intial attempt in BBBBBBBoston
            opts.neighbors = false;
            opts.tol = 1e-12;
            opts.checkpts = [];
            opts.ifcoeffs = true; 
            opts.ifstorecoeffs = true;
            
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

            % initialize vals, rint, coeffs... wrap this up?
            nalias = f.n;
            nd = numel(func);
            if nd == 1 && ~iscell(func)
              func1 = []; 
              func1{1} = func;
              func = func1;
            end
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
            vals = cell(1,nd);
            for k = 1:nd
              vals{k} = func{k}(xx,yy,zz);
            end
            vals = cat(4,vals{:});
            coeffs = treefun3.vals2coeffs(vals);
            rint = max(squeeze(sqrt((sclx*scly*sclz)*sum(vals.^2.*ww0, [1 2 3]))),1e-16); % initialze l2
            
            % 
            f.domain(:,1)   = dom(:);
            f.level(1)      = 0;
            f.height(1)     = 0;
            f.id(1)         = 1;
            f.parent(1)     = 0;
            f.children(:,1) = zeros(8, 1);
            f.coeffs{1}     = [];
            f.col           = uint64(0);
            f.row           = uint64(0);
            f.dep           = uint64(0); % is this the correct generalization?
            f.coeffs{1}     = coeffs(1:f.n,1:f.n,1:f.n,:); 
            f.rint0(:,1)    = zeros(nd,1); 
            f.rint(:,1)     = rint;
            f.vmax(:,1)     = squeeze(max(abs(vals),[],[1 2 3]));
            
            

            % f = buildDepthFirst(f, func);
            f = buildBreadthFirst(f, func, opts.tol, opts.checkpts, opts.ifcoeffs, opts.ifstorecoeffs);
            f.morton = cartesian2morton(f.col, f.row, f.dep);

            % Now do level restriction
            if ( opts.balance )
                f = balance(f);
                % f = balancef(f);
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

        f = refineBox(f, id, func);

        function f = buildBreadthFirst(f, func, tol, checkpts, ifcoeffs, ifstorecoeffs)

            % initialization a 2nd tree that's the same as the 1st one
            rint = f.rint(:,1);
            [~,~,~,nd] = size(f.coeffs{1});
            resolved = isResolved(f.coeffs{1}, f.domain(:,1), f.n, tol, f.vmax(:,1), func, checkpts, ifcoeffs, rint);
            f.rint0(:,1) = rint;
            if resolved
              return
            else
              refinecnt = 1;
              idl_start = 1; 
              idl = 1; % current level box id
            end
            if ~ifstorecoeffs  % save memory, once isResolved done
              f.coeffs{1} = [];
            end

            while(refinecnt)
            
              % process previous level unresolved boxes, idl(1), idl(2),....
              domc    = zeros(6,8*refinecnt); % allocate space for all children info of this level's boxes
              coeffsc = cell(1,8*refinecnt);
              rintc   = zeros(nd,8*refinecnt);
              rint0c  = zeros(nd,8*refinecnt);
              vmaxc   = zeros(nd,8*refinecnt);
              parentc = zeros(1,8*refinecnt);
              idc     = zeros(1,8*refinecnt);
              levelc  = zeros(1,8*refinecnt);
              colc    = uint64(zeros(1,8*refinecnt));
              rowc    = uint64(zeros(1,8*refinecnt));
              depc    = uint64(zeros(1,8*refinecnt));
              for k = 1:refinecnt % 1 by 1
                id = idl(k);
                [domck,coeffsck,rintck,vmaxck,colck,rowck,depck] = refineBoxv2(f.domain(:,id), f.n, func, f.col(id), f.row(id), f.dep(id)); % just refine the box...
                f.children(:,id) = idl_start+8*(k-1)+(1:8); % start to know 
                f.height(id)     = 1;
                % 1 to 8
                cidx             = 8*(k-1)+(1:8);
                domc(:,cidx)     = domck;
                coeffsc(cidx)    = coeffsck;
                rintc(:,cidx)    = rintck;
                vmaxc(:,cidx)    = vmaxck;
                parentc(cidx)    = id*ones(1,8); % these children's parent id
                idc(cidx)        = idl_start+8*(k-1)+(1:8); % children self id
                levelc(cidx)     = (f.level(id)+1)*ones(1,8);
                colc(cidx)       = colck;
                rowc(cidx)       = rowck;
                depc(cidx)       = depck;
              end
              childrenc  = zeros(8,8*refinecnt); % don't know yet, a bunch of 0s to be concatenate to the end of f.children, will know when next level
              heightc    = zeros(1,8*refinecnt); % tmp leaf box
              f.domain   = cat(2,f.domain,domc); % concatenate children info
              f.coeffs   = cat(2,f.coeffs,coeffsc);
              f.rint     = cat(2,f.rint,rintc);
              f.rint0    = cat(2,f.rint0,rint0c);
              f.vmax     = cat(2,f.vmax,vmaxc);
              f.parent   = cat(2,f.parent,parentc);
              f.id       = cat(2,f.id,idc);
              f.level    = cat(2,f.level,levelc);
              f.children = cat(2,f.children,childrenc);
              f.height   = cat(2,f.height,heightc);
              f.col      = cat(2,f.col,colc);
              f.row      = cat(2,f.row,rowc);
              f.dep      = cat(2,f.dep,depc);
            
              % update rint, - num of idl parents contribution + 8*num of idl children contribution
              rint = sqrt(rint.^2 - sum(f.rint(:,idl).^2,2) + sum(f.rint(:,(idl_start+1):(idl_start+8*refinecnt)).^2,2)); 
            
              % see if newly created boxes need to be refined next
              unresolvedl = true(1,8*length(idl)); % potentially unresolved
              idl0 = (idl_start+1):(idl_start+8*length(idl)); % current level id, from start to num of idl refined
              for k = 1:8*length(idl)
                id = idl_start + k;
                resolved = isResolved(f.coeffs{id}, f.domain(:,id), f.n, tol, f.vmax(:,id), func, checkpts, ifcoeffs, rint);
                f.rint0(:,id) = rint;
                if resolved
                  unresolvedl(k) = false;
                end
                if ~ifstorecoeffs  % save memory, once isResolved done
                  f.coeffs{id} = [];
                end
              end
            
              % for next loop
              idl_start = idl_start + 8*refinecnt; % next level starts from here
              idl = idl0(unresolvedl); % select subset of all potential boxes, based on unresolved
              refinecnt = length(idl);
              
            end
          
            % % Note: the length changes at each iteration here
            % id = 1;
            % rint = f.rint(:,1);
            % while ( id <= length(f.id) )
            %     resolved = isResolved(f.coeffs{id}, f.domain(:,id), f.n, tol, f.vmax(:,id), func, checkpts, ifcoeffs, rint);
            %     if ( resolved )
            %         f.height(id) = 0;
            %     else
            %         % Split into eight child boxes
            %         f = refineBox(f, id, func);
            %         f.height(id) = 1;
            %         f.coeffs{id} = []; % delete coeffs
            %         rint = sqrt(rint.^2 - f.rint(:,id).^2 + sum(f.rint(:,end-7:end).^2,2)); % whenever split, update rint to be more accurate ...
            %     end
            %     id = id + 1;
            % end

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
        vals = coeffs2refvals(coeffs);
        % refvals = chebvals2refvals(chebvals);
        checkvals = coeffs2checkvals(coeffs,x,y,z);

    end
    

end

function resolved = isResolved(coeffs, dom, n, tol, vmax, f, checkpts, ifcoeffs, rint)
% need f if checkpts

persistent xxx0 yyy0 zzz0 www0 nstored

% nalias = n;
nrefpts = n; % Sample at equispaced points to test error

if ( isempty(xxx0) || n ~= nstored )
    nstored = n;
    [xxx0, yyy0, zzz0] = ndgrid(linspace(0, 1, nrefpts));
    [wxxx0, wyyy0, wzzz0] = ndgrid(1/nstored*ones(nstored,1));
    www0 = wxxx0.*wyyy0.*wzzz0; % weight
end

sclx = diff(dom(1:2));
scly = diff(dom(3:4));
sclz = diff(dom(5:6));

h = sclx;
eta = 0;
[~,~,~,nd] = size(coeffs); % or rint...

if ~ifcoeffs
  resolved = 1;
  xxx = sclx*xxx0 + dom(1); 
  yyy = scly*yyy0 + dom(3); 
  zzz = sclz*zzz0 + dom(5); 
  vals = cell(1,nd);
  for k = 1:nd
    vals{k} = f{k}(xxx,yyy,zzz);
  end
  F = cat(4,vals{:});
  G = treefun3.coeffs2refvals(coeffs); % structured grid
  for k = 1:nd
    erra = sqrt(sum(squeeze(G(:,:,:,k) - F(:,:,:,k)).^2.*www0,'all'));
    resolved = resolved && ( erra < tol * sqrt(1/(sclx*scly*sclz)) * rint(k) );
  end

elseif ifcoeffs
  resolved = 1;
  erra = zeros(nd,1);
  for k = 1:nd
    erra(k) = sqrt(  sum(coeffs(end,:,:,k).^2,'all') ...
                + sum(coeffs(1:end-1,end,:,k).^2,'all') ...
                + sum(coeffs(1:end-1,1:end-1,end,k).^2,'all')) / sqrt(n^3 - (n-1)^3); % only last slice
    % erra = sqrt(  sum(coeffs(end-1:end,:,:,k).^2,'all') ...
    %             + sum(coeffs(1:end-2,end-1:end,:,k).^2,'all') ...
    %             + sum(coeffs(1:end-2,1:end-2,end-1:end,k).^2,'all')) / sqrt(n^3 - (n-2)^3); 
    resolved = resolved && ( erra(k) < tol* sqrt(1/(sclx*scly*sclz)) * rint(k) );
  end

else % did not check
  vmax = max(abs(vals(:))); % if needed, compute and store when refineBox

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
end

if ( ~isempty(checkpts) ) % check if func values @ checkpts agree
  % func = @(x,y,z) exp(-(x.^2+y.^2+z.^2)*5000000) + exp(-((x-1/2).^2+(y-1/3).^2+(z-3/5).^2)*1000000) + exp(-((x+1/2).^2+(y+1/3).^2+(z+3/5).^2)*2000000);
  % f = treefun3(func,[-2 2 -2 2 -2 2],10,struct('checkpts',[0 1/2 -1/2;0 1/3 -1/3;0 3/5 -3/5])); vs f = treefun3(func,[-2 2 -2 2 -2 2],10);
  % plot(f,func)
  % for later, reduce checkpts if resolved
  xxx = 2 * ((checkpts(1,:)' - dom(1))/sclx) - 1;
  yyy = 2 * ((checkpts(2,:)' - dom(3))/scly) - 1;
  zzz = 2 * ((checkpts(3,:)' - dom(5))/sclz) - 1;
  in = ( xxx>=-1 & xxx<=1 & ...
         yyy>=-1 & yyy<=1 & ...
         zzz>=-1 & zzz<=1);
  if ( any(in) )
    F = zeros(nd,sum(in));
    for k = 1:nd
      F(k,:) = f{k}(checkpts(1,in)',checkpts(2,in)',checkpts(3,in)')'; % nd x n checkpts
    end
    G = treefun3.coeffs2checkvals(coeffs,xxx(in),yyy(in),zzz(in));
    err_checkvals = max(abs(F - G),[],2);
    for k = 1:nd
      resolved = resolved && ( err_checkvals(k) * h^eta < tol * max(vmax(k), 1) );
    end
  end
end

end


function [domain,coeffs,rint,vmax,ccol,crow,cdep] = refineBoxv2(dom, n, func, col, row, dep)

% if nargin==0, test_refineBox3dv2; return; end

persistent xx0 yy0 zz0 ww0 nstored
nalias = n;
if ( isempty(xx0) || n ~= nstored )
    nstored = n;
    x0 = (1-cos(pi*(2*(1:nalias)'-1)/(2*nalias)))/2;
    [xx0, yy0, zz0] = ndgrid(x0);
    l = floor(nalias/2)+1;
    v = [2*exp(1i*pi*(0:nalias-l)/nalias)./(1-4*(0:nalias-l).^2)  zeros(1,l)];
    w0 = real(ifft(v(1:nalias) + conj(v(nalias+1:-1:2))))'/2;
    [wx0, wy0, wz0] = ndgrid(w0); 
    ww0 = wx0.*wy0.*wz0;
end

% Split into eight child boxes
xmid = mean(dom(1:2));
ymid = mean(dom(3:4));
zmid = mean(dom(5:6));

% domain and coeffs
domain = zeros(6,8);
coeffs = cell(1,8);
nd = numel(func);

% morton related
ccol = uint64(zeros(1,8));
crow = uint64(zeros(1,8));
cdep = uint64(zeros(1,8));


cdom1              = [dom(1) xmid dom(3) ymid dom(5) zmid];
csclx1             = diff(cdom1(1:2));
cscly1             = diff(cdom1(3:4));
csclz1             = diff(cdom1(5:6));
cxx1               = csclx1*xx0 + cdom1(1); 
cyy1               = cscly1*yy0 + cdom1(3); 
czz1               = csclz1*zz0 + cdom1(5); 
cvals1             = cell(1,nd);
for k = 1:nd
  cvals1{k}        = func{k}(cxx1,cyy1,czz1);
end
cvals1             = cat(4,cvals1{:});
ccoeffs1           = treefun3.vals2coeffs(cvals1);
ccol(1)            = 2*col;
crow(1)            = 2*row;
cdep(1)            = 2*dep;
% f.morton(cid1)     = cartesian2morton(f.col(cid1), f.row(cid1));
domain(:,1)        = cdom1;
coeffs{1}          = ccoeffs1(1:n,1:n,1:n,:); % to replace f.coeffs{cid1} = [];

cdom2              = [xmid dom(2) dom(3) ymid dom(5) zmid];
csclx2             = diff(cdom2(1:2));
cscly2             = diff(cdom2(3:4));
csclz2             = diff(cdom2(5:6));
cxx2               = csclx2*xx0 + cdom2(1); 
cyy2               = cscly2*yy0 + cdom2(3); 
czz2               = csclz2*zz0 + cdom2(5); 
cvals2             = cell(1,nd);
for k = 1:nd
  cvals2{k}        = func{k}(cxx2,cyy2,czz2);
end
cvals2             = cat(4,cvals2{:});
ccoeffs2           = treefun3.vals2coeffs(cvals2);
ccol(2)            = 2*col + 1;
crow(2)            = 2*row;
cdep(2)            = 2*dep;
% f.morton(cid2)     = cartesian2morton(f.col(cid2), f.row(cid2)); 
domain(:,2)        = cdom2;
coeffs{2}          = ccoeffs2(1:n,1:n,1:n,:); 

cdom3              = [dom(1) xmid ymid dom(4) dom(5) zmid];
csclx3             = diff(cdom3(1:2));
cscly3             = diff(cdom3(3:4));
csclz3             = diff(cdom3(5:6));
cxx3               = csclx3*xx0 + cdom3(1); 
cyy3               = cscly3*yy0 + cdom3(3); 
czz3               = csclz3*zz0 + cdom3(5); 
cvals3             = cell(1,nd);
for k = 1:nd
  cvals3{k}        = func{k}(cxx3,cyy3,czz3);
end
cvals3             = cat(4,cvals3{:});
ccoeffs3           = treefun3.vals2coeffs(cvals3);
ccol(3)            = 2*col;
crow(3)            = 2*row + 1;
cdep(3)            = 2*dep;
% f.morton(cid3)     = cartesian2morton(f.col(cid3), f.row(cid3));
domain(:,3)        = cdom3;
coeffs{3}          = ccoeffs3(1:n,1:n,1:n,:); 

cdom4              = [xmid dom(2) ymid dom(4) dom(5) zmid];
csclx4             = diff(cdom4(1:2));
cscly4             = diff(cdom4(3:4));
csclz4             = diff(cdom4(5:6));
cxx4               = csclx4*xx0 + cdom4(1); 
cyy4               = cscly4*yy0 + cdom4(3); 
czz4               = csclz4*zz0 + cdom4(5); 
cvals4             = cell(1,nd);
for k = 1:nd
  cvals4{k}        = func{k}(cxx4,cyy4,czz4);
end
cvals4             = cat(4,cvals4{:});
ccoeffs4           = treefun3.vals2coeffs(cvals4);
ccol(4)            = 2*col + 1;
crow(4)            = 2*row + 1;
cdep(4)            = 2*dep;
% f.morton(cid4)     = cartesian2morton(f.col(cid4), f.row(cid4));
domain(:,4)        = cdom4;
coeffs{4}          = ccoeffs4(1:n,1:n,1:n,:); 

cdom5              = [dom(1) xmid dom(3) ymid zmid dom(6)];
csclx5             = diff(cdom5(1:2));
cscly5             = diff(cdom5(3:4));
csclz5             = diff(cdom5(5:6));
cxx5               = csclx5*xx0 + cdom5(1); 
cyy5               = cscly5*yy0 + cdom5(3); 
czz5               = csclz5*zz0 + cdom5(5); 
cvals5             = cell(1,nd);
for k = 1:nd
  cvals5{k}        = func{k}(cxx5,cyy5,czz5);
end
cvals5             = cat(4,cvals5{:});
ccoeffs5           = treefun3.vals2coeffs(cvals5);
ccol(5)            = 2*col;
crow(5)            = 2*row;
cdep(5)            = 2*dep + 1;
%        
domain(:,5)        = cdom5;
coeffs{5}          = ccoeffs5(1:n,1:n,1:n,:); 

cdom6              = [xmid dom(2) dom(3) ymid zmid dom(6)];
csclx6             = diff(cdom6(1:2));
cscly6             = diff(cdom6(3:4));
csclz6             = diff(cdom6(5:6));
cxx6               = csclx6*xx0 + cdom6(1); 
cyy6               = cscly6*yy0 + cdom6(3); 
czz6               = csclz6*zz0 + cdom6(5); 
cvals6             = cell(1,nd);
for k = 1:nd
  cvals6{k}        = func{k}(cxx6,cyy6,czz6);
end
cvals6             = cat(4,cvals6{:});
ccoeffs6           = treefun3.vals2coeffs(cvals6);
ccol(6)            = 2*col + 1;
crow(6)            = 2*row;
cdep(6)            = 2*dep + 1;
%
domain(:,6)        = cdom6;
coeffs{6}          = ccoeffs6(1:n,1:n,1:n,:); 

cdom7              = [dom(1) xmid ymid dom(4) zmid dom(6)];
csclx7             = diff(cdom7(1:2));
cscly7             = diff(cdom7(3:4));
csclz7             = diff(cdom7(5:6));
cxx7               = csclx7*xx0 + cdom7(1); 
cyy7               = cscly7*yy0 + cdom7(3); 
czz7               = csclz7*zz0 + cdom7(5); 
cvals7             = cell(1,nd);
for k = 1:nd
  cvals7{k}        = func{k}(cxx7,cyy7,czz7);
end
cvals7             = cat(4,cvals7{:});
ccoeffs7           = treefun3.vals2coeffs(cvals7);
ccol(7)            = 2*col;
crow(7)            = 2*row + 1;
cdep(7)            = 2*dep + 1;
%        
domain(:,7)        = cdom7;
coeffs{7}          = ccoeffs7(1:n,1:n,1:n,:); 

cdom8              = [xmid dom(2) ymid dom(4) zmid dom(6)];
csclx8             = diff(cdom8(1:2));
cscly8             = diff(cdom8(3:4));
csclz8             = diff(cdom8(5:6));
cxx8               = csclx8*xx0 + cdom8(1); 
cyy8               = cscly8*yy0 + cdom8(3); 
czz8               = csclz8*zz0 + cdom8(5); 
cvals8             = cell(1,nd);
for k = 1:nd
  cvals8{k}        = func{k}(cxx8,cyy8,czz8);
end
cvals8             = cat(4,cvals8{:});
ccoeffs8           = treefun3.vals2coeffs(cvals8);
ccol(8)            = 2*col + 1;
crow(8)            = 2*row + 1;
cdep(8)            = 2*dep + 1;
domain(:,8)        = cdom8;
coeffs{8}          = ccoeffs8(1:n,1:n,1:n,:); 

% a few other things
rint = [squeeze(sqrt((csclx1*cscly1*csclz1)*sum(cvals1.^2.*ww0, [1 2 3]))),...
        squeeze(sqrt((csclx2*cscly2*csclz2)*sum(cvals2.^2.*ww0, [1 2 3]))),...
        squeeze(sqrt((csclx3*cscly3*csclz3)*sum(cvals3.^2.*ww0, [1 2 3]))),...
        squeeze(sqrt((csclx4*cscly4*csclz4)*sum(cvals4.^2.*ww0, [1 2 3]))),...
        squeeze(sqrt((csclx5*cscly5*csclz5)*sum(cvals5.^2.*ww0, [1 2 3]))),...
        squeeze(sqrt((csclx6*cscly6*csclz6)*sum(cvals6.^2.*ww0, [1 2 3]))),...
        squeeze(sqrt((csclx7*cscly7*csclz7)*sum(cvals7.^2.*ww0, [1 2 3]))),...
        squeeze(sqrt((csclx8*cscly8*csclz8)*sum(cvals8.^2.*ww0, [1 2 3])))];
vmax = [squeeze(max(abs(cvals1),[],[1 2 3])),...
        squeeze(max(abs(cvals2),[],[1 2 3])),...
        squeeze(max(abs(cvals3),[],[1 2 3])),...
        squeeze(max(abs(cvals4),[],[1 2 3])),...
        squeeze(max(abs(cvals5),[],[1 2 3])),...
        squeeze(max(abs(cvals6),[],[1 2 3])),...
        squeeze(max(abs(cvals7),[],[1 2 3])),...
        squeeze(max(abs(cvals8),[],[1 2 3]))];

end

