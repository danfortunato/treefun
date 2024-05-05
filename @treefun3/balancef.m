function g = balancef(f)
% fortran vol_tree_fix_lr
% assume square boxes? therefore scale is computed from dom(1:2)
% nd needs to be fixed later...
% flatNeighbors and leafNeighbors also need to be updated in the future
%

nd = 1; % 
dpars = zeros(1000,1);
ipars = zeros(100,1);
zpars = zeros(10,1);

if ( isempty(f) )
    g = f;
    return
end

% dada
dom = f.domain(:,1)';
scale = (dom(2)-dom(1))/2;
norder = f.n;
ndim = numel(dom(1:end/2));

% convert to vol_tree_build output format
iper = 0;
npbox = norder^ndim;
grid = zeros(ndim,norder^2);
nbmax = 2^ndim*f.id(end); % is 2^ndim enough, since some box needs to be refined more than once...
nlmax = f.height(1); % is this ok?
centers0 = 1/2*(f.domain(1:2:end,:)+f.domain(2:2:end,:));
centers = zeros(ndim,nbmax); centers(:,1:size(centers0,2)) = centers0;
nlevels = f.height(1);
nboxes = f.id(end);
boxsize = 2*scale./2.^(0:nlmax); 
laddr = zeros(2,nlmax+1); 
for ilev=0:nlmax
  laddr(1,ilev+1) = sum(f.level<ilev)+1;
  laddr(2,ilev+1) = sum(f.level<=ilev);
end
ilevel = zeros(nbmax,1); ilevel(1:nboxes) = f.level;
iparent = zeros(nbmax,1); iparent(1:nboxes) = f.parent;
nchild = zeros(nbmax,1); nchild(1:nboxes) = sum(f.children>0);
ichild = zeros(2^ndim,nbmax); ichild(:,1:nboxes) = f.children; ichild(ichild==0) = -1; % -1 convention in fortran tree

% this does not belong to fortran tree, but needed for treefun plot
coeffs = cell(1,nbmax); coeffs(1:nboxes) = f.coeffs(1:nboxes);

% call computecoll to get vol_tree_fix_lr coll info
nnbors0 = zeros(nboxes,1); nbors0 = zeros(3^ndim,nboxes);
[nnbors0,nbors0] = ...
  computecoll(ndim,nlevels,nboxes,laddr,boxsize,...
                  centers(:,1:nboxes),iparent(1:nboxes),nchild(1:nboxes),ichild(:,1:nboxes),iper,...
                  nnbors0,nbors0);
nnbors = zeros(nbmax,1); nbors = zeros(3^ndim,nbmax);
nnbors(1:nboxes) = nnbors0; nbors(:,1:nboxes) = nbors0;

% now fix lr (with coeffs, which is different from)
[centers,nboxes,laddr,ilevel,iparent,nchild,ichild,nnbors,nbors,coeffs] = ...
      vol_tree_fix_lr_treefun(ndim,iper,nd,dpars,zpars,ipars,...
                               norder,npbox,grid,centers,nlevels,nboxes,boxsize,...
                               nbmax,nlmax,laddr,ilevel,iparent,nchild,ichild,nnbors,nbors,...
                               coeffs);

centers = centers(:,1:nboxes);
ilevel = ilevel(1:nboxes);
iparent = iparent(1:nboxes);
nchild = nchild(1:nboxes);
ichild = ichild(:,1:nboxes);
nnbors = nnbors(1:nboxes);
nbors = nbors(:,1:nboxes);
coeffs = coeffs(1:nboxes);

% convert back to treefun format, try to follow treefun naming convention
g = f;
scl                 = boxsize(ilevel+1); % x,y,z all the same
g.domain            = zeros(2*ndim,nboxes);
g.domain(1:2:end,:) = centers - 1/2*scl;
g.domain(2:2:end,:) = centers + 1/2*scl;
g.level             = ilevel(:)';
g.id                = 1:nboxes;
g.parent            = iparent(:)';
g.children          = ichild;
g.height            = zeros(1,nboxes);
g.height(nchild>0)  = 1;
for k = length(g.id):-1:1
  if ( ~(g.height(k)==0) )
    g.height(k) = 1 + max(g.height(g.children(:,k)));
  end
end
% redo some stuff
g.col               = uint64(zeros(1,nboxes));
g.row               = uint64(zeros(1,nboxes));
g.dep               = uint64(zeros(1,nboxes));
% g.morton            = uint64(zeros(1,nboxes));
mc = 2^ndim;
for ilev = 0:nlevels-1
  for ibox = laddr(1,ilev+1):laddr(2,ilev+1)
    if nchild(ibox) > 0
      if ndim == 3 % column row stuff
        j              = 1;
        jbox           = g.children(j,ibox);
        g.col(jbox)    = 2*g.col(ibox);
        g.row(jbox)    = 2*g.row(ibox);
        g.dep(jbox)    = 2*g.dep(ibox);
        % g.morton(jbox) = cartesian2morton(g.col(jbox), g.row(jbox));
        j              = 2;
        jbox           = g.children(j,ibox);
        g.col(jbox)    = 2*g.col(ibox) + 1;
        g.row(jbox)    = 2*g.row(ibox);
        g.dep(jbox)    = 2*g.dep(ibox);
        % g.morton(jbox) = cartesian2morton(g.col(jbox), g.row(jbox));
        j              = 3;
        jbox           = g.children(j,ibox);
        g.col(jbox)    = 2*g.col(ibox);
        g.row(jbox)    = 2*g.row(ibox) + 1;
        g.dep(jbox)    = 2*g.dep(ibox);
        % g.morton(jbox) = cartesian2morton(g.col(jbox), g.row(jbox));
        j              = 4;
        jbox           = g.children(j,ibox);
        g.col(jbox)    = 2*g.col(ibox) + 1;
        g.row(jbox)    = 2*g.row(ibox) + 1;
        g.dep(jbox)    = 2*g.dep(ibox);
        % g.morton(jbox) = cartesian2morton(g.col(jbox), g.row(jbox));
        j              = 5;
        jbox           = g.children(j,ibox);
        g.col(jbox)    = 2*g.col(ibox);
        g.row(jbox)    = 2*g.row(ibox);
        g.dep(jbox)    = 2*g.dep(ibox) + 1;
        
        j              = 6;
        jbox           = g.children(j,ibox);
        g.col(jbox)    = 2*g.col(ibox) + 1;
        g.row(jbox)    = 2*g.row(ibox);
        g.dep(jbox)    = 2*g.dep(ibox) + 1;
        
        j              = 7;
        jbox           = g.children(j,ibox);
        g.col(jbox)    = 2*g.col(ibox);
        g.row(jbox)    = 2*g.row(ibox) + 1;
        g.dep(jbox)    = 2*g.dep(ibox) + 1;
        
        j              = 8;
        jbox           = g.children(j,ibox);
        g.col(jbox)    = 2*g.col(ibox) + 1;
        g.row(jbox)    = 2*g.row(ibox) + 1;
        g.dep(jbox)    = 2*g.dep(ibox) + 1;
        
      end
    end
  end
end
g.coeffs = coeffs;
end

function [centers,nboxes,laddr,ilevel,iparent,nchild,ichild,nnbors,nbors,coeffs] = vol_tree_fix_lr_treefun(ndim,iperiod,nd,dpars,zpars,ipars, ...
                                                                                    norder,npbox,grid,centers,nlevels,nboxes,boxsize, ...
                                                                                    nbmax,nlmax,laddr,ilevel,iparent,nchild,ichild,nnbors,nbors, ...
                                                                                    coeffs)

laddrtail=zeros(2,nlmax+1);
% boxsize(0:nlmax),laddr(2,0:nlmax),laddrtail(2,0:nlmax)
idx1 = 1;
%
bs0=boxsize(0+idx1);
%
mc=2^ndim;
mnbors=3^ndim;

iflag = zeros(nbmax,1); iflag2 = zeros(nbmax,1);

for i=1:nboxes
   iflag(i) = 0;
end

% Flag boxes that violate level restriction by "1" Violatioin refers to any
% box that is directly touching a box that is more than one level finer
% Method:
% 1) Carry out upward pass. For each box B, look at the colleagues of B's
% grandparent 
% 2) See if any of those colleagues are childless and in contact with B.
%  
% Note that we only need to get up to level two, as we will not find a
% violation at level 0 and level 1 For such boxes, we set iflag(i) = 1
% figure(1), clf,
for ilev=nlevels:-1:2
  distest = 1.05d0*(boxsize(ilev-1+idx1) + boxsize(ilev-2+idx1))/2.0d0;
  for ibox = laddr(1,ilev+idx1):laddr(2,ilev+idx1) % look at current box
    idad = iparent(ibox); % its parent box
    igranddad = iparent(idad); % its grand parent box
    
    for i=1:nnbors(igranddad) % look at neighbor of grand parent
       jbox = nbors(i,igranddad);
       if((nchild(jbox)==0) && (iflag(jbox)==0)) % if neighbor of grand parent has no child, and has not been flagged
          ict = 0; % count overlapping dimension between ibox & its grand parent's neighbor
          for k=1:ndim
             dis = centers(k,jbox) - centers(k,idad);
             if (iperiod==1)
                dp1=dis+bs0;
                dm1=dis-bs0;
                if (abs(dis).gt.abs(dp1)), dis=dp1; end
                if (abs(dis).gt.abs(dm1)), dis=dm1; end
             end
             if(abs(dis)<=distest), ict = ict + 1; end
          end
          if(ict==ndim)
            iflag(jbox) = 1;
            % figure(1),
            % plotibox(centers(:,ibox),boxsize(ilev-1+idx1)/2);
            % plotibox(centers(:,jbox),4*boxsize(ilev-1+idx1)/2);
            % axis equal
          end
       end
    end
  end
end

% Find all boxes that need to be given a flag+ A flag+ box will be denoted
% by setting iflag(box) = 2 This refers to any box that is not already
% flagged and is bigger than and is contacting a flagged box or another box
% that has already been given a flag +. It is found by performing an upward
% pass and looking at the flagged box's parents colleagues and a flag+
% box's parents colleagues and seeing if they are childless and present the
% case where a bigger box is contacting a flagged or flag+ box.
% figure(2), clf,
%% if secondary violator, 1st send it downstairs during Step 4, then wait to be processed in Step 5 or no need to process at all, seems to be later, from the "worst-case" example 
for ilev = nlevels:-1:1
  distest = 1.05d0*(boxsize(ilev+idx1) + boxsize(ilev-1+idx1))/2.0d0;
  for ibox = laddr(1,ilev+idx1):laddr(2,ilev+idx1)
    if((iflag(ibox)==1) || (iflag(ibox)==2)) % look at current box
      idad = iparent(ibox); % its parent box
      for i=1:nnbors(idad) % look at neighbor of parent
        jbox = nbors(i,idad);
        if((nchild(jbox)==0) && (iflag(jbox)==0)) % if neighbor of parent has no child, and has not been flagged
          ict = 0; % count overlapping dimension between ibox & its parent's neighbor
          for k=1:ndim
            dis = centers(k,jbox) - centers(k,ibox);
            if (iperiod==1)
              dp1=dis+bs0;
              dm1=dis-bs0;
              if (abs(dis)>abs(dp1)), dis=dp1; end
              if (abs(dis)>abs(dm1)), dis=dm1; end
            end
            if(abs(dis)<=distest), ict = ict + 1; end
          end
          if(ict==ndim)
            iflag(jbox) = 2;
            % figure(2),
            % plotibox(centers(:,ibox),boxsize(ilev-1+idx1)/2);
            % plotibox(centers(:,jbox),4*boxsize(ilev-1+idx1)/2);
            % axis equal
          end
        end
      end
    end
  end
end

% Subdivide all flag and flag+ boxes. Flag all the children of flagged
% boxes as flag++. Flag++ boxes are denoted by setting iflag(box) = 3. The
% flag++ boxes need to be checked later to see which of them need further
% refinement. While creating new boxes, we will need to update all the tree
% structures as well. Note that all the flagged boxes live between levels 1
% and nlevels - 2. We process the boxes via a downward pass. We first
% determine the number of boxes that are going to be subdivided at each
% level and everything else accordingly
for ilev = 0:nlevels
   laddrtail(1,ilev+idx1) = 0;
   laddrtail(2,ilev+idx1) = -1;
end

nboxes0 = nboxes;
%% Frank's thesis, Appendix A, Step 4, divide all of the primary and secondary violators once. I think
% figure out iflag = 3 case
for ilev = 1:nlevels-2

  laddrtail(1,ilev+1+idx1) = nboxes+1;

  nbloc = laddr(2,ilev+idx1)-laddr(1,ilev+idx1)+1;
  
  [iflag,centers,nboxes,ilevel,iparent,nchild,ichild,coeffs]  = vol_tree_refine_boxes_flag_treefun(ndim,iflag,nd,npbox,...
                                   dpars,zpars,ipars,grid,nbmax,laddr(1,ilev+idx1),nbloc, ...
                                   centers,boxsize(ilev+1+idx1),nboxes,ilev,ilevel,iparent,nchild,ichild, ...
                                   coeffs);
  
  laddrtail(2,ilev+1+idx1) = nboxes;

end

ishuffle1 = zeros(nboxes,1);
centerstmp = centers(:,1:nboxes); ileveltmp = ilevel(1:nboxes); iparenttmp = iparent(1:nboxes); nchildtmp = nchild(1:nboxes); ichildtmp = ichild(:,1:nboxes); iflagtmp = iflag(1:nboxes);
coeffstmp = coeffs(1:nboxes);
[centerstmp,laddr,ileveltmp,ishuffle1,iparenttmp,nchildtmp,ichildtmp,coeffstmp,iflagtmp] = ... 
              vol_tree_reorg_treefun(ndim,nboxes,nd,npbox,centerstmp,nlevels,laddr,...
                                      laddrtail,ileveltmp,ishuffle1,iparenttmp,nchildtmp,ichildtmp,coeffstmp,iflagtmp);
centers(:,1:nboxes) = centerstmp; ilevel(1:nboxes) = ileveltmp; iparent(1:nboxes) = iparenttmp; nchild(1:nboxes) = nchildtmp; ichild(:,1:nboxes) = ichildtmp; iflag(1:nboxes) = iflagtmp;
coeffs(1:nboxes) = coeffstmp;

nnborstmp = nnbors(1:nboxes); nborstmp = nbors(:,1:nboxes);
[nnborstmp,nborstmp] = computecoll(ndim,nlevels,nboxes,laddr,boxsize,...
                                      centers(:,1:nboxes),iparent(1:nboxes),nchild(1:nboxes),...
                                      ichild(:,1:nboxes),iperiod,nnborstmp,nborstmp);
nnbors(1:nboxes) = nnborstmp; nbors(:,1:nboxes) = nborstmp;


% Processing of flag and flag+ boxes is done Start processing flag++ boxes.
% We will use a similar strategy as before. We keep checking the flag++
% boxes that require subdivision if they still violate the level
% restriction criterion, create the new boxes, append them to the end of
% the list to begin with and in the end reorganize the tree structure. We
% shall accomplish this via a downward pass as new boxes that get added in
% the downward pass will also be processed simultaneously. We shall
% additionally also need to keep on updating the colleague information as
% we proceed in the downward pass
% 
% Reset the flags array to remove all the flag and flag+ cases. This is to
% ensure reusability of the subdivide_flag routine to handle the flag++ case

for ibox=1:nboxes
  if(iflag(ibox)~=3), iflag(ibox) = 0; end
end

for ilev = 0:nlevels
  laddrtail(1,ilev+idx1) = 0;
  laddrtail(2,ilev+idx1) = -1;
end

% boxsize(0:nlmax),laddr(2,0:nlmax),laddrtail(2,0:nlmax)
%% Frank's thesis, Appendix A, Step 5, at each level, for a descendant of a primary violator, test whether it is still in violation
for ilev = 2:nlevels-2

  % Step 1
  iflagtmp = iflag(1:nboxes);
  iflagtmp = vol_updateflags(ndim,iperiod,ilev,nboxes,nlevels, ...
            laddr,nchild(1:nboxes),ichild(:,1:nboxes),nnbors(1:nboxes),nbors(:,1:nboxes),centers(:,1:nboxes),boxsize,iflagtmp);
  iflag(1:nboxes) = iflagtmp;

  iflagtmp = iflag(1:nboxes);
  iflagtmp = vol_updateflags(ndim,iperiod,ilev,nboxes,nlevels, ...
            laddrtail,nchild(1:nboxes),ichild(:,1:nboxes),nnbors(1:nboxes),nbors(:,1:nboxes),centers(:,1:nboxes),boxsize,iflagtmp);
  iflag(1:nboxes) = iflagtmp;

  % Step 2
  laddrtail(1,ilev+1+idx1) = nboxes + 1;

  nbloc = laddr(2,ilev+idx1)-laddr(1,ilev+idx1)+1;
  [iflag,centers,nboxes,ilevel,iparent,nchild,ichild,coeffs] = vol_tree_refine_boxes_flag_treefun(ndim,iflag,nd,npbox, ...
                                                       dpars,zpars,ipars,grid,nbmax,laddr(1,ilev+idx1),nbloc, ...
                                                       centers,boxsize(ilev+1+idx1),nboxes,ilev,ilevel,iparent,nchild,ichild, ...
                                                       coeffs);
  
  nbloc = laddrtail(2,ilev+idx1)-laddrtail(1,ilev+idx1)+1;
  [iflag,centers,nboxes,ilevel,iparent,nchild,ichild,coeffs] = vol_tree_refine_boxes_flag_treefun(ndim,iflag,nd,npbox, ...
                                                       dpars,zpars,ipars,grid,nbmax,laddrtail(1,ilev+idx1),nbloc, ...
                                                       centers,boxsize(ilev+1+idx1),nboxes,ilev,ilevel,iparent,nchild,ichild, ...
                                                       coeffs);
  laddrtail(2,ilev+1+idx1) = nboxes;      

  % Step 3
  for ibox = laddrtail(1,ilev+1+idx1):laddrtail(2,ilev+1+idx1)
    nnbors(ibox) = 0;
    idad = iparent(ibox);
    for i=1:nnbors(idad)
      jbox = nbors(i,idad);
      for j=1:mc
        kbox = ichild(j,jbox);
        if(kbox>0) 
          ifnbor=1;
          for k=1:ndim
            dis=centers(k,kbox)-centers(k,ibox);
            if (iperiod==1) 
              dp1=dis+bs0;
              dm1=dis-bs0;
              if (abs(dis)>abs(dp1)), dis=dp1; end
              if (abs(dis)>abs(dm1)), dis=dm1; end
            end
            if((abs(dis)>1.05*boxsize(ilev+1+idx1)))
              ifnbor=0;
              break
            end
          end
          if (ifnbor==1)
            nnbors(ibox) = nnbors(ibox)+1;
            nbors(nnbors(ibox),ibox) = kbox;
          end
        end
      end
    end
    % End of computing colleagues of box i
  end
  
end

ishuffle2 = zeros(nboxes,1);
centerstmp = centers(:,1:nboxes); ileveltmp = ilevel(1:nboxes); iparenttmp = iparent(1:nboxes); nchildtmp = nchild(1:nboxes); ichildtmp = ichild(:,1:nboxes); iflagtmp = iflag(1:nboxes);
coeffstmp = coeffs(1:nboxes);
[centerstmp,laddr,ileveltmp,ishuffle2,iparenttmp,nchildtmp,ichildtmp,coeffstmp,iflagtmp] = ... 
              vol_tree_reorg_treefun(ndim,nboxes,nd,npbox,centerstmp,nlevels,laddr,...
                                      laddrtail,ileveltmp,ishuffle2,iparenttmp,nchildtmp,ichildtmp,coeffstmp,iflagtmp);
centers(:,1:nboxes) = centerstmp; ilevel(1:nboxes) = ileveltmp; iparent(1:nboxes) = iparenttmp; nchild(1:nboxes) = nchildtmp; ichild(:,1:nboxes) = ichildtmp; iflag(1:nboxes) = iflagtmp;
coeffs(1:nboxes) = coeffstmp;

nnborstmp = nnbors(1:nboxes); nborstmp = nbors(:,1:nboxes);
[nnborstmp,nborstmp] = computecoll(ndim,nlevels,nboxes,laddr,boxsize,...
                                      centers(:,1:nboxes),iparent(1:nboxes),nchild(1:nboxes),...
                                      ichild(:,1:nboxes),iperiod,nnborstmp,nborstmp);
nnbors(1:nboxes) = nnborstmp; nbors(:,1:nboxes) = nborstmp;

end

function [iflag,centers,nbctr,ilevel,iparent,nchild,ichild,coeffs] = vol_tree_refine_boxes_flag_treefun(ndim,iflag,nd,npbox,dpars,zpars,ipars,grid,nboxes,ifirstbox,nbloc,centers,bs,nbctr,nlctr,ilevel,iparent,nchild,ichild,coeffs)

% myfun = @(x,y,z) sin(x).*cos(y).*sin(z)+exp(x.*y.*z);
% myvals = myfun(2*(xx0-1/2),2*(yy0-1/2),2*(zz0-1/2));
% mycoeffs = treefun3.vals2coeffs(myvals);
% mycvals1_2 = myfun(cxx1,cyy1,czz1);
% mycvals1 = reshape(treefun3.coeffs2checkvals(mycoeffs,cxx1(:),cyy1(:),czz1(:)),[norder norder norder]);
% mycvals2 = reshape(treefun3.coeffs2checkvals(mycoeffs,cxx2(:),cyy2(:),czz2(:)),[norder norder norder]);
% mycvals3 = reshape(treefun3.coeffs2checkvals(mycoeffs,cxx3(:),cyy3(:),czz3(:)),[norder norder norder]);
% mycvals4 = reshape(treefun3.coeffs2checkvals(mycoeffs,cxx4(:),cyy4(:),cyy4(:)),[norder norder norder]);
% mycvals5 = reshape(treefun3.coeffs2checkvals(mycoeffs,cxx5(:),cyy5(:),czz5(:)),[norder norder norder]);
% mycvals6 = reshape(treefun3.coeffs2checkvals(mycoeffs,cxx6(:),cyy6(:),czz6(:)),[norder norder norder]);
% mycvals7 = reshape(treefun3.coeffs2checkvals(mycoeffs,cxx7(:),cyy7(:),czz7(:)),[norder norder norder]);
% mycvals8 = reshape(treefun3.coeffs2checkvals(mycoeffs,cxx8(:),cyy8(:),cyy8(:)),[norder norder norder]);
% 
% figure(2), clf, scatter3(2*(xx0(:)-1/2),2*(yy0(:)-1/2),2*(zz0(:)-1/2),[],log10(abs(vals(:))))
% hold on, scatter3(cxx1(:),cyy1(:),czz1(:),[],log10(abs(cvals1(:))))
% hold on, scatter3(cxx2(:),cyy2(:),czz2(:),[],log10(abs(cvals2(:))))
% hold on, scatter3(cxx3(:),cyy3(:),czz3(:),[],log10(abs(cvals3(:))))
% hold on, scatter3(cxx4(:),cyy4(:),czz4(:),[],log10(abs(cvals4(:))))
% hold on, scatter3(cxx5(:),cyy5(:),czz5(:),[],log10(abs(cvals5(:))))
% hold on, scatter3(cxx6(:),cyy6(:),czz6(:),[],log10(abs(cvals6(:))))
% hold on, scatter3(cxx7(:),cyy7(:),czz7(:),[],log10(abs(cvals7(:))))
% hold on, scatter3(cxx8(:),cyy8(:),czz8(:),[],log10(abs(cvals8(:))))

% additional
persistent xx0 yy0 xxx0 yyy0 zz0 zzz0 nstored
norder=round((npbox)^(1/ndim));
n=norder;

if ndim == 3
  nalias = n;
  nrefpts = n; % Sample at equispaced points to test error

  if ( isempty(xx0) || isempty(xxx0) || n ~= nstored )
      nstored = n;
      x0 = (1-cos(pi*(2*(1:nalias)'-1)/(2*nalias)))/2;
      [xx0, yy0, zz0] = ndgrid(x0);
      [xxx0, yyy0, zzz0] = ndgrid(linspace(0, 1, nrefpts));
  end
  sclx = 2*bs;
  scly = 2*bs;
  sclz = 2*bs;
end

%
isgn=zeros(ndim,2^ndim);
isgn=get_child_box_sign(ndim,isgn);

mc=2^ndim;

ilastbox = ifirstbox+nbloc-1;

bsh = bs/2;

isum=zeros(nbloc,1);

isum = cumsum(iflag(ifirstbox:(ifirstbox+nbloc-1))>0); % there are iflag = 2 entries

for i = 1:nbloc
  ibox = ifirstbox + i-1;
  if(iflag(ibox)>0)
    nbl = nbctr + (isum(i)-1)*mc;
    nchild(ibox) = mc;
    % additional
    if ndim == 3
      xx = sclx*xx0 + (centers(1,ibox)-1/2*sclx);
      yy = scly*yy0 + (centers(2,ibox)-1/2*scly);
      zz = sclz*zz0 + (centers(3,ibox)-1/2*sclz);
      vals = treefun3.coeffs2vals(coeffs{ibox});
      % split into eight child boxes
      dom = [-1 1 -1 1 -1 1];
      xmid = mean(dom(1:2));
      ymid = mean(dom(3:4));
      zmid = mean(dom(5:6));
      j=1;
      jbox = nbl+j;
      cdom1              = [dom(1) xmid dom(3) ymid dom(5) zmid];
      csclx1             = diff(cdom1(1:2));
      cscly1             = diff(cdom1(3:4));
      csclz1             = diff(cdom1(5:6));
      cxx1               = csclx1*xx0 + cdom1(1); 
      cyy1               = cscly1*yy0 + cdom1(3); 
      czz1               = csclz1*zz0 + cdom1(5); 
      cvals1             = reshape(treefun3.coeffs2checkvals(coeffs{ibox},cxx1(:),cyy1(:),czz1(:)),[norder norder norder]);
      ccoeffs1           = treefun3.vals2coeffs(cvals1);
      coeffs{jbox}       = ccoeffs1;
      j=2;
      jbox = nbl+j;
      cdom2              = [xmid dom(2) dom(3) ymid dom(5) zmid];
      csclx2             = diff(cdom2(1:2));
      cscly2             = diff(cdom2(3:4));
      csclz2             = diff(cdom2(5:6));
      cxx2               = csclx2*xx0 + cdom2(1); 
      cyy2               = cscly2*yy0 + cdom2(3); 
      czz2               = csclz2*zz0 + cdom2(5); 
      cvals2             = reshape(treefun3.coeffs2checkvals(coeffs{ibox},cxx2(:),cyy2(:),czz2(:)),[norder norder norder]);
      ccoeffs2           = treefun3.vals2coeffs(cvals2);
      coeffs{jbox}       = ccoeffs2;
      j=3;
      jbox = nbl+j;
      cdom3              = [dom(1) xmid ymid dom(4) dom(5) zmid];
      csclx3             = diff(cdom3(1:2));
      cscly3             = diff(cdom3(3:4));
      csclz3             = diff(cdom3(5:6));
      cxx3               = csclx3*xx0 + cdom3(1); 
      cyy3               = cscly3*yy0 + cdom3(3); 
      czz3               = csclz3*zz0 + cdom3(5); 
      cvals3             = reshape(treefun3.coeffs2checkvals(coeffs{ibox},cxx3(:),cyy3(:),czz3(:)),[norder norder norder]);
      ccoeffs3           = treefun3.vals2coeffs(cvals3);
      coeffs{jbox}       = ccoeffs3;
      j=4;
      jbox = nbl+j;
      cdom4              = [xmid dom(2) ymid dom(4) dom(5) zmid];
      csclx4             = diff(cdom4(1:2));
      cscly4             = diff(cdom4(3:4));
      csclz4             = diff(cdom4(5:6));
      cxx4               = csclx4*xx0 + cdom4(1); 
      cyy4               = cscly4*yy0 + cdom4(3); 
      czz4               = csclz4*zz0 + cdom4(5); 
      cvals4             = reshape(treefun3.coeffs2checkvals(coeffs{ibox},cxx4(:),cyy4(:),czz4(:)),[norder norder norder]);
      ccoeffs4           = treefun3.vals2coeffs(cvals4);
      coeffs{jbox}       = ccoeffs4;
      j=5;
      jbox = nbl+j;
      cdom5              = [dom(1) xmid dom(3) ymid zmid dom(6)];
      csclx5             = diff(cdom5(1:2));
      cscly5             = diff(cdom5(3:4));
      csclz5             = diff(cdom5(5:6));
      cxx5               = csclx5*xx0 + cdom5(1); 
      cyy5               = cscly5*yy0 + cdom5(3); 
      czz5               = csclz5*zz0 + cdom5(5); 
      cvals5             = reshape(treefun3.coeffs2checkvals(coeffs{ibox},cxx5(:),cyy5(:),czz5(:)),[norder norder norder]);
      ccoeffs5           = treefun3.vals2coeffs(cvals5);
      coeffs{jbox}       = ccoeffs5;
      j=6;
      jbox = nbl+j;
      cdom6              = [xmid dom(2) dom(3) ymid zmid dom(6)];
      csclx6             = diff(cdom6(1:2));
      cscly6             = diff(cdom6(3:4));
      csclz6             = diff(cdom6(5:6));
      cxx6               = csclx6*xx0 + cdom6(1); 
      cyy6               = cscly6*yy0 + cdom6(3); 
      czz6               = csclz6*zz0 + cdom6(5); 
      cvals6             = reshape(treefun3.coeffs2checkvals(coeffs{ibox},cxx6(:),cyy6(:),czz6(:)),[norder norder norder]);
      ccoeffs6           = treefun3.vals2coeffs(cvals6);
      coeffs{jbox}       = ccoeffs6;
      j=7;
      jbox = nbl+j;
      cdom7              = [dom(1) xmid ymid dom(4) zmid dom(6)];
      csclx7             = diff(cdom7(1:2));
      cscly7             = diff(cdom7(3:4));
      csclz7             = diff(cdom7(5:6));
      cxx7               = csclx7*xx0 + cdom7(1); 
      cyy7               = cscly7*yy0 + cdom7(3); 
      czz7               = csclz7*zz0 + cdom7(5); 
      cvals7             = reshape(treefun3.coeffs2checkvals(coeffs{ibox},cxx7(:),cyy7(:),czz7(:)),[norder norder norder]);
      ccoeffs7           = treefun3.vals2coeffs(cvals7);
      coeffs{jbox}       = ccoeffs7;
      j=8;
      jbox = nbl+j;
      cdom8              = [xmid dom(2) ymid dom(4) zmid dom(6)];
      csclx8             = diff(cdom8(1:2));
      cscly8             = diff(cdom8(3:4));
      csclz8             = diff(cdom8(5:6));
      cxx8               = csclx8*xx0 + cdom8(1); 
      cyy8               = cscly8*yy0 + cdom8(3); 
      czz8               = csclz8*zz0 + cdom8(5); 
      cvals8             = reshape(treefun3.coeffs2checkvals(coeffs{ibox},cxx8(:),cyy8(:),czz8(:)),[norder norder norder]);
      ccoeffs8           = treefun3.vals2coeffs(cvals8);
      coeffs{jbox}       = ccoeffs8;
      %
      coeffs{ibox}       = [];
    end
    for j=1:mc
      jbox = nbl+j;
      for k=1:ndim
        centers(k,jbox) = centers(k,ibox)+isgn(k,j)*bsh;
      end
      iparent(jbox) = ibox;
      nchild(jbox) = 0;
      for l=1:mc
        ichild(l,jbox) = -1;
      end
      ichild(j,ibox) = jbox;
      ilevel(jbox) = nlctr+1; 
      if(iflag(ibox)==1), iflag(jbox) = 3; end
      if(iflag(ibox)==2), iflag(jbox) = 0; end
    end
  end
end

if(nbloc>0), nbctr = nbctr + isum(nbloc)*mc; end

end

function iflag = vol_updateflags(ndim,iperiod,curlev,nboxes,nlevels,laddr,nchild,ichild,nnbors,nbors,centers,boxsize,iflag)
%
%

idx1 = 1;
bs0=boxsize(0+idx1);

mc=2^ndim;
distest = 1.05d0*(boxsize(curlev+idx1) + boxsize(curlev+1+idx1))/2.0d0;

for ibox = laddr(1,curlev+idx1):laddr(2,curlev+idx1)
  if(iflag(ibox)==3)
    iflag(ibox) = 0;
    for i=1:nnbors(ibox)
      jbox = nbors(i,ibox);
      
      for j=1:mc
        kbox = ichild(j,jbox);
        if(kbox>0)
          if(nchild(kbox)>0)
            ict = 0;
            for k=1:ndim
              dis = centers(k,kbox) - centers(k,ibox);
              if (iperiod==1) 
                dp1=dis+bs0;
                dm1=dis-bs0;
                if (abs(dis)>abs(dp1)), dis=dp1; end
                if (abs(dis)>abs(dm1)), dis=dm1; end
              end
              if(abs(dis)<=distest), ict = ict + 1; end
            end
            if(ict==ndim)
               iflag(ibox) = 1;
               break;
            end
          end
        end
      end
      if iflag(ibox) == 1
       break;
      end
    end
    % 1111       continue   
  end
end

end

function [nnbors,nbors] = computecoll(ndim,nlevels,nboxes,laddr,boxsize,centers,iparent,nchild,ichild,iperiod,nnbors,nbors)

% laddr(2,0:nlevels)
% boxsize(0:nlevels)

idx1 = 1;

bs0=boxsize(0+idx1);

mc= 2^ndim;
mnbors=3^ndim;

for i=1:nboxes
  nnbors(i) = 0;
  for j=1:mnbors
    nbors(j,i) = -1;
  end
end

nnbors(1) = 1;
nbors(1,1) = 1;
for ilev = 1:nlevels
  ifirstbox = laddr(1,ilev+idx1);
  ilastbox = laddr(2,ilev+idx1);
  
  for ibox = ifirstbox:ilastbox
    dad = iparent(ibox);
    for i=1:nnbors(dad)
      jbox = nbors(i,dad);
      for j=1:mc
        kbox = ichild(j,jbox);
        if(kbox>0)
          ifnbor=1;
          for k=1:ndim
            dis=abs(centers(k,kbox)-centers(k,ibox));
            if (iperiod==1)
              dp1=bs0-dis;
              if (dp1<dis), dis=dp1; end
            end
            if (dis>1.05*boxsize(ilev+idx1))
              ifnbor=0;
              break
            end
          end
             
          if(ifnbor==1)
            nnbors(ibox) = nnbors(ibox)+1;
            nbors(nnbors(ibox),ibox) = kbox;
          end
        end
      end
    end
  end
end

end

function [centers,laddr,ilevel,ishuffle,iparent,nchild,ichild,coeffs,iflag] = vol_tree_reorg_treefun(ndim,nboxes,nd,npbox,centers,nlevels,laddr,laddrtail,ilevel,ishuffle,iparent,nchild,ichild,coeffs,iflag)
%
%

% laddrtail(2,0:nlevels)
% laddr(2,0:nlevels)
% tladdr(2,0:nlevels)

tladdr = zeros(2,nlevels+1);

idx1 = 1;

mc=2^ndim;

tilevel=zeros(nboxes,1); tiparent=zeros(nboxes,1); tnchild=zeros(nboxes,1);
tichild=zeros(mc,nboxes); tiflag=zeros(nboxes,1); iboxtocurbox=zeros(nboxes,1);
tcenters=zeros(ndim,nboxes);

for ilev = 0:nlevels
  tladdr(1,ilev+idx1) = laddr(1,ilev+idx1);
  tladdr(2,ilev+idx1) = laddr(2,ilev+idx1);
end

tcenters = centers;
tilevel = ilevel;
tiparent = iparent;
tnchild = nchild;
tichild = ichild;
tcoeffs = coeffs;

for ibox=1:nboxes
  tiflag(ibox) = iflag(ibox);
end

for ilev = 0:1
  for ibox = laddr(1,ilev+idx1):laddr(2,ilev+idx1)
    iboxtocurbox(ibox) = ibox;
  end
end

ilevptr=zeros(nlevels+1,1); ilevptr2=zeros(nlevels,1);

ilevptr(2) = laddr(1,2+idx1);

for ilev=2:nlevels
  nblev = laddr(2,ilev+idx1)-laddr(1,ilev+idx1)+1; % original num of boxes at ilev
  ilevptr2(ilev) = ilevptr(ilev) + nblev;
  nblev = laddrtail(2,ilev+idx1)-laddrtail(1,ilev+idx1)+1; % additional number of boxes at ilev due to refinement
  ilevptr(ilev+1) = ilevptr2(ilev) + nblev;
end

curbox = laddr(1,2+idx1);
ishuffle = zeros(nboxes,1);
for ilev=2:nlevels
  laddr(1,ilev+idx1) = curbox;
  for ibox = tladdr(1,ilev+idx1):tladdr(2,ilev+idx1)
    ilevel(curbox) = tilevel(ibox);
    nchild(curbox) = tnchild(ibox);
    ishuffle(ibox) = curbox; % ibox info copied to curbox
    for k=1:ndim
      centers(k,curbox) = tcenters(k,ibox);
    end
    coeffs{curbox} = tcoeffs{ibox};
    iflag(curbox) = tiflag(ibox);
    iboxtocurbox(ibox) = curbox;
    
    curbox = curbox + 1;
  end
  for ibox = laddrtail(1,ilev+idx1):laddrtail(2,ilev+idx1)
    ilevel(curbox) = tilevel(ibox);
    ishuffle(ibox) = curbox; % ibox info copied to curbox
    for k=1:ndim
      centers(k,curbox) = tcenters(k,ibox);
    end
    coeffs{curbox} = tcoeffs{ibox};
    nchild(curbox) = tnchild(ibox);
    iflag(curbox) = tiflag(ibox);
    iboxtocurbox(ibox) = curbox;
    
    curbox = curbox + 1;
  end
  laddr(2,ilev+idx1) = curbox-1;
end

for ibox=1:nboxes
  if(tiparent(ibox)==-1), iparent(iboxtocurbox(ibox)) = -1; end
  if(tiparent(ibox)>0)
    iparent(iboxtocurbox(ibox)) = iboxtocurbox(tiparent(ibox)); end
  for i=1:mc
    if(tichild(i,ibox)==-1), ichild(i,iboxtocurbox(ibox)) = -1; end
    if(tichild(i,ibox)>0)
      ichild(i,iboxtocurbox(ibox)) = iboxtocurbox(tichild(i,ibox)); end
  end
end

end

function isgn = get_child_box_sign(ndim,isgn)

mc = 2^ndim;
for j=1:ndim
  isgn(j,1)=-1;
end

for j=1:ndim
  for i=1:2^(j-1):mc
    if (i>1), isgn(j,i)=-isgn(j,i-2^(j-1)); end
    for k=1:2^(j-1)-1
      isgn(j,i+k)=isgn(j,i);
    end
  end
end

end