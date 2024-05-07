function g = balance(f)
% based on @treefun2.balance
% initial attempt in BBBBBBoston

if ( isempty(f) )
    g = f;
    return
end

% Note: We assume IDs are stored in level order

% Unpack data because MATLAB seems to be slow for class member access...
fmorton  = f.morton(:);
fcol     = f.col;
frow     = f.row;
fdep     = f.dep;
fxy      = [fcol(:) frow(:) fdep(:)];
flevel   = f.level;
fheight  = f.height;
fcoeffs  = f.coeffs;
fnlevels = max(fheight) + 1;
fnboxes  = length(f.id);

% Generate level indices for f
flevelIdx = zeros(fnlevels+1, 1);
idx = 1;
flevelIdx(idx) = 1;
currentLevel = 0;
for k = 2:fnboxes
    if ( currentLevel ~= flevel(k) )
        idx = idx+1;
        flevelIdx(idx) = k;
        currentLevel = flevel(k);
    end
end
flevelIdx(end) = fnboxes+1;

% Generate a set of 2:1-balanced Morton IDs level by level from the bottom
% (leaves) to the top (root)
S = [];
T = cell(fnlevels, 1);
for l = fnlevels-1:-1:1
    idx = flevelIdx(l):flevelIdx(l+1)-1;
    fmortonl = fmorton(idx);
    N = unique([S ; fmortonl(fheight(idx) > 0)]);
    S = mortonparent(mortoncolleagues(N, l));
    [child1, child2, child3, child4, child5, child6, child7, child8] = mortonchildren(N);
    T{l+1} = unique([child1 ; child2 ; child3 ; child4 ; ...
                     child5 ; child6 ; child7 ; child8]);

    % fxyl = fxy(idx,:);
    % N = unique([S ; fxyl(fheight(idx)>0, :)], 'rows', 'stable');
    % S = xyparent(xycolleagues(N, l));
    % [childxy1, childxy2, childxy3, childxy4] = xychildren(N);
    % %childxy = [childxy1 childxy2 childxy3 childxy4].';
    % cxy = zeros(4*size(childxy1,1), 2, 'uint64');
    % cxy(1:4:end,:) = childxy1;
    % cxy(2:4:end,:) = childxy2;
    % cxy(3:4:end,:) = childxy3;
    % cxy(4:4:end,:) = childxy4;
    % T{l+1} = unique(cxy, 'rows', 'stable');
    % T{l+1} = mortonsort(T{l+1});
end
T{1} = uint64(0);
%T{1} = uint64([0 0]);

% Generate level indices for g
gnlevels = fnlevels;
gnboxes = 0;
glevelIdx = zeros(gnlevels+1, 1);
glevelIdx(1) = 1;
for l = 1:gnlevels
    gnboxes = gnboxes + size(T{l}, 1);
    glevelIdx(l+1) = gnboxes + 1;
end

% First, make the tree from the Morton IDs
gmorton = cell2mat(T).';
[gcol, grow, gdep] = morton2cartesian(gmorton);
%gmorton = [];
%xy = cell2mat(T);
%gcol = xy(:,1).';
%grow = xy(:,2).';
glevel    = zeros(1, gnboxes);
gheight   = zeros(1, gnboxes);
gparent   = zeros(1, gnboxes);
gchildren = zeros(8, gnboxes);
gcoeffs   = cell(1, gnboxes);

for l = gnlevels:-1:2
    id = glevelIdx(l);
    jparent = 1;
    ps = mortonparent(T{l});
    %pxys = xyparent(T{l});
    for k = 1:8:size(T{l},1)
        ids = id:id+7;
        glevel(ids) = l-1;
        p = ps(k);
        %pxy = pxys(k,:);
        while ( T{l-1}(jparent) ~= p && jparent <= length(T{l-1}) )
        %while ( ~all(T{l-1}(jparent,:) == pxy) && jparent <= size(T{l-1}, 1) )
            jparent = jparent + 1;
        end
        if ( jparent > size(T{l-1}, 1) )
            error('Couldn''t find parent node.');
        end
        pid = glevelIdx(l-1) + jparent - 1;
        gparent(ids) = pid;
        gchildren(:,pid) = ids;
        gheight(pid) = 1 + max(gheight(ids));
        id = id + 8;
    end
end

dom = f.domain(:,1);
scl = 1./2.^glevel;
offx = scl*(dom(2)-dom(1));
offy = scl*(dom(4)-dom(3));
offz = scl*(dom(6)-dom(5));
xo = offx.*double(gcol);
yo = offy.*double(grow);
zo = offz.*double(gdep);
gdomain = [dom(1) + xo; dom(1)+offx + xo; ...
           dom(3) + yo; dom(3)+offy + yo; ...
           dom(5) + zo; dom(5)+offz + zo];

% Finally, fill in the coefficients
gid = 1;
for fid = 1:fnboxes
    if ( fheight(fid) == 0 )
        fm = fmorton(fid);
        %fc = fcol(fid);
        %fr = frow(fid);
        fl = flevel(fid);
        while ( ~(gmorton(gid)==fm && glevel(gid)==fl) && gid <= gnboxes )
        %while ( ~(gcol(gid)==fc && grow(gid)==fr && glevel(gid)==fl) && gid <= gnboxes )
            gid = gid + 1;
        end
        if ( gid > gnboxes )
            error('Couldn''t find corresponding tree node.');
        end
        gcoeffs{gid} = fcoeffs{fid};
    end
end

% Propagate the coefficients to the leaves, since they may be stored in
% intermediate nodes of g (which used to be leaves of f)
for l = 1:gnlevels-1
    for gid = glevelIdx(l):glevelIdx(l+1)-1
        if ( gheight(gid) > 0 && ~isempty(gcoeffs{gid}) )
            [LLD, LRD, ULD, URD, LLT, LRT, ULT, URT] = coeffs2children(gcoeffs{gid});
            gcoeffs{gchildren(1,gid)} = LLD;
            gcoeffs{gchildren(2,gid)} = LRD;
            gcoeffs{gchildren(3,gid)} = ULD;
            gcoeffs{gchildren(4,gid)} = URD;
            gcoeffs{gchildren(5,gid)} = LLT;
            gcoeffs{gchildren(6,gid)} = LRT;
            gcoeffs{gchildren(7,gid)} = ULT;
            gcoeffs{gchildren(8,gid)} = URT;
        end
    end
end

% Make a new treefun2
g = treefun3();
g.n        = f.n;
g.domain   = gdomain;
g.level    = glevel;
g.height   = gheight;
g.id       = 1:gnboxes;
g.parent   = gparent;
g.children = gchildren;
g.coeffs   = gcoeffs;
g.col      = gcol;
g.row      = grow;
g.morton   = gmorton;

end

%% Column/row routines
function nbrxy = xycolleagues(xy, level)
    x = xy(:,1);
    y = xy(:,2);
    z = xy(:,3);
    %bnd = 2^(level-1)-1;
    bnd = bitshift(uint64(1), level-1) - 1;
    % Order is important here:
    %nbrxy = max(uint64(0), min(bnd, [x-1 y ; x+1 y ; x y-1 ; x y+1 ; x-1 y-1 ; x+1 y-1 ; x-1 y+1 ; x+1 y+1]));
    nbrxy = max(uint64(0), min(bnd, [[x-1 y-1 z-1; x y-1 z-1; x+1 y-1 z-1; x-1 y z-1; x y z-1; x+1 y z-1; x-1 y+1 z-1; x y+1 z-1; x+1 y+1 z-1]; ...
                                     [x-1 y-1 z  ; x y-1 z  ; x+1 y-1 z  ; x-1 y z  ;          x+1 y z  ; x-1 y+1 z  ; x y+1 z  ; x+1 y+1 z  ]; ...
                                     [x-1 y-1 z+1; x y-1 z+1; x+1 y-1 z+1; x-1 y z+1; x y z+1; x+1 y z+1; x-1 y+1 z+1; x y+1 z+1; x+1 y+1 z+1]]));
    %nbrxy = max(uint64(0), min(bnd, [x-1 y-1; x y-1; x+1 y-1 ; x-1 y ; x+1 y ; x-1 y+1 ; x y+1 ; x+1 y+1]));
end

function pxy = xyparent(xy)
    pxy = bitshift(xy, -1);
end

function [cxy1, cxy2, cxy3, cxy4, cxy5, cxy6, cxy7, cxy8] = xychildren(pxy)
    cxy1 = bitshift(pxy, 2);
    cxy2 = cxy1 + uint64([1 0 0]);
    cxy3 = cxy1 + uint64([0 1 0]);
    cxy4 = cxy1 + uint64([1 1 0]);
    cxy5 = cxy1 + uint64([0 0 1]);
    cxy6 = cxy1 + uint64([1 0 1]);
    cxy7 = cxy1 + uint64([0 1 1]);
    cxy8 = cxy1 + uint64([1 1 1]);
end


%% Morton routines
function [xy, idx] = mortonsort(xy)
    m = xy2morton(xy);
    [m, idx] = sort(m);
    xy = morton2xy(m);
end

function morton = xy2morton(xy)
    partxy = Part1By2_64(xy);
    morton = bitshift(partxy(:,3), 2) + bitshift(partxy(:,2), 1) + partxy(:,1);
end

function xy = morton2xy(morton)
    xy = Compact1By2_64([morton bitshift(morton, -1) bitshift(morton, -2)]);
end

function coll = mortoncolleagues(m, level)
    %bnd = 2^(level-1)-1;
    bnd = bitshift(uint64(1), level-1) - 1;
    %[x, y] = morton2cartesian(m);
    %nbrx = max(uint64(0), min(bnd, [x-1 ; x+1 ; x   ; x   ; x-1 ; x+1 ; x-1 ; x+1]));
    %nbry = max(uint64(0), min(bnd, [y   ; y   ; y-1 ; y+1 ; y-1 ; y-1 ; y+1 ; y+1]));
    %coll = cartesian2morton(nbrx, nbry);
    xy = morton2xy(m);
    x = xy(:,1);
    y = xy(:,2);
    z = xy(:,3);
    nbrxy = min(bnd, [[x-1 y-1 z-1; x y-1 z-1; x+1 y-1 z-1; x-1 y z-1; x y z-1; x+1 y z-1; x-1 y+1 z-1; x y+1 z-1; x+1 y+1 z-1]; ...
                      [x-1 y-1 z  ; x y-1 z  ; x+1 y-1 z  ; x-1 y z  ;          x+1 y z  ; x-1 y+1 z  ; x y+1 z  ; x+1 y+1 z  ]; ...
                      [x-1 y-1 z+1; x y-1 z+1; x+1 y-1 z+1; x-1 y z+1; x y z+1; x+1 y z+1; x-1 y+1 z+1; x y+1 z+1; x+1 y+1 z+1]]);
    coll = xy2morton(nbrxy);
end

function p = mortonparent(c)
    p = bitshift(c, -3);
end

function [c1, c2, c3, c4, c5, c6, c7, c8] = mortonchildren(p)
    c1 = bitshift(p, 3);
    c2 = c1 + 1;
    c3 = c1 + 2;
    c4 = c1 + 3;
    c5 = c1 + 4;
    c6 = c1 + 5;
    c7 = c1 + 6;
    c8 = c1 + 7;
end

function [x, y, z] = morton2cartesian(morton)
    x = Compact1By2_64(morton);
    y = Compact1By2_64(bitshift(morton, -1));
    z = Compact1By2_64(bitshift(morton, -2));
end

function morton = cartesian2morton(x, y)
    morton = bitshift(Part1By2_64(y), 1) + Part1By2_64(x);
end

function x = Part1By1(x)
% "Insert" a 0 bit after each of the 16 low bits of x
    x = bitand(x, 0x0000ffff);
    x = bitand(bitxor(x, bitshift(x, 8)), 0x00ff00ff);
    x = bitand(bitxor(x, bitshift(x, 4)), 0x0f0f0f0f);
    x = bitand(bitxor(x, bitshift(x, 2)), 0x33333333);
    x = bitand(bitxor(x, bitshift(x, 1)), 0x55555555);
end

function x = Compact1By1(x)
% Inverse of Part1By1 - "delete" all odd-indexed bits
    x = bitand(x, 0x55555555);
    x = bitand(bitxor(x, bitshift(x, -1)), 0x33333333);
    x = bitand(bitxor(x, bitshift(x, -2)), 0x0f0f0f0f);
    x = bitand(bitxor(x, bitshift(x, -4)), 0x00ff00ff);
    x = bitand(bitxor(x, bitshift(x, -8)), 0x0000ffff);
end

function x = Part1By2_64(x)
% "Insert" two 0 bits after each of the 16 low bits of x
% https://github.com/trevorprater/pymorton/blob/f248683c19d90e904193c58bbd03e77bc2c43768/pymorton/pymorton.py#L20

x = bitand(x, uint64(0x1fffff));
x = bitand(bitxor(x, bitshift(x, 32)), uint64(0x1f00000000ffff));
x = bitand(bitxor(x, bitshift(x, 16)), uint64(0x1f0000ff0000ff));
x = bitand(bitxor(x, bitshift(x, 8)),  uint64(0x100f00f00f00f00f));
x = bitand(bitxor(x, bitshift(x, 4)),  uint64(0x10c30c30c30c30c3));
x = bitand(bitxor(x, bitshift(x, 2)),  uint64(0x1249249249249249));

end

function x = Compact1By2_64(x)
% Inverse of Part1By2 - "delete" all odd-indexed bits
    x = bitand(x, uint64(0x1249249249249249));
    x = bitand(bitxor(x, bitshift(x, -2)),  uint64(0x10c30c30c30c30c3));
    x = bitand(bitxor(x, bitshift(x, -4)),  uint64(0x100f00f00f00f00f));
    x = bitand(bitxor(x, bitshift(x, -8)),  uint64(0x1f0000ff0000ff));
    x = bitand(bitxor(x, bitshift(x, -16)), uint64(0x1f00000000ffff));
    x = bitand(bitxor(x, bitshift(x, -32)), uint64(0x1fffff));
end