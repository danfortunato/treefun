function u = poisson(f, isource)

%%% Point FMM plus local correction

n = f.n;
domain = f.domain;
coeffs = f.coeffs;

% Convert second-kind Chebyshev points to Legendre points on each box
Timer.tic();
ids = leaves(f);
xx     = cell(length(ids), 1);
yy     = cell(length(ids), 1);
ww     = cell(length(ids), 1);
ff     = cell(length(ids), 1);
ff_cfs = cell(length(ids), 1);
[x0, w0] = legpts(n, [0 1]);
[xx0, yy0] = meshgrid(x0);
ww0 = w0(:) * w0(:).';
CC2LV = chebcoeffs2legvals(eye(n));
CC2LC = chebcoeffs2legcoeffs(eye(n));
for k = 1:length(ids)
    id = ids(k);
    dom = domain(:,id);
    sclx = diff(dom(1:2));
    scly = diff(dom(3:4));
    xx{k} = sclx*xx0 + dom(1);
    yy{k} = scly*yy0 + dom(3);
    ww{k} = sclx*scly*ww0;
    ff{k}     = CC2LV * coeffs{id} * CC2LV.';
    ff_cfs{k} = CC2LC * coeffs{id} * CC2LC.';
end
Timer.toc('Vals to coeffs');

if ( nargin < 2 )
    isource = true(size(ids));
end
nisource = length(find(isource));

xx_ = [xx{isource}]; xx_ = xx_(:);
yy_ = [yy{isource}]; yy_ = yy_(:);
ww_ = [ww{isource}]; ww_ = ww_(:);
ff_ = [ff{isource}]; ff_ = ff_(:);

iprec = 4;
nsource = nisource * n^2;
source = [xx_ yy_].';
ifcharge = 1;
charge = -ww_.*ff_/(2*pi);
ifdipole = 0;
dipstr = zeros(1, nsource);
dipvec = zeros(2, nsource);
ifpot  = 1;
ifgrad = 0;
ifhess = 0;
target = [];
ntarget = length(target);
ifpottarg  = 0;
ifgradtarg = 0;
ifhesstarg = 0;

% Call a point FMM
% Note: This will not work if box points overlap (e.g., if we use
% second-kind Chebyshev points)
Timer.tic();
out = rfmm2dpart(iprec,                    ...
                 nsource, source,          ...
                 ifcharge, charge,         ...
                 ifdipole, dipstr, dipvec, ...
                 ifpot, ifgrad, ifhess,    ...
                 ntarget, target,          ...
                 ifpottarg, ifgradtarg, ifhesstarg);
Timer.toc('FMM');

if ( out.ier ~= 0 )
    error('FMM exited with error code %d\n.', out.ier);
end

uu = reshape(out.pot, n, n, nisource);
temp = zeros(n, n, length(ids));
temp(:,:,isource) = uu;
uu = temp;

Timer.tic();
uu = correct_close_batched(f, uu, ff, ff_cfs, ww, isource);
Timer.toc('Close correction');

u = f;
coeffs = u.coeffs;
LV2CC = legvals2chebcoeffs(eye(n));
for k = 1:length(ids)
    id = ids(k);
    coeffs{id} = LV2CC * uu(:,:,k) * LV2CC.';
end
u.coeffs = coeffs;

end

function uu = correct_close(f, uu, ff, ff_cfs, ww)

persistent naiveMats closeMats

if ( isempty(naiveMats) || isempty(closeMats) )
    n = size(ff{1}, 1);
    file = sprintf([treefunroot '/quadrature/laplace_%d.mat'], n);
    load(file, 'naiveMats', 'closeMats');
end

leafIDs = leaves(f);
idToLinearIdx = zeros(max(leafIDs), 1);
for k = 1:length(leafIDs)
    id = leafIDs(k);
    idToLinearIdx(id) = k;
end

for k = 1:length(leafIDs)
    id  = leafIDs(k);
    dom = f.domain(:,id);
    sclx = diff(dom(1:2));
    scly = diff(dom(3:4));

    tau = ff{k}.*ww{k};
    C = ff_cfs{k} * sclx * scly;

    % Correct self
    naive_u = reshape(naiveMats(:,:,1) * tau(:), f.n, f.n);
    naive_u = naive_u - log(sclx)/(2*pi) * (sum(tau, 'all') - tau);
    accurate_u = reshape(reshape(C, [], f.n^2) * closeMats(:,:,1), f.n, f.n);
    accurate_u = accurate_u.' - C(1,1)*log(sclx)/(2*pi);
    uu{k} = uu{k} + accurate_u - naive_u;

    % Correct neighbors
    sourcedom = f.domain(:,id);
    sourcelevel = f.level(id);
    neighborIDs = f.leafNeighbors(:,id);
    for nid = cat(1,neighborIDs{:}).'
        kn = idToLinearIdx(nid);
        neighbordom = f.domain(:,leafIDs(kn));
        neighborlevel = f.level(leafIDs(kn));
        code = computeCode(sourcedom, sourcelevel, neighbordom, neighborlevel);
        naive_u = reshape(naiveMats(:,:,code) * tau(:), f.n, f.n);
        naive_u = naive_u - log(sclx)/(2*pi) * sum(tau, 'all');
        accurate_u = reshape(reshape(C, [], f.n^2) * closeMats(:,:,code), f.n, f.n);
        accurate_u = accurate_u.' - C(1,1)*log(sclx)/(2*pi);
        uu{kn} = uu{kn} + accurate_u - naive_u;
    end
end

end

function uu = correct_close_batched(f, uu, ff, ff_cfs, ww, isource)

persistent naiveMats closeMats

if ( isempty(naiveMats) || isempty(closeMats) )
    n = size(ff{1}, 1);
    file = sprintf([treefunroot '/quadrature/laplace_%d.mat'], n);
    load(file, 'naiveMats', 'closeMats');
end

n = f.n;
domain = f.domain;
level = f.level;
leafNeighbors = f.leafNeighbors;

leafIDs = leaves(f);
nleaf = length(leafIDs);
idToLinearIdx = zeros(max(leafIDs), 1);
for k = 1:nleaf
    idToLinearIdx(leafIDs(k)) = k;
end

tau = zeros(numel(ff{1}), nleaf);
C   = zeros(numel(ff{1}), nleaf);

ncodes = 32 + 1;
naive_u     = cell(ncodes, 1);
accurate_u  = cell(ncodes, 1);
codeToLeaf  = cell(ncodes, 1);
codeToLeaf{1} = 1:nleaf; % All leaves will have a self correction
codeToNeighbor = zeros(ncodes, nleaf);
codeToIdx      = zeros(ncodes, nleaf);

for k = 1:nleaf
    if ( ~isource(k) )
        continue
    end

    id  = leafIDs(k);
    dom = domain(:,id);
    sclx = diff(dom(1:2));
    scly = diff(dom(3:4));

    tau(:,k) = reshape(ff{k}.*ww{k}, [], 1);
    C(:,k) = reshape(ff_cfs{k}*sclx*scly, [], 1);

    sourcedom = domain(:,id);
    sourcelevel = level(id);
    neighborIDs = leafNeighbors(:,id);
    for nid = cat(1,neighborIDs{:}).'
        kn = idToLinearIdx(nid);
        if ( isource(kn) )
            neighbordom = domain(:,leafIDs(kn));
            neighborlevel = level(leafIDs(kn));
            code = computeCode(sourcedom, sourcelevel, neighbordom, neighborlevel);
            codeToLeaf{code}(end+1) = k;
            codeToIdx(code,k) = length(codeToLeaf{code});
            codeToNeighbor(code,k) = kn;
        end
    end
end

% Do the corrections in a batch
for code = 1:ncodes
    naive_u{code}    = reshape(naiveMats(:,:,code) * tau(:,codeToLeaf{code}), n, n, []);
    accurate_u{code} = reshape(closeMats(:,:,code) *   C(:,codeToLeaf{code}), n, n, []);
end

for k = 1:nleaf
    if ( ~isource(k) )
        continue
    end

    id = leafIDs(k);
    dom = domain(:,id);
    sclx = diff(dom(1:2));

    % Correct self
    naive_u{1}(:,:,k) = naive_u{1}(:,:,k) - log(sclx)/(2*pi) * reshape(sum(tau(:,k)) - tau(:,k), n, n);
    accurate_u{1}(:,:,k) = accurate_u{1}(:,:,k) - C(1,k)*log(sclx)/(2*pi);
    uu(:,:,k) = uu(:,:,k) + accurate_u{1}(:,:,k) - naive_u{1}(:,:,k);

    % Correct neighbors
    codes = find(codeToNeighbor(:,k));
    for code = codes(:).'
        idx = codeToIdx(code,k);
        nu = naive_u{code}(:,:,idx);
        au = accurate_u{code}(:,:,idx);
        nu = nu - log(sclx)/(2*pi) * sum(tau(:,k), 'all');
        au = au - C(1,k)*log(sclx)/(2*pi);
        kn = codeToNeighbor(code,k);
        uu(:,:,kn) = uu(:,:,kn) + au - nu;
    end
end

end

function code = computeCode(dom, level, ndom, nlevel)

% dom = domain of source box
% level = level of source box
% ndom = domain of neighbor box
% nlevel = level of neighbor box

% Codes go clockwise starting from the left neighbor
switch level
    case nlevel
        % Neighbor is the same level
        if ( ndom(1) < dom(1) )
            if ( ndom(3) < dom(3) )
                code = 8; % Left down
            elseif ( abs(ndom(3) - dom(3)) < 1e-14 )
                code = 1; % Left
            else
                code = 2; % Left up
            end
        elseif ( abs(ndom(1) - dom(1)) < 1e-14 )
            if ( ndom(3) < dom(3) )
                code = 7; % Down
            else
                code = 3; % Up
            end
        else
            if ( ndom(3) < dom(3) )
                code = 6; % Right down
            elseif ( abs(ndom(3) - dom(3)) < 1e-14 )
                code = 5; % Right
            else
                code = 4; % Right up
            end
        end
    case nlevel+1
        % Neighbor is coarser
        if ( ndom(1) < dom(1) )
            if ( abs(ndom(4) - dom(3)) < 1e-14 )
                if ( abs(ndom(2) - dom(1)) < 1e-14 )
                    code = 29; % Left down
                else
                    code = 24; % Down with rights aligning
                end
            elseif ( abs(ndom(4) - dom(4)) < 1e-14 )
                code = 25; % Left with tops aligning
            elseif ( abs(ndom(3) - dom(3)) < 1e-14 )
                code = 21; % Left with bottoms aligning
            else
                if ( abs(ndom(2) - dom(1)) < 1e-14 )
                    code = 30; % Left up
                else
                    code = 26; % Up with rights aligning
                end
            end
        elseif ( ndom(2) > dom(2) )
            if ( abs(ndom(4) - dom(3)) < 1e-14 )
                if ( abs(ndom(1) - dom(2)) < 1e-14 )
                    code = 32; % Right down
                else
                    code = 28; % Down with lefts aligning
                end
            elseif ( abs(ndom(4) - dom(4)) < 1e-14 )
                code = 23; % Right with tops aligning
            elseif ( abs(ndom(3) - dom(3)) < 1e-14 )
                code = 27; % Right with bottoms aligning
            else
                if ( abs(ndom(1) - dom(2)) < 1e-14 )
                    code = 31; % Right up
                else
                    code = 22; % Up with lefts aligning
                end
            end
        end
    case nlevel-1
        % Neighbor is finer
        if ( ndom(1) < dom(1) )
            if ( ndom(3) < dom(3) )
                code = 20; % Left down
            elseif ( abs(ndom(3) - dom(3)) < 1e-14 )
                code = 9; % Left with bottoms aligning
            elseif ( abs(ndom(4) - dom(4)) < 1e-14 )
                code = 10; % Left with tops aligning
            else
                code = 11; % Left up
            end
        elseif ( abs(ndom(1) - dom(1)) < 1e-14 )
            if ( ndom(3) < dom(3) )
                code = 19; % Down with lefts aligning
            else
                code = 12; % Up with lefts aligning
            end
        elseif ( abs(ndom(2) - dom(2)) < 1e-14 )
            if ( ndom(3) < dom(3) )
                code = 18; % Down with rights aligning
            else
                code = 13; % Up with rights aligning
            end
        else
            if ( ndom(3) < dom(3) )
                code = 17; % Right down
            elseif ( abs(ndom(3) - dom(3)) < 1e-14 )
                code = 16; % Right with bottoms aligning
            elseif ( abs(ndom(4) - dom(4)) < 1e-14 )
                code = 15; % Right with tops aligning
            else
                code = 14; % Right up
            end
        end
    otherwise
        error(['No lookup entry for neighbor found. Is your treefun2 ' ...
            'level-restricted?']);
end

% Add one because "self" is first
code = code + 1;

end
