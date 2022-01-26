function u = poisson(f)

%% Travis's box code

% lvls = levels(f);
% nlev = length(lvls)-1;
% nblevel = cellfun(@length, lvls);
% nboxes = length(f.boxes);
% 
% % TODO: Child box ordering?
% ichildbox  = zeros(4, nboxes);
% iboxlev    = zeros(nboxes, 1);
% istartlev  = [1; cumsum(nblevel(1:end-1))+1];
% 
% levelbox   = [f.boxes.level];   % Intrinsic f.boxes order
% iparentbox = [f.boxes.parent];  % Intrinsic f.boxes order
% icolbox    = [f.boxes.col];     % Intrinsic f.boxes order
% irowbox    = [f.boxes.row];     % Intrinsic f.boxes order
% for k = 1:nboxes
%     if ( f.boxes(k).parent == 0 )
%         iparentbox(k) = -1;
%     else
%         iparentbox(k) = f.boxes(k).parent;
%     end
%     if ( isempty(f.boxes(k).children) )
%         ichildbox(:,k) = -1;
%     else
%         ichildbox(:,k) = f.boxes(k).children([3 4 2 1]);
%         %ichildbox(:,k) = f.boxes(k).children;
%     end
% end
% 
% idx = 1;
% for l = 1:length(lvls)
%     for k = 1:length(lvls{l})
%         box = lvls{l}(k);
%         iboxlev(idx) = box.id;
%         idx = idx + 1;
%     end
% end
% 
% % These are the points that the box code will use to sample the function on
% % each leaf box:
% %
% %    x0 = linspace(0+1/8, 1-1/8, 4);
% %    [xx0, yy0] = meshgrid(x0);
% %
% % The order of values starts from the lower left corner of the box and
% % proceeds upwards row by row from left to right.
% ff = zeros(16, nboxes);
% % leafs = leaves(f);
% % for k = 1:length(leafs)
% %     id = leafs(k).id;
% %     vals = coeffs2boxcodevals(leafs(k).coeffs);
% %     vals = flipud(vals);
% %     vals = vals.';
% %     ff(:,id) = vals(:);
% % end
% 
% for k = 1:nboxes
%     ibox = iboxlev(k);
%     if ( ichildbox(1,ibox) < 0 )
%         coeffs = f.boxes(ibox).coeffs;
%         coeffs = coeffs(1:4,1:4);
%         vals = coeffs2boxcodevals(coeffs);
%         %vals = vals.';
%         %vals = flipud(vals);
%         ff(:,ibox) = vals(:);
%         
% %         dom = f.boxes(ibox).domain;
% %         sclx = diff(dom(1:2));
% %         scly = diff(dom(3:4));
% %         x0 = linspace(dom(1)+sclx/8, dom(2)-sclx/8, 4);
% %         y0 = linspace(dom(3)+scly/8, dom(4)-scly/8, 4);
% %         [xx0, yy0] = meshgrid(x0, y0);
% %         surf(xx0, yy0, vals)
% % %         [xx, yy] = chebpts2(4, 4, dom);
% % %         surf(xx, yy, chebvals)
% %         hold on
% %         drawnow
% %         shg
%     end
% end
% 
% iprec = 3;
% ifgrad = 0;
% ifhess = 0;
% 
% % Domain info
% zll = f.domain([1 3]);         % Lower left corner of domain
% blength = diff(f.domain(1:2)); % Width of domain
% %zll = [-0.5 -0.5];
% %blength = 1;
% 
% %ff = 1+0*ff;
% out = pfmm2d4pw(nlev, levelbox, iparentbox, ichildbox, icolbox, ...
%     irowbox, nboxes, nblevel, iboxlev, istartlev, ifgrad, ifhess, ...
%     iprec, ff, zll, blength);
% 
% u = f;
% % for k = 1:nboxes
% %     if ( ~isempty(u.boxes(k).children) )
% %         %u.boxes(k).children = u.boxes(k).children([3 4 2 1]);
% %         %u.boxes(k).children = u.boxes(k).children([4 3 1 2]);
% %     end
% % end
% 
% % leafs = leaves(f);
% % for k = 1:length(leafs)
% %     id = leafs(k).id;
% %     vals = out.pot(:,id);
% %     if ( max(abs(vals)) > 1e9 )
% %         k
% %     end
% %     vals = reshape(vals, 4, 4);
% %     vals = flipud(vals);
% %     vals = vals.';
% %     coeffs = zeros(f.n, f.n);
% %     coeffs(1:4,1:4) = boxcodevals2coeffs(vals);
% %     u.boxes(id).coeffs = coeffs;
% % end
% 
% normSol = 0;
% 
% for k = 1:nboxes
%     ibox = iboxlev(k);
%     if ( ichildbox(1,ibox) < 0 )
%         vals = out.pot(:,ibox);
%         if ( max(abs(vals)) > 1e9 )
%             ibox
%         end
%         vals = reshape(vals, 4, 4);
%         %vals = flipud(vals);
%         vals = vals.';
%         coeffs = zeros(f.n, f.n);
%         [coeffs(1:4,1:4), chebvals] = boxcodevals2coeffs(vals);
%         u.boxes(ibox).coeffs = coeffs;
%         
%         %err = max(abs(vals(:) - truevals(:)));
%         %normSol = max(err, normSol);
% 
% %         dom = u.boxes(ibox).domain;
% %         sclx = diff(dom(1:2));
% %         scly = diff(dom(3:4));
% %         x0 = linspace(dom(1)+sclx/8, dom(2)-sclx/8, 4);
% %         y0 = linspace(dom(3)+scly/8, dom(4)-scly/8, 4);
% %         [xx0, yy0] = meshgrid(x0, y0);
% %         %surf(xx0, yy0, vals)
% %         [xx, yy] = chebpts2(4, 4, dom);
% %         surf(xx, yy, chebvals)
% %         hold on
% %         drawnow
% %         shg
%     end
% end

%% Point FMM plus local correction

% Convert second-kind Chebyshev points to Legendre points on each box
leaf = leaves(f);
xx = cell(length(leaf), 1);
yy = cell(length(leaf), 1);
ww = cell(length(leaf), 1);
ff = leafvals(f);
[x0, w0] = legpts(f.n, [0 1]);
[xx0, yy0] = meshgrid(x0);
ww0 = w0(:) * w0(:).';
CV2LV = chebvals2legvals(eye(f.n));
for k = 1:length(leaf)
    dom = leaf(k).domain;
    sclx = diff(dom(1:2));
    scly = diff(dom(3:4));
    xx{k} = sclx*xx0 + dom(1);
    yy{k} = scly*yy0 + dom(3);
    ww{k} = sclx*scly*ww0;
    ff{k} = CV2LV * ff{k} * CV2LV.';
end

xx_ = [xx{:}]; xx_ = xx_(:);
yy_ = [yy{:}]; yy_ = yy_(:);
ww_ = [ww{:}]; ww_ = ww_(:);
ff_ = [ff{:}]; ff_ = ff_(:);

iprec = 5;
nsource = length(f) * f.n^2;
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
out = rfmm2dpart(iprec,                    ...
                 nsource, source,          ...
                 ifcharge, charge,         ...
                 ifdipole, dipstr, dipvec, ...
                 ifpot, ifgrad, ifhess,    ...
                 ntarget, target,          ...
                 ifpottarg, ifgradtarg, ifhesstarg);

uu = reshape(out.pot, f.n, f.n, length(leaf));
uu = squeeze(mat2cell(uu, f.n, f.n, ones(length(leaf), 1)));
ff_cfs = cell(length(leaf), 1);
LV2LC = legvals2legcoeffs(eye(f.n));
for k = 1:length(leaf)
    ff_cfs{k} = LV2LC * ff{k} * LV2LC.';
end
%uu = correct_close(f, uu, ff, ff_cfs, ww);
uu = correct_close_batched(f, uu, ff, ff_cfs, ww);

u = f;
LV2CC = legvals2chebcoeffs(eye(f.n));
for k = 1:length(leaf)
    uu{k} = LV2CC * uu{k} * LV2CC.';
    id = leaf(k).id;
    u.boxes(id).coeffs = uu{k};
end

end

function uu = correct_close(f, uu, ff, ff_cfs, ww)

persistent naiveMats closeMats

if ( isempty(naiveMats) || isempty(closeMats) )
    load('quadrature/laplace_16.mat', 'naiveMats', 'closeMats');
end

leaf = leaves(f);
idToLinearIdx = zeros(max([leaf.id]), 1);
for k = 1:length(leaf)
    idToLinearIdx(leaf(k).id) = k;
end

for k = 1:length(leaf)
    id  = leaf(k).id;
    dom = leaf(k).domain;
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
    neighborIDs = f.leafNeighbors(:, id);
    for nid = [neighborIDs{:}]
        kn = idToLinearIdx(nid);
        code = computeCode(leaf(k), leaf(kn));
        naive_u = reshape(naiveMats(:,:,code) * tau(:), f.n, f.n);
        naive_u = naive_u - log(sclx)/(2*pi) * sum(tau, 'all');
        accurate_u = reshape(reshape(C, [], f.n^2) * closeMats(:,:,code), f.n, f.n);
        accurate_u = accurate_u.' - C(1,1)*log(sclx)/(2*pi);
        uu{kn} = uu{kn} + accurate_u - naive_u;
    end
end

end

function uu = correct_close_batched(f, uu, ff, ff_cfs, ww)

persistent naiveMats closeMats
%naiveMats = [];
%closeMats = [];

if ( isempty(naiveMats) || isempty(closeMats) )
    n = size(ff{1}, 1);
    file = sprintf([treefunroot '/quadrature/laplace_%d.mat'], n);
    load(file, 'naiveMats', 'closeMats');
end

leaf = leaves(f);
nleaf = length(leaf);
idToLinearIdx = zeros(max([leaf.id]), 1);
for k = 1:nleaf
    idToLinearIdx(leaf(k).id) = k;
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
    id  = leaf(k).id;
    dom = leaf(k).domain;
    sclx = diff(dom(1:2));
    scly = diff(dom(3:4));

    tau(:,k) = reshape(ff{k}.*ww{k}, [], 1);
    C(:,k) = reshape(ff_cfs{k}*sclx*scly, [], 1);

    neighborIDs = f.leafNeighbors(:,id);
    for nid = [neighborIDs{:}]
        kn = idToLinearIdx(nid);
        code = computeCode(leaf(k), leaf(kn));
        codeToLeaf{code}(end+1) = k;
        codeToIdx(code,k) = length(codeToLeaf{code});
        codeToNeighbor(code,k) = kn;
    end
end

% Do the corrections in a batch
for code = 1:ncodes
    naive_u{code}    = reshape(naiveMats(:,:,code) * tau(:,codeToLeaf{code}), f.n, f.n, []);
    accurate_u{code} = reshape(closeMats(:,:,code) *   C(:,codeToLeaf{code}), f.n, f.n, []);
end

for k = 1:nleaf
    dom = leaf(k).domain;
    sclx = diff(dom(1:2));

    % Correct self
    naive_u{1}(:,:,k) = naive_u{1}(:,:,k) - log(sclx)/(2*pi) * reshape(sum(tau(:,k)) - tau(:,k), f.n, f.n);
    accurate_u{1}(:,:,k) = accurate_u{1}(:,:,k) - C(1,k)*log(sclx)/(2*pi);
    uu{k} = uu{k} + accurate_u{1}(:,:,k) - naive_u{1}(:,:,k);

    % Correct neighbors
    codes = find(codeToNeighbor(:,k));
    for code = codes(:).'
        idx = codeToIdx(code,k);
        nu = naive_u{code}(:,:,idx);
        au = accurate_u{code}(:,:,idx);
        nu = nu - log(sclx)/(2*pi) * sum(tau(:,k), 'all');
        au = au - C(1,k)*log(sclx)/(2*pi);
        kn = codeToNeighbor(code,k);
        uu{kn} = uu{kn} + au - nu;
    end
end

end

function code = computeCode(sourcebox, neighborbox)

dom  = sourcebox.domain;
ndom = neighborbox.domain;

% Codes go clockwise starting from the left neighbor
switch sourcebox.level
    case neighborbox.level
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
    case neighborbox.level+1
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
    case neighborbox.level-1
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
