function u = boxcode2d(f, eps)

lvls = levels(f);
nlev = length(lvls)-1;
nblevel = cellfun(@length, lvls);
nboxes = length(f.id);
n = f.n;
npts = n^2;
nd = 2;
pttype = 'F';
polytype = 'T';
zll = f.domain([1 3]);         % Lower left corner of domain
blength = diff(f.domain(1:2)); % Width of domain

% TODO: Child box ordering?
levelbox   = f.level;
iparentbox = f.parent;
ichildbox  = f.children;
icolbox    = f.col;
irowbox    = f.row;
istartlev  = [1; cumsum(nblevel(1:end-1))+1];
iparentbox(iparentbox == 0) = -1;
ichildbox(ichildbox == 0) = -1;
iboxlev = [lvls{:}];

% iboxlev = zeros(nboxes, 1);
% idx = 1;
% for l = 1:length(lvls)
%     for k = 1:length(lvls{l})
%        box = lvls{l}(k);
%        iboxlev(idx) = box.id;
%        idx = idx + 1;
%     end
% end

ifcolleag = false;
ifnbr = false;
iicolleagbox = -1;
ineighbors = -1;
innbrs = -1;

ilevelbox   = 21;
iicolbox    = ilevelbox + nboxes;
iirowbox    = iicolbox + nboxes;
iichildbox  = iirowbox + nboxes;
iiparentbox = iichildbox + 4*nboxes;
iiboxlev    = iiparentbox + nboxes;
inblevel    = iiboxlev + nboxes;
iistartlev  = inblevel + nlev + 1;
ltot        = iistartlev + nlev + 1;

% Pack into itree
%litree = 10*nboxes + 5*(nlev+1) + 100
litree = ltot;
itree = zeros(litree, 1);
itree(1)  = nboxes;
itree(2)  = nlev;
itree(3)  = ifcolleag;
itree(4)  = ifnbr;
itree(5)  = ilevelbox;
itree(6)  = iicolbox;
itree(7)  = iirowbox;
itree(8)  = iichildbox;
itree(9)  = iiparentbox;
itree(10) = iiboxlev;
itree(11) = inblevel;
itree(12) = iistartlev;
itree(13) = iicolleagbox;
itree(14) = ineighbors;
itree(15) = innbrs;

itree(ilevelbox:ilevelbox+nboxes-1)     = levelbox;
itree(iicolbox:iicolbox+nboxes-1)       = icolbox;
itree(iirowbox:iirowbox+nboxes-1)       = irowbox;
itree(iichildbox:iichildbox+4*nboxes-1) = ichildbox;
itree(iiparentbox:iiparentbox+nboxes-1) = iparentbox;
itree(iiboxlev:iiboxlev+nboxes-1)       = iboxlev;
itree(inblevel:inblevel+nlev)           = nblevel;
itree(iistartlev:iistartlev+nlev)       = istartlev;

I = eye(n);
CC2LV = chebcoeffs2legvals(I);
LV2CC = legvals2chebcoeffs(I);

fvals = zeros(npts, nboxes);
leaf = leaves(f);
for id = leaf(:).'
    cheb_cfs = f.coeffs{id};
    leg_vals = CC2LV * cheb_cfs.' * CC2LV.';
    fvals(:,id) = leg_vals(:);
end

ifpot  = true;
ifgrad = false;
ifhess = false;

U = boxcode2d(eps, itree, blength, nd, npts, n, pttype, polytype, fvals, ifpot, ifgrad, ifhess);

pot = U.pot;
u = f;
for id = leaf(:).'
    leg_vals = reshape(pot(:,id), n, n);
    u.coeffs{id} = LV2CC * leg_vals.' * LV2CC.';
end

end
