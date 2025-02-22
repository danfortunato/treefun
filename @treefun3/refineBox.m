function f = refineBox(f, id, func)

persistent xx0 yy0 zz0 ww0 nstored
nalias = f.n;
nd = numel(func);
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

% Split into eight child boxes
dom = f.domain(:,id);
xmid = mean(dom(1:2));
ymid = mean(dom(3:4));
zmid = mean(dom(5:6));

cid1 = length(f.id)+1;
f.domain(:,cid1)   = [dom(1) xmid dom(3) ymid dom(5) zmid];
f.id(cid1)         = cid1;
f.parent(cid1)     = id;
f.children(:,cid1) = 0;
f.level(cid1)      = f.level(id)+1;
f.height(cid1)     = 0;
cdom1              = f.domain(:,cid1);
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
cvals1             = cat(4,cvals1{:}); % additional for nd
ccoeffs1           = treefun3.vals2coeffs(cvals1);
f.coeffs{cid1}     = ccoeffs1(1:f.n,1:f.n,1:f.n,:); % to replace f.coeffs{cid1} = [];
f.rint(:,cid1)     = squeeze(sqrt((csclx1*cscly1*csclz1)*sum(cvals1.^2.*ww0, [1 2 3])));
f.vmax(:,cid1)     = squeeze(max(abs(cvals1),[],[1 2 3]));
f.col(cid1)        = 2*f.col(id);
f.row(cid1)        = 2*f.row(id);
f.dep(cid1)        = 2*f.dep(id);
% f.morton(cid1)     = cartesian2morton(f.col(cid1), f.row(cid1));

cid2 = length(f.id)+1;
f.domain(:,cid2)   = [xmid dom(2) dom(3) ymid dom(5) zmid];
f.id(cid2)         = cid2;
f.parent(cid2)     = id;
f.children(:,cid2) = 0;
f.level(cid2)      = f.level(id)+1;
f.height(cid2)     = 0;
cdom2              = f.domain(:,cid2);
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
cvals2             = cat(4,cvals2{:}); % additional for nd
ccoeffs2           = treefun3.vals2coeffs(cvals2);
f.coeffs{cid2}     = ccoeffs2(1:f.n,1:f.n,1:f.n,:); 
f.rint(:,cid2)     = squeeze(sqrt((csclx2*cscly2*csclz2)*sum(cvals2.^2.*ww0, [1 2 3])));
f.vmax(:,cid2)     = squeeze(max(abs(cvals2),[],[1 2 3]));
f.col(cid2)        = 2*f.col(id) + 1;
f.row(cid2)        = 2*f.row(id);
f.dep(cid2)        = 2*f.dep(id);
% f.morton(cid2)     = cartesian2morton(f.col(cid2), f.row(cid2));

cid3 = length(f.id)+1;
f.domain(:,cid3)   = [dom(1) xmid ymid dom(4) dom(5) zmid];
f.id(cid3)         = cid3;
f.parent(cid3)     = id;
f.children(:,cid3) = 0;
f.level(cid3)      = f.level(id)+1;
f.height(cid3)     = 0;
cdom3              = f.domain(:,cid3);
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
cvals3             = cat(4,cvals3{:}); % additional for nd
ccoeffs3           = treefun3.vals2coeffs(cvals3);
f.coeffs{cid3}     = ccoeffs3(1:f.n,1:f.n,1:f.n,:); 
f.rint(:,cid3)     = squeeze(sqrt((csclx3*cscly3*csclz3)*sum(cvals3.^2.*ww0, [1 2 3])));
f.vmax(:,cid3)     = squeeze(max(abs(cvals3),[],[1 2 3]));
f.col(cid3)        = 2*f.col(id);
f.row(cid3)        = 2*f.row(id) + 1;
f.dep(cid3)        = 2*f.dep(id);
% f.morton(cid3)     = cartesian2morton(f.col(cid3), f.row(cid3));

cid4 = length(f.id)+1;
f.domain(:,cid4)   = [xmid dom(2) ymid dom(4) dom(5) zmid];
f.id(cid4)         = cid4;
f.parent(cid4)     = id;
f.children(:,cid4) = 0;
f.level(cid4)      = f.level(id)+1;
f.height(cid4)     = 0;
cdom4              = f.domain(:,cid4);
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
cvals4             = cat(4,cvals4{:}); % additional for nd
ccoeffs4           = treefun3.vals2coeffs(cvals4);
f.coeffs{cid4}     = ccoeffs4(1:f.n,1:f.n,1:f.n,:); 
f.rint(:,cid4)     = squeeze(sqrt((csclx4*cscly4*csclz4)*sum(cvals4.^2.*ww0, [1 2 3])));
f.vmax(:,cid4)     = squeeze(max(abs(cvals4),[],[1 2 3]));
f.col(cid4)        = 2*f.col(id) + 1;
f.row(cid4)        = 2*f.row(id) + 1;
f.dep(cid4)        = 2*f.dep(id);
% f.morton(cid4)     = cartesian2morton(f.col(cid4), f.row(cid4));

cid5 = length(f.id)+1;
f.domain(:,cid5)   = [dom(1) xmid dom(3) ymid zmid dom(6)];
f.id(cid5)         = cid5;
f.parent(cid5)     = id;
f.children(:,cid5) = 0;
f.level(cid5)      = f.level(id)+1;
f.height(cid5)     = 0;
cdom5              = f.domain(:,cid5);
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
cvals5             = cat(4,cvals5{:}); % additional for nd
ccoeffs5           = treefun3.vals2coeffs(cvals5);
f.coeffs{cid5}     = ccoeffs5(1:f.n,1:f.n,1:f.n,:); 
f.rint(:,cid5)     = squeeze(sqrt((csclx5*cscly5*csclz5)*sum(cvals5.^2.*ww0, [1 2 3])));
f.vmax(:,cid5)     = squeeze(max(abs(cvals5),[],[1 2 3]));
f.col(cid5)        = 2*f.col(id);
f.row(cid5)        = 2*f.row(id);
f.dep(cid5)        = 2*f.dep(id) + 1;

cid6 = length(f.id)+1;
f.domain(:,cid6)   = [xmid dom(2) dom(3) ymid zmid dom(6)];
f.id(cid6)         = cid6;
f.parent(cid6)     = id;
f.children(:,cid6) = 0;
f.level(cid6)      = f.level(id)+1;
f.height(cid6)     = 0;
cdom6              = f.domain(:,cid6);
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
cvals6             = cat(4,cvals6{:}); % additional for nd
ccoeffs6           = treefun3.vals2coeffs(cvals6);
f.coeffs{cid6}     = ccoeffs6(1:f.n,1:f.n,1:f.n,:); 
f.rint(:,cid6)     = squeeze(sqrt((csclx6*cscly6*csclz6)*sum(cvals6.^2.*ww0, [1 2 3])));
f.vmax(:,cid6)     = squeeze(max(abs(cvals6),[],[1 2 3]));
f.col(cid6)        = 2*f.col(id) + 1;
f.row(cid6)        = 2*f.row(id);
f.dep(cid6)        = 2*f.dep(id) + 1;

cid7 = length(f.id)+1;
f.domain(:,cid7)   = [dom(1) xmid ymid dom(4) zmid dom(6)];
f.id(cid7)         = cid7;
f.parent(cid7)     = id;
f.children(:,cid7) = 0;
f.level(cid7)      = f.level(id)+1;
f.height(cid7)     = 0;
cdom7              = f.domain(:,cid7);
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
cvals7             = cat(4,cvals7{:}); % additional for nd
ccoeffs7           = treefun3.vals2coeffs(cvals7);
f.coeffs{cid7}     = ccoeffs7(1:f.n,1:f.n,1:f.n,:); 
f.rint(:,cid7)     = squeeze(sqrt((csclx7*cscly7*csclz7)*sum(cvals7.^2.*ww0, [1 2 3])));
f.vmax(:,cid7)     = squeeze(max(abs(cvals7),[],[1 2 3]));
f.col(cid7)        = 2*f.col(id);
f.row(cid7)        = 2*f.row(id) + 1;
f.dep(cid7)        = 2*f.dep(id) + 1;

cid8 = length(f.id)+1;
f.domain(:,cid8)   = [xmid dom(2) ymid dom(4) zmid dom(6)];
f.id(cid8)         = cid8;
f.parent(cid8)     = id;
f.children(:,cid8) = 0;
f.level(cid8)      = f.level(id)+1;
f.height(cid8)     = 0;
cdom8              = f.domain(:,cid8);
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
cvals8             = cat(4,cvals8{:}); % additional for nd
ccoeffs8           = treefun3.vals2coeffs(cvals8);
f.coeffs{cid8}     = ccoeffs8(1:f.n,1:f.n,1:f.n,:); 
f.rint(:,cid8)     = squeeze(sqrt((csclx8*cscly8*csclz8)*sum(cvals8.^2.*ww0, [1 2 3])));
f.vmax(:,cid8)     = squeeze(max(abs(cvals8),[],[1 2 3]));
f.col(cid8)        = 2*f.col(id) + 1;
f.row(cid8)        = 2*f.row(id) + 1;
f.dep(cid8)        = 2*f.dep(id) + 1;

f.children(:,id) = [cid1 cid2 cid3 cid4 cid5 cid6 cid7 cid8];
f.height(id) = 1;
f.coeffs{id} = [];

end