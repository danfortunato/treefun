function f = refineBox(f, id)

% Split into four child boxes
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
f.coeffs{cid1}     = [];
% f.col(cid1)        = 2*f.col(id);
% f.row(cid1)        = 2*f.row(id);
% f.morton(cid1)     = cartesian2morton(f.col(cid1), f.row(cid1));

cid2 = length(f.id)+1;
f.domain(:,cid2)   = [xmid dom(2) dom(3) ymid dom(5) zmid];
f.id(cid2)         = cid2;
f.parent(cid2)     = id;
f.children(:,cid2) = 0;
f.level(cid2)      = f.level(id)+1;
f.height(cid2)     = 0;
f.coeffs{cid2}     = [];
% f.col(cid2)        = 2*f.col(id) + 1;
% f.row(cid2)        = 2*f.row(id);
% f.morton(cid2)     = cartesian2morton(f.col(cid2), f.row(cid2));

cid3 = length(f.id)+1;
f.domain(:,cid3)   = [dom(1) xmid ymid dom(4) dom(5) zmid];
f.id(cid3)         = cid3;
f.parent(cid3)     = id;
f.children(:,cid3) = 0;
f.level(cid3)      = f.level(id)+1;
f.height(cid3)     = 0;
f.coeffs{cid3}     = [];
% f.col(cid3)        = 2*f.col(id);
% f.row(cid3)        = 2*f.row(id) + 1;
% f.morton(cid3)     = cartesian2morton(f.col(cid3), f.row(cid3));

cid4 = length(f.id)+1;
f.domain(:,cid4)   = [xmid dom(2) ymid dom(4) dom(5) zmid];
f.id(cid4)         = cid4;
f.parent(cid4)     = id;
f.children(:,cid4) = 0;
f.level(cid4)      = f.level(id)+1;
f.height(cid4)     = 0;
f.coeffs{cid4}     = [];
% f.col(cid4)        = 2*f.col(id) + 1;
% f.row(cid4)        = 2*f.row(id) + 1;
% f.morton(cid4)     = cartesian2morton(f.col(cid4), f.row(cid4));

cid5 = length(f.id)+1;
f.domain(:,cid5)   = [dom(1) xmid dom(3) ymid zmid dom(6)];
f.id(cid5)         = cid5;
f.parent(cid5)     = id;
f.children(:,cid5) = 0;
f.level(cid5)      = f.level(id)+1;
f.height(cid5)     = 0;
f.coeffs{cid5}     = [];

cid6 = length(f.id)+1;
f.domain(:,cid6)   = [xmid dom(2) dom(3) ymid zmid dom(6)];
f.id(cid6)         = cid6;
f.parent(cid6)     = id;
f.children(:,cid6) = 0;
f.level(cid6)      = f.level(id)+1;
f.height(cid6)     = 0;
f.coeffs{cid6}     = [];

cid7 = length(f.id)+1;
f.domain(:,cid7)   = [dom(1) xmid ymid dom(4) zmid dom(6)];
f.id(cid7)         = cid7;
f.parent(cid7)     = id;
f.children(:,cid7) = 0;
f.level(cid7)      = f.level(id)+1;
f.height(cid7)     = 0;
f.coeffs{cid7}     = [];

cid8 = length(f.id)+1;
f.domain(:,cid8)   = [xmid dom(2) ymid dom(4) zmid dom(6)];
f.id(cid8)         = cid8;
f.parent(cid8)     = id;
f.children(:,cid8) = 0;
f.level(cid8)      = f.level(id)+1;
f.height(cid8)     = 0;
f.coeffs{cid8}     = [];

f.children(:,id) = [cid1 cid2 cid3 cid4 cid5 cid6 cid7 cid8];
f.height(id) = 1;
f.coeffs{id} = [];

end
