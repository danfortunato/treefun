function f = refineBox(f, id)

% Split into four child boxes
dom = f.domain(:,id);
xmid = mean(dom(1:2));
ymid = mean(dom(3:4));

cid1 = length(f.id)+1;
f.domain(:,cid1)   = [dom(1) xmid dom(3) ymid];
f.id(cid1)         = cid1;
f.parent(cid1)     = id;
f.children(:,cid1) = 0;
f.level(cid1)      = f.level(id)+1;
f.height(cid1)     = 0;
f.coeffs{cid1}     = [];
f.col(cid1)        = 2*(f.col(id)-1) + 1;
f.row(cid1)        = 2*(f.row(id)-1) + 1;

cid2 = length(f.id)+1;
f.domain(:,cid2)   = [xmid dom(2) dom(3) ymid];
f.id(cid2)         = cid2;
f.parent(cid2)     = id;
f.children(:,cid2) = 0;
f.level(cid2)      = f.level(id)+1;
f.height(cid2)     = 0;
f.coeffs{cid2}     = [];
f.col(cid2)        = 2*(f.col(id)-1) + 2;
f.row(cid2)        = 2*(f.row(id)-1) + 1;

cid3 = length(f.id)+1;
f.domain(:,cid3)   = [dom(1) xmid ymid dom(4)];
f.id(cid3)         = cid3;
f.parent(cid3)     = id;
f.children(:,cid3) = 0;
f.level(cid3)      = f.level(id)+1;
f.height(cid3)     = 0;
f.coeffs{cid3}     = [];
f.col(cid3)        = 2*(f.col(id)-1) + 1;
f.row(cid3)        = 2*(f.row(id)-1) + 2;

cid4 = length(f.id)+1;
f.domain(:,cid4)   = [xmid dom(2) ymid dom(4)];
f.id(cid4)         = cid4;
f.parent(cid4)     = id;
f.children(:,cid4) = 0;
f.level(cid4)      = f.level(id)+1;
f.height(cid4)     = 0;
f.coeffs{cid4}     = [];
f.col(cid4)        = 2*(f.col(id)-1) + 2;
f.row(cid4)        = 2*(f.row(id)-1) + 2;

f.children(:,id) = [cid1 cid2 cid3 cid4];
f.height(id) = 1;
f.coeffs{id} = [];

end
