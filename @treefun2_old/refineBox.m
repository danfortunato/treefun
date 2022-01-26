function f = refineBox(f, id)

% Split into four child boxes
dom = f.boxes(id).domain;
xmid = mean(dom(1:2));
ymid = mean(dom(3:4));
parent = f.boxes(id);

child1 = struct();
child1.domain   = [dom(1) xmid dom(3) ymid];
child1.id       = length(f.boxes)+1;
child1.parent   = id;
child1.children = [];
child1.level    = parent.level+1;
child1.height   = 0;
child1.coeffs   = [];
child1.col      = 2*(parent.col-1) + 1;
child1.row      = 2*(parent.row-1) + 1;
f.boxes(end+1) = child1;

child2 = struct();
child2.domain   = [xmid dom(2) dom(3) ymid];
child2.id       = length(f.boxes)+1;
child2.parent   = id;
child2.children = [];
child2.level    = parent.level+1;
child2.height   = 0;
child2.coeffs   = [];
child2.col      = 2*(parent.col-1) + 2;
child2.row      = 2*(parent.row-1) + 1;
f.boxes(end+1) = child2;

child3 = struct();
child3.domain   = [dom(1) xmid ymid dom(4)];
child3.id       = length(f.boxes)+1;
child3.parent   = id;
child3.children = [];
child3.level    = parent.level+1;
child3.height   = 0;
child3.coeffs   = [];
child3.col      = 2*(parent.col-1) + 1;
child3.row      = 2*(parent.row-1) + 2;
f.boxes(end+1) = child3;

child4 = struct();
child4.domain   = [xmid dom(2) ymid dom(4)];
child4.id       = length(f.boxes)+1;
child4.parent   = id;
child4.children = [];
child4.level    = parent.level+1;
child4.height   = 0;
child4.coeffs   = [];
child4.col      = 2*(parent.col-1) + 2;
child4.row      = 2*(parent.row-1) + 2;
f.boxes(end+1) = child4;

f.boxes(id).children = [child1.id, child2.id, ...
                        child3.id, child4.id];

%f.boxes(id).children = [child3.id, child4.id, ...
%                        child2.id, child1.id];

f.boxes(id).height = 1;
f.boxes(id).coeffs = [];

end
