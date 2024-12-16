function g = permute(f, perm)
%PERMUTE   Permute a TREEFUN2.

g = f;
g.domain = f.domain(:,perm);
g.level  = f.level(perm);
g.height = f.height(perm);
g.coeffs = f.coeffs(perm);
g.col    = f.col(perm);
g.row    = f.row(perm);
g.morton = f.morton(perm);
g.id     = perm(f.id(perm));
g.root   = perm(f.root);

g.parent = f.parent(perm);
idx = g.parent ~= 0;
g.parent(idx) = perm(g.parent(idx));

g.children = f.children(:,perm);
idx = g.children ~= 0;
g.children(idx) = perm(g.children(idx));

g.flatNeighbors = f.flatNeighbors(:,perm);
idx = ~isnan(g.flatNeighbors);
g.flatNeighbors(idx) = perm(g.flatNeighbors(idx));

g.leafNeighbors = cellfun(@(x) perm(x), f.leafNeighbors(:,perm), 'UniformOutput', false);

end
