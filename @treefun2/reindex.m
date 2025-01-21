function f = reindex(f)
%REINDEX   Reindex a TREEFUN2.

% Inverse permutation
nboxes = length(f.id);
iperm(f.id) = 1:nboxes;

% Now reindex
f.id = 1:nboxes;
f.root = iperm(f.root);

idx = f.parent ~= 0;
f.parent(idx) = iperm(f.parent(idx));

idx = f.children ~= 0;
f.children(idx) = iperm(f.children(idx));

idx = ~isnan(f.flatNeighbors);
f.flatNeighbors(idx) = iperm(f.flatNeighbors(idx));

f.leafNeighbors = cellfun(@(x) iperm(x), f.leafNeighbors, 'UniformOutput', false);


end
