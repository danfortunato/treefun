function levels = levels(f)
%LEVELS   Get the boxes from a TREEFUN2 in level order.
%   LEVELS(F) returns the box IDs from the TREEFUN2 F in level order. The
%   result is a cell array whose K-th entry is an array of box IDs from
%   level K in the tree.

nlevels = height(f);
levels = cell(nlevels, 1);
for l = 0:nlevels-1
    ids = [];
    for k = 1:length(f.id)
        if ( f.level(k) == l )
            ids(end+1) = f.id(k); %#ok<AGROW>
        end
    end
    levels{l+1} = ids;
end

end
