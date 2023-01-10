function f = refine(f, ids)
%REFINE   Refine a TREEFUN2.

if ( nargin < 2 )
    ids = leaves(f);
end

for id = ids(:).'
    % This was a leaf, so we'll use its coeffs to evaluate on the new
    % children
    coeffs = f.coeffs{id};
    % Split into four child boxes
    f = refineBox(f, id);
    children = f.children(:,id);
    [LL, LR, UL, UR] = coeffs2children(coeffs);
    f.coeffs{children(1)} = LL; % Lower left
    f.coeffs{children(2)} = LR; % Lower right
    f.coeffs{children(3)} = UL; % Upper left
    f.coeffs{children(4)} = UR; % Upper right
end

% Do a cumulative sum in reverse to correct the heights
for k = length(f.id):-1:1
    if ( ~isLeaf(f, k) )
        f.height(k) = 1 + max(f.height(f.children(:,k)));
    end
end

[f.flatNeighbors, f.leafNeighbors] = generateNeighbors(f);

end
