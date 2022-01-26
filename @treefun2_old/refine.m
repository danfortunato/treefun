function f = refine(f)
%REFINE   Refine a TREEFUN2.

leaf = leaves(f);
for k = 1:length(leaf)
    id = leaf(k).id;
    % This was a leaf, so we'll use its coeffs to evaluate on
    % the new children
    coeffs = f.boxes(id).coeffs;
    % Split into four child boxes
    f = refineBox(f, id);
    children = f.boxes(id).children;
    [LL, LR, UL, UR] = coeffs2children(coeffs);
    f.boxes(children(1)).coeffs = LL; % Lower left
    f.boxes(children(2)).coeffs = LR; % Lower right
    f.boxes(children(3)).coeffs = UL; % Upper left
    f.boxes(children(4)).coeffs = UR; % Upper right
end

f = balance2(f);

end
