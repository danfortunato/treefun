function f = refine(f)
%REFINE   Refine a TREEFUN2.

ids = leaves(f);
for id = ids(:).'
    % This was a leaf, so we'll use its coeffs to evaluate on
    % the new children
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

f = balance2(f);

end
