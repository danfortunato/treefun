function f = refineUnresolved(f, tol)
%REFINEUNRESOLVED   Refine the unresolved boxes of a TREEFUN2.

if ( nargin < 2 )
    tol = 1e-12;
end

ids = leaves(f);

for id = ids(:).'

    if ( isResolved(f.coeffs{id}, tol) )
        % Don't split
        continue
    end

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

f = balance(f);
[f.flatNeighbors, f.leafNeighbors] = generateNeighbors(f);

end

function [resolved, coeffs] = isResolved(coeffs, tol)

n = size(coeffs, 1);
Ex = sum(abs(coeffs(end-1:end,:)), 'all') / (2*n);
Ey = sum(abs(coeffs(:,end-1:end)), 'all') / (2*n);
err = (Ex + Ey) / 2;
vals = treefun2.coeffs2vals(coeffs);
vmax = max(abs(vals(:)));
resolved = ( err < tol * max(vmax, 1) );

end
