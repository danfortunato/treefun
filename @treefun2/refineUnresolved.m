function f = refineUnresolved(f, tol)
%REFINEUNRESOLVED   Refine the unresolved boxes of a TREEFUN2.

if ( nargin < 2 )
    tol = 1e-12;
end

ids = leaves(f);
unresolved = [];
for id = ids(:).'
    if ( ~isResolved(f.coeffs{id}, tol) )
        unresolved(end+1) = id; %#ok<AGROW>
    end
end
f = refine(f, unresolved);

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
