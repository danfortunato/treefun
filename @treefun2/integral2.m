function I = integral2(f)
%INTEGRAL2   Double integral of a TREEFUN2.
%   I = INTEGRAL2(F) returns the double integral of the TREEFUN2 F over its
%   domain.
%
%   See also SUM2.

I = 0;

% Empty check:
if ( isempty(f) )
    return
end

% If a patch uses an N x N discretization, then quadrature is performed on
% that patch using N points.

ids = leaves(f);
[~, w0] = chebpts(f.n, [0 1]);
ww0 = w0(:) * w0(:).';
[ny, nx] = size(f.coeffs{ids(1)});
qx = nx;
qy = ny;

for id = ids(:).'
    dom = f.domain(:,id);
    sclx = diff(dom(1:2));
    scly = diff(dom(3:4));
    U = zeros(qy,qx);
    U(1:ny,1:nx) = f.coeffs{id};
    U = U(1:qy,1:qx);
    V = treefun2.coeffs2vals(U);
    I = I + sum(sum(V .* ww0 * sclx * scly));
end

end
