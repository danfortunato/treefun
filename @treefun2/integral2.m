function I = integral2(f)
%INTEGRAL2   Double integral of an ULTRASEM.SOL.
%   I = INTEGRAL2(F) returns the double integral of the ULTRASEM.SOL F over
%   its domain.
%
%   See also SUM2.

I = 0;

% Empty check:
if ( isempty(f) )
    return
end

% If a patch uses an N x N discretization, then quadrature is performed on
% that patch using N points.

leaf = leaves(f);
[~, w0] = chebpts(f.n, [0 1]);
ww0 = w0(:) * w0(:).';
[ny, nx] = size(leaf(1).coeffs);
qx = nx;
qy = ny;

for k = 1:length(leaf)
    dom = leaf(k).domain;
    sclx = diff(dom(1:2));
    scly = diff(dom(3:4));
    U = zeros(qy,qx);
    U(1:ny,1:nx) = leaf(k).coeffs;
    U = U(1:qy,1:qx);
    V = util.coeffs2vals( util.coeffs2vals(U).' ).';
    I = I + sum(sum(V .* ww0 * sclx * scly));
end

end
