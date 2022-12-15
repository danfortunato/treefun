function [dx, dy] = fevald(f, x, y)
%FEVALD   Evaluate the derivative of a TREEFUN2.

persistent D DT
if ( isempty(D) || size(D,1) ~= f.n)
    S = ultraS.convertmat(f.n, 0, 0);
    Dultra = ultraS.diffmat(f.n, 1);
    D = full(S \ Dultra);
    DT = D.';
end

dx = zeros(size(x));
dy = zeros(size(x));
[ids, xs, ys, idxs] = pt2leaf(f, x, y);

for k = 1:length(ids)
    id = ids(k);
    dom = f.domain(:,id);
    sclx = 2/diff(dom(1:2));
    scly = 2/diff(dom(3:4));
    xm = sclx*(xs{k}-dom(1))-1;
    ym = scly*(ys{k}-dom(3))-1;
    if ( ~isempty(f.coeffs{id}) )
        dx_cfs = sclx * f.coeffs{id} * DT;
        dy_cfs = scly * D * f.coeffs{id};
        dx(idxs{k}) = clenshaw2d(dx_cfs, xm, ym);
        dy(idxs{k}) = clenshaw2d(dy_cfs, xm, ym);
    end
end

end

function v = clenshaw2d(C, x, y)
%CLENSHAW2D   Evaluate a 2D Chebyshev expansion at the given points.

v = 0*x;
Cy = chebtech2.clenshaw(y, C).';
for k = 1:numel(x)
    v(k) = chebtech2.clenshaw(x(k), Cy(:,k));
end

end
