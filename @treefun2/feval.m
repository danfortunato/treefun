function vals = feval(f, x, y)
%FEVAL   Evaluate a TREEFUN2.

vals = zeros(size(x));
[ids, xs, ys, idxs] = pt2leaf(f, x, y);

for k = 1:length(ids)
    id = ids(k);
    dom = f.domain(:,id);
    sclx = 2/diff(dom(1:2));
    scly = 2/diff(dom(3:4));
    xm = sclx*(xs{k}-dom(1))-1;
    ym = scly*(ys{k}-dom(3))-1;
    if ( ~isempty(f.coeffs{id}) )
        vals(idxs{k}) = clenshaw2d(f.coeffs{id}, xm, ym);
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
