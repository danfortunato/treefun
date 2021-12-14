function v = clenshaw2d(C, x, y)
%CLENSHAW2D   Evaluate a 2D Chebyshev expansion at the given points.

v = 0*x;
Cy = chebtech2.clenshaw(y, C).';
for k = 1:numel(x)
    v(k) = chebtech2.clenshaw(x(k), Cy(:,k));
end

end
