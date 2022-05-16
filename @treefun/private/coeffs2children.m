function [L, R] = coeffs2children(coeffs)
%COEFFS2CHILDREN   Convert Chebyshev coefficients on a parent to 2D
% Chebyshev coefficients on its four children.

persistent EvalL EvalR nstored
n = size(coeffs, 1);

if ( isempty(EvalL) || isempty(EvalR) || n ~= nstored )
    nstored = n;
    x  = chebpts(n, [-1 1]);
    xL = chebpts(n, [-1 0]);
    xR = chebpts(n, [ 0 1]);
    baryL = barymat(xL, x);
    baryR = barymat(xR, x);
    C2V = chebtech2.coeffs2vals(eye(n));
    EvalL = chebtech2.vals2coeffs(baryL * C2V);
    EvalR = chebtech2.vals2coeffs(baryR * C2V);
end

L = EvalL * coeffs; % Left
R = EvalR * coeffs; % Right

end
