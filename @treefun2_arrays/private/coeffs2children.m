function [LL, LR, UL, UR] = coeffs2children(coeffs)
%COEFFS2CHILDREN   Convert 2D Chebyshev coefficients on a parent to 2D
% Chebyshev coefficients on its four children.

persistent EvalL EvalR pstored
p = size(coeffs, 1);

if ( isempty(EvalL) || isempty(EvalR) || p ~= pstored )
    pstored = p;
    x  = chebpts(p, [-1 1]);
    xL = chebpts(p, [-1 0]);
    xR = chebpts(p, [ 0 1]);
    baryL = barymat(xL, x);
    baryR = barymat(xR, x);
    C2V = chebtech2.coeffs2vals(eye(p));
    EvalL = chebtech2.vals2coeffs(baryL * C2V);
    EvalR = chebtech2.vals2coeffs(baryR * C2V);
end

LL = EvalL * coeffs * EvalL.'; % Lower left
LR = EvalL * coeffs * EvalR.'; % Lower right
UL = EvalR * coeffs * EvalL.'; % Upper left
UR = EvalR * coeffs * EvalR.'; % Upper right

end
