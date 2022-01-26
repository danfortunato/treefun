function vals = coeffs2refvals(coeffs)
%COEFFS2REFVALS   Convert 2D Chebyshev coefficients to values at reference
% points.

persistent Eval pstored
p = size(coeffs, 1);
nrefpts = 2*p;

if ( isempty(Eval) || p ~= pstored )
    pstored = p;
    x = linspace(-1, 1, nrefpts).';
    Eval = chebtech2.clenshaw(x, eye(p));
end

vals = Eval * coeffs * Eval.';

end
