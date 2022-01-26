function vals = coeffs2refvals(coeffs)
%COEFFS2REFVALS   Convert 2D Chebyshev coefficients to values at reference
% points.

persistent Eval pstored
p = size(coeffs, 1);
nrefpts = 2*p;

if ( isempty(Eval) || p ~= pstored )
    pstored = p;
    Eval = zeros(nrefpts, p);
    x = linspace(-1, 1, nrefpts).';
    c = zeros(p, 1);
    for k = 1:p
        c(k) = 1;
        Eval(:,k) = chebtech2.clenshaw(x, c);
        c(k) = 0;
    end
end

vals = Eval * coeffs * Eval.';

end
