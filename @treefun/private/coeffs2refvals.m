function vals = coeffs2refvals(coeffs)
%COEFFS2REFVALS   Convert Chebyshev coefficients to values at reference
% points.

persistent Eval nstored
n = size(coeffs, 1);
nrefpts = 2*n;

if ( isempty(Eval) || n ~= nstored )
    nstored = n;
    x = linspace(-1, 1, nrefpts).';
    Eval = chebtech2.clenshaw(x, eye(n));
end

vals = Eval * coeffs;

end
