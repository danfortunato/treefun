function vals = coeffs2plotvals(coeffs)
%COEFFS2PLOTVALS   Convert 2D Chebyshev coefficients to values at plot
% points.

persistent Eval nstored
nplotpts = 100;

n = size(coeffs, 1);
if ( isempty(Eval) || n ~= nstored )
    nstored = n;
    x = linspace(-1, 1, nplotpts).';
    Eval = chebtech2.clenshaw(x, eye(n));
end

vals = Eval * coeffs * Eval.';

end
