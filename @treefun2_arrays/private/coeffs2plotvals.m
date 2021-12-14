function vals = coeffs2plotvals(coeffs)
%COEFFS2PLOTVALS   Convert 2D Chebyshev coefficients to values
% at plot points.

persistent Eval
nplotpts = 100;
pmax = 20; % Maximum p to precompute

if ( isempty(Eval) )
    Eval = zeros(nplotpts, pmax);
    x = linspace(-1, 1, nplotpts).';
    c = zeros(pmax, 1);
    for k = 1:nplotpts
        c(k) = 1;
        Eval(:,k) = util.clenshaw(c, x);
        c(k) = 0;
    end
end

[px, py] = size(coeffs);
vals = Eval(:,1:py) * coeffs * Eval(:,1:px)';

end
