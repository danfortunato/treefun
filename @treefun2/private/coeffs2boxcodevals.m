function vals = coeffs2boxcodevals(coeffs)
%COEFFS2BOXCODEVALS   Convert 2D Chebyshev coefficients to values at box
% code points.

% These are the points that the box code will use to sample the function on
% each leaf box:
%
%    x0 = linspace(0+1/8, 1-1/8, 4);
%    [xx0, yy0] = meshgrid(x0);
%
% The order of values starts from the lower left corner of the box and
% proceeds upwards row by row from left to right.

persistent Eval
nboxpts = 4;
pmax = 20; % Maximum p to precompute

if ( isempty(Eval) )
    x = linspace(-1+1/4, 1-1/4, nboxpts).';
    [xc, ~, vc] = chebpts(nboxpts);
    Eval = barymat(x, xc, vc);
%     Eval = zeros(nboxpts, pmax);
%     x = linspace(0+1/8, 1-1/8, nboxpts).';
%     c = zeros(pmax, 1);
%     for k = 1:pmax
%         c(k) = 1;
%         Eval(:,k) = chebtech2.clenshaw(x, c);
%         c(k) = 0;
%     end
end

%[px, py] = size(coeffs);
%vals = Eval(:,1:py) * coeffs * Eval(:,1:px).';
chebvals = coeffs2vals(coeffs);
vals = Eval * chebvals * Eval.';

end
