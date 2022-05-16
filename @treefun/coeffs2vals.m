function vals = coeffs2vals(coeffs)
%COEFFS2VALS   Convert Chebyshev coefficients to values at Chebyshev
% points.
%
%   See also VALS2COEFFS.

% Store the coeffs2vals matrix for sizes < cutoff
persistent F nstored
cutoff = 30;

% Get the length of the input:
n = size(coeffs, 1);

if ( n <= 1 )
    % Trivial case (constant):
    vals = coeffs;
elseif ( n < cutoff )
    % Use matrix multiplication for small problems
    if ( isempty(F) || n ~= nstored )
        nstored = n;
        F = chebtech2.coeffs2vals(eye(n));
    end
    vals = F * coeffs;
else
    % Use fast transform
    vals = chebtech2.coeffs2vals(coeffs);
end

end