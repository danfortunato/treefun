function coeffs = vals2coeffs(vals)
%VALS2COEFFS   Convert values at Chebyshev points to Chebyshev
% coefficients.
%
%   See also COEFFS2VALS.

% Store the vals2coeffs matrix for sizes < cutoff
persistent F nstored
cutoff = 30;

% Get the length of the input:
n = size(vals, 1);

if ( n <= 1 )
    % Trivial case (constant):
    coeffs = vals;
elseif ( n < cutoff )
    % Use matrix multiplication for small problems
    if ( isempty(F) || n ~= nstored )
        nstored = n;
        F = chebtech2.vals2coeffs(eye(n));
    end
    coeffs = F * vals;
else
    % Use fast transform
    coeffs = chebtech2.vals2coeffs(vals);
end

end