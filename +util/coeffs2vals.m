function values = coeffs2vals(coeffs)
%COEFFS2VALS   Convert Chebyshev coefficients to values at Chebyshev points.
%
%   See also VALS2COEFFS.

% Store the coeffs2vals matrices for sizes < cutoff
persistent F
cutoff = 50;

% Get the length of the input:
n = size(coeffs, 1);

if ( n <= 1 )
    % Trivial case (constant):
    values = coeffs;
elseif ( n < cutoff )
    % Use matrix multiplication for small problems
    if ( isempty(F) )
        F = cell(cutoff, 1);
    end
    if ( isempty(F{n}) )
        F{n} = chebtech2.coeffs2vals(eye(n));
    end
    values = F{n} * coeffs;
else
    % Use fast transform
    values = chebtech2.coeffs2vals(coeffs);
end

end
