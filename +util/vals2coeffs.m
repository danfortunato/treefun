function coeffs = vals2coeffs(values)
%VALS2COEFFS   Convert values at Chebyshev points to Chebyshev coefficients.
%
%   See also COEFFS2VALS.

% Store the vals2coeffs matrices for sizes < cutoff
persistent F
cutoff = 50;

% Get the length of the input:
n = size(values, 1);

if ( n <= 1 )
    % Trivial case (constant):
    coeffs = values;
elseif ( n < cutoff )
    % Use matrix multiplication for small problems
    if ( isempty(F) )
        F = cell(cutoff, 1);
    end
    if ( isempty(F{n}) )
        F{n} = chebtech2.vals2coeffs(eye(n));
    end
    coeffs = F{n} * values;
else
    % Use fast transform
    coeffs = chebtech2.vals2coeffs(values);
end

end
