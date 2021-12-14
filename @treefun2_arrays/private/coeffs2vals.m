function vals = coeffs2vals(coeffs)
%COEFFS2VALS   Convert Chebyshev coefficients to values at Chebyshev points.
%
%   See also VALS2COEFFS.

% Store the coeffs2vals matrices for sizes < cutoff
persistent F
cutoff = 30;

% Get the length of the input:
p = size(coeffs, 1);

if ( p <= 1 )
    % Trivial case (constant):
    vals = coeffs;
elseif ( p < cutoff )
    % Use matrix multiplication for small problems
    if ( isempty(F) )
        F = cell(cutoff, 1);
    end
    if ( isempty(F{p}) )
        F{p} = chebtech2.coeffs2vals(eye(p));
    end
    vals = F{p} * coeffs * F{p}.';
else
    % Use fast transform
    vals = chebtech2.coeffs2vals( chebtech2.coeffs2vals(coeffs).' ).';
end

end
