function coeffs = vals2coeffs(vals)
%VALS2COEFFS   Convert values at Chebyshev points to Chebyshev coefficients.
%
%   See also COEFFS2VALS.

% Store the vals2coeffs matrices for sizes < cutoff
persistent F
cutoff = 30;

% Get the length of the input:
p = size(vals, 1);

if ( p <= 1 )
    % Trivial case (constant):
    coeffs = vals;
elseif ( p < cutoff )
    % Use matrix multiplication for small problems
    if ( isempty(F) )
        F = cell(cutoff, 1);
    end
    if ( isempty(F{p}) )
        F{p} = chebtech2.vals2coeffs(eye(p));
    end
    coeffs = F{p} * vals * F{p}.';
else
    % Use fast transform
    coeffs = chebtech2.vals2coeffs( chebtech2.vals2coeffs(vals).' ).';
end

end
