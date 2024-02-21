function coeffs = vals2coeffs(vals)
%VALS2COEFFS   Convert values at Chebyshev points to Chebyshev
% coefficients.
%
%   See also COEFFS2VALS.

% Store the vals2coeffs matrices for sizes < cutoff
persistent F
cutoff = 1e+03; % ... later

% Get the length of the input:
[p,~,~,nd] = size(vals); % vals of dim p^3 x nd

if ( p <= 1 )
    % Trivial case (constant):
    coeffs = vals;
elseif ( p < cutoff )
    % Use matrix multiplication for small problems
    if ( isempty(F) )
        F = cell(cutoff, 1);
    end
    if ( isempty(F{p}) )
        F{p} = 2*cos(pi*((1:p)-1)'*(2*(p:-1:1)-1)/(2*p))/p;
        F{p}(1,:) = 1/2*F{p}(1,:);
    end
    tmp1hat = permute(tensorprod(F{p},vals,2,1),[2 3 1 4]);
    tmp2hat = permute(tensorprod(F{p},tmp1hat,2,1),[2 3 1 4]);
    coeffs  = permute(tensorprod(F{p},tmp2hat,2,1),[2 3 1 4]);
else
    % Use fast transform ... not yet
    % coeffs = chebtech2.vals2coeffs( chebtech2.vals2coeffs(vals).' ).';
end

end
