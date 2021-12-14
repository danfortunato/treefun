function f = diff(f, n)
%DIFF   Differentiate a TREEFUN.
%   DIFF(F, N) is the N-th derivative of the TREEFUN F.

% Empty check:
if ( isempty(f) )
    return
end

% Default to first derivative:
if ( nargin < 2 || isempty(n) )
    n = 1;
elseif ( n == 0 )
    return
end

m = f.n;
boxes = leaves(f);
for k = 1:length(boxes)
    coeffs = boxes(k).coeffs;
    dom = boxes(k).domain;
    for l = 1:n
        coeffs = [ cdiff(coeffs); zeros(m, 1) ];
    end
    % Rescale from [-1,1]:
    sclx = (2/diff(dom))^n;
    % Assign to the true leaf box via its box ID:
    id = boxes(k).id;
    f.boxes(id).coeffs = sclx*coeffs;
end

end

function dC = cdiff(C)
%CDIFF   Recurrence relation for coefficients of derivative.
%   CDIFF(C) returns the matrix of Chebyshev coefficients whose columns are
%   the derivatives of the columns of C.

[n, m] = size(C);
dC = zeros(n-1, m);                        % Initialize vector {c_r}
w = repmat(2*(1:n-1)', 1, m);
v = w.*C(2:end,:);                         % Temporal vector
dC(n-1:-2:1,:) = cumsum(v(n-1:-2:1,:), 1); % Compute c_{n-2}, c_{n-4}, ...
dC(n-2:-2:1,:) = cumsum(v(n-2:-2:1,:), 1); % Compute c_{n-3}, c_{n-5}, ...
dC(1,:) = .5*dC(1,:);                      % Adjust the value for c_0

end
