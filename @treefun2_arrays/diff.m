function f = diff(f, n, dim)
%DIFF   Differentiate a TREEFUN2
%   DIFF(F, N, DIM) is the N-th derivative of the TREEFUN2 F in the
%   dimension DIM.
%      DIM = 1 (default) is the derivative in the y direction.
%      DIM = 2 is the derivative in the x direction.
%
%   DIFF(F, [NX NY]) takes the NX-th derivative of F in the x direction
%   and the NY-th derivative of F in the y direction.
%
%   See also DIFFX, DIFFY.

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

% Default to partial derivative in y:
if ( nargin < 3 )
    dim = 1;
elseif ( numel(dim) ~= 1 )
    error('TREEFUN2:diff:dim', 'Dimension should be either 1 or 2.');
end

% Make sure N is in [NX NY] format
if ( isscalar(n) )
    if ( dim == 1 )
        n = [0 n];
    elseif ( dim == 2 )
        n = [n 0];
    else
        error('TREEFUN2:diff:dim', 'Dimension should be either 1 or 2.');
    end
end
dx = n(1);
dy = n(2);
m = f.n;

boxes = leaves(f);
for k = 1:length(boxes)
    coeffs = boxes(k).coeffs;
    dom = boxes(k).domain;
    for l = 1:dy, coeffs = [ cdiff(coeffs);     zeros(1, m) ]; end
    for l = 1:dx, coeffs = [ cdiff(coeffs.').', zeros(m, 1) ]; end
    % Rescale from [-1,1]^2:
    sclx = (2/diff(dom(1:2)))^dx;
    scly = (2/diff(dom(3:4)))^dy;
    % Assign to the true leaf box via its box ID:
    id = boxes(k).id;
    f.boxes(id).coeffs = sclx*scly*coeffs;
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
