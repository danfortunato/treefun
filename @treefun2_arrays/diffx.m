function f = diffx(f, n)
%DIFFX   Differentiate a TREEFUN2 with respect to x.
%   DIFFX(F) returns a TREEFUN2 representing the derivative of the TREEFUN2
%   F in its first argument.
%
%   DIFFX(F, N) returns a TREEFUN2 representing the N-th derivative of the
%   TREEFUN2 F in its first argument.
%
%   See also DIFFY, DIFF.

% Default to first derivative:
if ( nargin == 1 )
    n = 1;
end

f = diff(f, n, 2);

end
