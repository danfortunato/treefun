function f = diffy(f)
%DIFFY   Differentiate a TREEFUN2 with respect to y.
%   DIFFY(F) returns a TREEFUN2 representing the derivative of the TREEFUN2
%   F in its second argument.
%
%   DIFFY(F, N) returns a TREEFUN2 representing the N-th derivative of the
%   TREEFUN2 F in its second argument.
%
%   See also DIFFX, DIFF.

% Default to first derivative:
if ( nargin == 1 )
    n = 1;
end

f = diff(f, n, 1);

end
