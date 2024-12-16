function I = sum2(f, varargin)
%SUM2   Double integral of a TREEFUN2.
%   I = SUM2(F) returns the double integral of the TREEFUN2 F over its domain.
%
%   I = SUM2(F, 'all') returns an array of double integrals over each leaf of F.
%
% See also INTEGRAL, INTEGRAL2.

I = integral2(f, varargin);

end
