function I = integral(f, varargin)
%INTEGRAL   Double integral of a TREEFUN2.
%   I = INTEGRAL(F) returns the double integral of the TREEFUN2 F over its
%   domain.
%
%   I = INTEGRAL(F, 'all') returns an array of double integrals over each leaf
%   of F.
%
% See also INTEGRAL2, SUM2.

I = integral2(f, varargin);

end
