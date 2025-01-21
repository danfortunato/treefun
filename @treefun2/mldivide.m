function f = mldivide(c, f)
%/   Left matrix divide for TREEFUN2V.
%   C\F divides the TREEFUN2 F by the scalar C.
%
% See also MRDIVIDE, LDIVIDE.

if ( isnumeric(c) && isscalar(c) )
    f = ldivide(c, f);
else
    error('TREEFUN2:mldivide:invalid', 'Not supported. Did you mean ''./''?');
end

end
