function f = mrdivide(f, c)
%/   Right matrix divide for TREEFUN2.
%   F/C divides the TREEFUN2 F by the scalar C.
%
% See also MLDIVIDE, RDIVIDE.

if ( isnumeric(c) && isscalar(c) )
    f = rdivide(f, c);
else
    error('TREEFUN2:mrdivide:invalid', 'Not supported. Did you mean ''./''?');
end

end
