function f = uminus(f)
%-   Unary minus for TREEFUN2.
%   -F negates the TREEFUN2 F.
%
%   G = uminus(F) is called for the syntax '-F'.
%
%   See also UPLUS.

ids = leaves(f);
coeffs = f.coeffs;
for id = ids(:).'
    coeffs{id} = -coeffs{id};
end
f.coeffs = coeffs;

end
