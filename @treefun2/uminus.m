function f = uminus(f)
%-   Unary minus for TREEFUN2.
%   -F negates the TREEFUN2 F.
%
%   G = uminus(F) is called for the syntax '-F'.
%
%   See also UPLUS.

ids = leaves(f);
for id = ids(:).'
    f.coeffs{id} = -f.coeffs{id};
end

end
