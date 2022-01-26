function f = uminus(f)
%-   Unary minus for TREEFUN2.
%   -F negates the TREEFUN2 F.
%
%   G = uminus(F) is called for the syntax '-F'.
%
%   See also UPLUS.

boxes = leaves(f);
for k = 1:length(boxes)
    id = boxes(k).id;
    f.boxes(id).coeffs = -f.boxes(id).coeffs;
end

end
