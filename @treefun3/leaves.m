function ids = leaves(f)
%LEAVES   Get the leaf IDs in a TREEFUN2.
%   LEAVES(F) returns the leaf IDs in the TREEFUN2 F.

ids = f.id(f.height == 0);

end
