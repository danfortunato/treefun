function L = laplacian(f)
%LAPLACIAN   Laplacian of a TREEFUN2.
%   LAPLACIAN(F) returns a TREEFUN2 representing the Laplacian of the
%   TREEFUN2 F.
%
%   See also LAP.

L = diff(f, 2, 2) + diff(f, 2, 1);

end
