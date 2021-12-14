function out = isempty(f)
%ISEMPTY   Test for empty TREEFUN.
%   ISEMPTY(F) returns logical true if F is an empty TREEFUN and false
%   otherwise.

out = isempty(f.boxes);

end
