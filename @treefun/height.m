function out = height(f)
%HEIGHT   Length of a TREEFUN.
%   HEIGHT(F) returns the number of levels in the TREEFUN.

if ( isempty(f) )
    out = 0;
else
    out = max([f.boxes.height]) + 1;
end

end
