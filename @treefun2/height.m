function out = height(f)
%HEIGHT   Length of a TREEFUN2.
%   HEIGHT(F) returns the number of levels in the TREEFUN2.

if ( isempty(f) )
    out = 0;
else
    out = max([f.boxes.height]) + 1;
end

end
