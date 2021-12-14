function leaves = leaves(f)
%LEAVES   Get the leaf boxes in a TREEFUN2.
%   LEAVES(F) returns the leaf boxes in the TREEFUN2 F.

ids = [];
for k = 1:length(f.boxes)
    if ( f.boxes(k).height == 0 )
        ids(end+1) = k; %#ok<AGROW>
    end
end
leaves = f.boxes(ids);

end
