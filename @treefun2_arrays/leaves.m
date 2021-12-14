function leaves = leaves(f)
%LEAVES   Get the leaf boxes in a TREEFUN2.
%   LEAVES(F) returns the leaf boxes in the TREEFUN2 F.

leaves = [];
for k = 1:length(f.id)
    if ( f.height(k) == 0 )
        leaves(end+1) = f.id(k); %#ok<AGROW>
    end
end

end
