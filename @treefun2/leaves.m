function ids = leaves(f)
%LEAVES   Get the leaf IDs in a TREEFUN2.
%   LEAVES(F) returns the leaf IDs in the TREEFUN2 F.

ids = [];
for k = 1:length(f.id)
    if ( f.height(k) == 0 )
        ids(end+1) = f.id(k); %#ok<AGROW>
    end
end

end
