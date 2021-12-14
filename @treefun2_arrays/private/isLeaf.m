function out = isLeaf(f, id)
%ISLEAF   Is this box a leaf?

out = ( f.height(id) == 0 );

end
