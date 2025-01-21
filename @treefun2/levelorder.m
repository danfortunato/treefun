function f = levelorder(f)
%LEVELORDER   Put a TREEFUN2 in level order.
%   LEVELORDER(F) returns a TREEFUN2 representing the same tree as F but whose
%   boxes are sorted into level order.

lvls = levels(f);
f = reindex(permute(f, [lvls{:}]));

end
