function out = isParent(f, a, b)
%ISPARENT   Is box B a parent of box A?

out = all(f.domain(a, [1 3]) >= f.domain(b, [1 3])) && ...
      all(f.domain(a, [2 4]) <= f.domain(b, [2 4])) && ...
     ~all(f.domain(a,:) == f.domain(b,:));

end
