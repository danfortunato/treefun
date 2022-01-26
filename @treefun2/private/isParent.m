function out = isParent(f, a, b)
%ISPARENT   Is box B a parent of box A?

out = all(f.domain([1 3], a) >= f.domain([1 3], b)) && ...
      all(f.domain([2 4], a) <= f.domain([2 4], b)) && ...
     ~all(f.domain(:,a) == f.domain(:,b));

end
