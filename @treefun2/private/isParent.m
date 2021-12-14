function out = isParent(a, b)
%ISPARENT   Is box B a parent of box A?

out = all(a.domain([1 3]) >= b.domain([1 3])) && ...
      all(a.domain([2 4]) <= b.domain([2 4])) && ...
     ~all(a.domain == b.domain);

end
