function out = isParent(a, b)
%ISPARENT   Is box B a parent of box A?

out = all(a.domain(1) >= b.domain(1)) && ...
      all(a.domain(2) <= b.domain(2)) && ...
     ~all(a.domain == b.domain);

end
