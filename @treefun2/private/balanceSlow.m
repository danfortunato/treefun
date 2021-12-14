function f = balanceSlow(f, func)

L = [leaves(f).id];
k = 1;
while ( k <= length(L) )
    id = L(k);
    neighborIDs = neighbors(f, id, 'leaf');
    neighborIDs = [neighborIDs{:}];
    if ( any([f.boxes(neighborIDs).level] > f.boxes(id).level+1) )
        f = refine(f, id);
        f.boxes(id).height = 1;
        L = [L f.boxes(id).children]; %#ok<AGROW>
        idx = [f.boxes(neighborIDs).level] < f.boxes(id).level;
        L = [L neighborIDs(idx)]; %#ok<AGROW>
    end
    k = k + 1;
end

% Do a cumulative sum in reverse to correct the heights
for k = length(f.boxes):-1:1
   if ( ~isLeaf(f.boxes(k).height) )
       f.boxes(k).height = 1 + max([f.boxes(f.boxes(k).children).height]);
   end
end

% Now fill in the leaf coefficients
L = [leaves(f).id];
[xx0, yy0] = chebpts2(f.n, f.n, [0 1 0 1]);
for k = 1:length(L)
    id = L(k);
    if ( isempty(f.boxes(id).coeffs) )
        dom = f.boxes(id).domain;
        sclx = diff(dom(1:2));
        scly = diff(dom(3:4));
        xx = sclx*xx0 + dom(1);
        yy = scly*yy0 + dom(3);
        vals = func(xx,yy);
        f.boxes(id).coeffs = vals2coeffs(vals);
    end
end

end
