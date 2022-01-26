function f = balance_v1(f, func)

L = leaves(f);
k = 1;
while ( k <= length(L) )
    id = L(k);
    neighborIDs = neighbors(f, id, 'leaf');
    neighborIDs = [neighborIDs{:}];
    if ( any(f.level(neighborIDs) > f.level(id)+1) )
        f = refine(f, id);
        f.height(id) = 1;
        L = [L f.children(:,id).']; %#ok<AGROW>
        idx = f.level(neighborIDs) < f.level(id);
        L = [L neighborIDs(idx)]; %#ok<AGROW>
    end
    k = k + 1;
end

% Do a cumulative sum in reverse to correct the heights
for k = length(f.id):-1:1
    if ( ~isLeaf(f, k) )
        f.height(k) = 1 + max(f.height(f.children(:,k)));
    end
end

% Now fill in the leaf coefficients
L = leaves(f);
[xx0, yy0] = chebpts2(f.n, f.n, [0 1 0 1]);
for k = 1:length(L)
    id = L(k);
    if ( isempty(f.coeffs{id}) )
        dom = f.domain(:,id);
        sclx = diff(dom(1:2));
        scly = diff(dom(3:4));
        xx = sclx*xx0 + dom(1);
        yy = scly*yy0 + dom(3);
        vals = func(xx,yy);
        f.coeffs{id} = treefun2.vals2coeffs(vals);
    end
end

end
