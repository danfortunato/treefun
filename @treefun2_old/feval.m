function vals = feval(f, x, y)
%FEVAL   Evaluate a TREEFUN2.

sx = size(x);
x = x(:);
y = y(:);
vals = 0*x(:);
found = false(size(vals));

k = 1;
boxesToCheck = 1;
while ( ~all(found) )

    if ( k > length(boxesToCheck) )
        keyboard
    end
    
    id = boxesToCheck(k);
    dom = f.boxes(id).domain;
    idx = ((x > dom(1) & x < dom(2)) | (abs(x - dom(1)) < 1e-14) | (abs(x - dom(2)) < 1e-14)) ...
        & ((y > dom(3) & y < dom(4)) | (abs(y - dom(3)) < 1e-14) | (abs(y - dom(4)) < 1e-14));

    if ( any(idx) )
        if ( ~isLeaf(f.boxes(id)) )
            boxesToCheck = [boxesToCheck f.boxes(id).children]; %#ok<AGROW>
        else
            % Evaluate at points
            sclx = 2/diff(dom(1:2));
            scly = 2/diff(dom(3:4));
            xm = sclx*(x(idx)-dom(1))-1;
            ym = scly*(y(idx)-dom(3))-1;
            vals(idx) = clenshaw2d(f.boxes(id).coeffs, xm, ym);
            found(idx) = true;
        end
    end

    k = k + 1;
end

vals = reshape(vals, sx);

end

function v = clenshaw2d(C, x, y)
%CLENSHAW2D   Evaluate a 2D Chebyshev expansion at the given points.

v = 0*x;
Cy = chebtech2.clenshaw(y, C).';
for k = 1:numel(x)
    v(k) = chebtech2.clenshaw(x(k), Cy(:,k));
end

end
