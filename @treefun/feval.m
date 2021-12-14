function vals = feval(f, x)
%FEVAL   Evaluate a TREEFUN.

sx = size(x);
x = x(:);
vals = 0*x(:);
found = false(size(vals));

k = 1;
boxesToCheck = 1;
while ( ~all(found) )
    
    id = boxesToCheck(k);
    dom = f.boxes(id).domain;
    idx = (x > dom(1) & x < dom(2)) | (abs(x - dom(1)) < 1e-14) | (abs(x - dom(2)) < 1e-14);

    if ( any(idx) )
        if ( ~isLeaf(f.boxes(id)) )
            boxesToCheck = [boxesToCheck f.boxes(id).children]; %#ok<AGROW>
        else
            % Evaluate at points
            sclx = 2/diff(dom);
            xm = sclx*(x(idx)-dom(1))-1;
            vals(idx) = chebtech2.clenshaw(xm, f.boxes(id).coeffs);
            found(idx) = true;
        end
    end

    k = k + 1;
end

vals = reshape(vals, sx);

end
