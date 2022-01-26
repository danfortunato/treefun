function [xx, yy, ww] = leafpts(f)

[xx0, yy0] = chebpts2(f.n, f.n, [0 1 0 1]);
[~, w0] = chebpts(f.n, [0 1]);
ww0 = w0(:) * w0(:).';

ids = leaves(f);
xx = cell(length(ids), 1);
yy = cell(length(ids), 1);
ww = cell(length(ids), 1);
for k = 1:length(ids)
    id = ids(k);
    dom = f.domain(:,id);
    sclx = diff(dom(1:2));
    scly = diff(dom(3:4));
    xx{k} = sclx*xx0 + dom(1);
    yy{k} = scly*yy0 + dom(3);
    ww{k} = sclx*scly*ww0;
end

end
