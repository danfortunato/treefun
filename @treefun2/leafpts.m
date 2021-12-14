function [xx, yy, ww] = leafpts(f)

[xx0, yy0] = chebpts2(f.n, f.n, [0 1 0 1]);
[~, w0] = chebpts(f.n, [0 1]);
ww0 = w0(:) * w0(:).';

leaf = leaves(f);
xx = cell(length(leaf), 1);
yy = cell(length(leaf), 1);
ww = cell(length(leaf), 1);
for k = 1:length(leaf)
    dom = leaf(k).domain;
    sclx = diff(dom(1:2));
    scly = diff(dom(3:4));
    xx{k} = sclx*xx0 + dom(1);
    yy{k} = scly*yy0 + dom(3);
    ww{k} = sclx*scly*ww0;
end

end
