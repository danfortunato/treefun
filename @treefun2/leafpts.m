function [xx, yy, ww] = leafpts(f, type)

if ( nargin < 2 )
    type = @(n) chebpts(n, [0 1], 2);
end

[x0, w0] = type(f.n);
[xx0, yy0] = meshgrid(x0);
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
