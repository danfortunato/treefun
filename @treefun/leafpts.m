function [x, w] = leafpts(f)

[x0, w0] = chebpts(f.n, [0 1]);

leaf = leaves(f);
x = cell(length(leaf), 1);
w = cell(length(leaf), 1);
for k = 1:length(leaf)
    dom = leaf(k).domain;
    sclx = diff(dom);
    x{k} = sclx*x0 + dom(1);
    w{k} = sclx*w0;
end

end
