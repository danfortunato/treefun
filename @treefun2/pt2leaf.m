function [ids, xs, ys, idxs] = pt2leaf(f, x, y)
%PT2LEAF   Find the leaves containing the given points.

[ids, xs, ys, idxs] = ipt2leaf(f, 1, x(:), y(:), (1:numel(x)).');

end

function [ids, xs, ys, idxs] = ipt2leaf(f, root, x, y, idx)

in = ( x>=f.domain(1,root) & x<=f.domain(2,root) & ...
       y>=f.domain(3,root) & y<=f.domain(4,root) );

if ( any(in) )
    xin = x(in);
    yin = y(in);
    idxin = idx(in);
    if ( f.height(root) == 0 )
        ids = root;
        xs = {xin};
        ys = {yin};
        idxs = {idxin};
    else
        [id1, xs1, ys1, idxs1] = ipt2leaf(f, f.children(1,root), xin, yin, idxin);
        [id2, xs2, ys2, idxs2] = ipt2leaf(f, f.children(2,root), xin, yin, idxin);
        [id3, xs3, ys3, idxs3] = ipt2leaf(f, f.children(3,root), xin, yin, idxin);
        [id4, xs4, ys4, idxs4] = ipt2leaf(f, f.children(4,root), xin, yin, idxin);
        ids  = [id1 ; id2 ; id3 ; id4];
        xs   = [xs1 ; xs2 ; xs3 ; xs4];
        ys   = [ys1 ; ys2 ; ys3 ; ys4];
        idxs = [idxs1 ; idxs2 ; idxs3 ; idxs4];
    end
else
    ids  = [];
    xs   = {};
    ys   = {};
    idxs = {};
end

end
