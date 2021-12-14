function val = mylegendre(n, x)
    npts = numel(x);
    val = zeros(npts, 1);
    mex_id_ = 'computeLegendre(i int, i int, i double[], io double[])';
[val] = gateway(mex_id_, n, npts, x, val);
    val = reshape(val, size(x));
end
