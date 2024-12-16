function f = resample(f, n)

m = min(f.n, n);
ids = leaves(f);

for id = ids(:).'
    cfs = zeros(n);
    cfs(1:m,1:m) = f.coeffs{id}(1:m,1:m);
    f.coeffs{id} = cfs;
end

f.n = n;

end
