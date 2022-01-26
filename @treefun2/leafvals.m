function vals = leafvals(f)

ids = leaves(f);
vals = cell(length(ids), 1);
for k = 1:length(ids)
    id = ids(k);
    vals{k} = treefun2.coeffs2vals(f.coeffs{id});
end

end
