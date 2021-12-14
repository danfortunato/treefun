function vals = leafvals(f)

leaf = leaves(f);
vals = cell(length(leaf), 1);
for k = 1:length(leaf)
    vals{k} = coeffs2vals(leaf(k).coeffs);
end

end
