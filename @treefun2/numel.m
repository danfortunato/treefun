function N = numel(f)
%NUMEL   Number of degrees of freedom in a TREEFUN2.
%   N = NUMEL(F) returns the total number of degrees of freedom in the TREEFUN2
%   object F.

N = 0;
for k = 1:length(f)
    N = N + numel(f.coeffs{k});
end

end
