function plotcoeffs(f)
%PLOTCOEFFS   Display coefficients graphically.
%   PLOTCOEFFS(F) plots the absolute values of the coefficients underlying
%   the representation of the TREEFUN2 F on a semilog scale.

holdState = ishold();
nplotpts = 100;

% Loop over the patches:
hold on
ids = leaves(f);
for k = 1:length(ids)
    id = ids(k);
    %[x, y] = chebpts2(f.n, f.n, f.domain(:,id));
    [x, y] = meshgrid(linspace(f.domain(1,id), f.domain(2,id), nplotpts), ...
                      linspace(f.domain(3,id), f.domain(4,id), nplotpts));
    vals = abs(f.coeffs{id});
    coeffs = treefun2.vals2coeffs(vals);
    u = coeffs2plotvals(coeffs);
    %stem3(x, y, abs(f.coeffs{id}), 'ok', 'MarkerFaceColor', 'k');
    stem3(x, y, u, 'o', 'MarkerFaceColor', 'k')
    %surf(x, y, log10(abs(u)))
end
shading interp
view(3)
%set(gca, 'ZScale', 'log')
%set(gca,'ColorScale','log')

if ( ~holdState )
    hold off;
end

end
