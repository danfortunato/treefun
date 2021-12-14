function plotcoeffs(f)
%PLOTCOEFFS   Display coefficients graphically.
%   PLOTCOEFFS(F) plots the absolute values of the coefficients underlying
%   the representation of the TREEFUN2 F on a semilog scale.

holdState = ishold();

% Loop over the patches:
hold on
boxes = leaves(f);
for k = 1:length(boxes)
    box = boxes(k);
    %[x, y] = chebpts2(f.n, f.n, box.domain);
    [x, y] = meshgrid(linspace(box.domain(1), box.domain(2), 200), ...
                      linspace(box.domain(3), box.domain(4), 200));
    vals = abs(box.coeffs);
    coeffs = vals2coeffs(vals);
    u = coeffs2plotvals(coeffs);
    %stem3(x, y, abs(box.coeffs), 'ok', 'MarkerFaceColor', 'k');
    %stem3(x, y, u, 'o', 'MarkerFaceColor', 'k')
    surf(x, y, log10(abs(u)))
end
shading interp
view(3)
%set(gca, 'ZScale', 'log')
%set(gca,'ColorScale','log')

if ( ~holdState )
    hold off;
end

end
