function generateQuadrature(n)

doSave = false;

x = legpts(n, [0 1]);
[xx, yy] = meshgrid(x);
xx = xx(:);  yy = yy(:);
xf = 0.5*xx; yf = 0.5*yy;
xc = 2*xx;   yc = 2*yy;

        % Neighbor is the same level
nbrs = {[xx-1   yy];     % Left
        [xx-1   yy+1];   % Left up
        [xx     yy+1];   % Up
        [xx+1   yy+1];   % Right up
        [xx+1   yy];     % Right
        [xx+1   yy-1];   % Right down
        [xx     yy-1];   % Down
        [xx-1   yy-1];   % Left down
        % Neighbor is finer
        [xf-0.5 yf];     % Left with bottoms aligning
        [xf-0.5 yf+0.5]; % Left with tops aligning
        [xf-0.5 yf+1];   % Left up
        [xf     yf+1];   % Up with lefts aligning
        [xf+0.5 yf+1];   % Up with rights aligning
        [xf+1   yf+1];   % Right up
        [xf+1   yf+0.5]; % Right with tops aligning
        [xf+1   yf];     % Right with bottoms aligning
        [xf+1   yf-0.5]; % Right down
        [xf+0.5 yf-0.5]; % Down with rights aligning
        [xf     yf-0.5]; % Down with lefts aligning
        [xf-0.5 yf-0.5]; % Left down
        % Neighbor is coarser
        [xc-2   yc];     % Left with bottoms aligning
        [xc     yc+1];   % Up with lefts aligning
        [xc+1   yc-1];   % Right with tops aligning
        [xc-1   yc-2];   % Down with rights aligning
        [xc-2   yc-1];   % Left with tops aligning
        [xc-1   yc+1];   % Up with rights aligning
        [xc+1   yc];     % Right with bottoms aligning
        [xc     yc-2];   % Down with lefts aligning
        [xc-2   yc-2];   % Left down
        [xc-2   yc+1];   % Left up
        [xc+1   yc+1];   % Right up
        [xc+1   yc-2]};  % Right down

naiveMats = zeros(n^2, n^2, length(nbrs));
naiveMats(:,:,1) = laplace_kernel_self(xx, yy);
for k = 1:length(nbrs)
    naiveMats(:,:,k+1) = laplace_kernel(xx, yy, nbrs{k}(:,1), nbrs{k}(:,2));
end

% Load David's close quadrature matrices
closeMats = readNPY('~/Research/Stein/send_to_dan/box_code_lookups/laplace/laplace_16.npy');
closeMats = permute(closeMats, [4 5 2 3 1]);
closeMats = reshape(closeMats, [n^2 n^2 33]);
for k = 1:length(nbrs)+1
    for j = 1:n^2
        closeMats(:,j,k) = reshape(reshape(closeMats(:,j,k), n, n).', [], 1);
    end
    for j = 1:n^2
        closeMats(j,:,k) = reshape(reshape(closeMats(j,:,k), n, n).', [], 1);
    end
end

closeMats1 = zeros(n^2, n^2, 33);

% fprintf('# Precomputing self quadrature\n');
% kern = @(sx,sy,tx,ty) -log(sqrt((tx-sx).^2 + (ty-sy).^2)) / (2*pi);
% t = tic;
% for j = 0:n-1
%     for k = 0:n-1
%         f = @(x,y) mylegendre(j, 2*x-1) .* mylegendre(k, 2*y-1);
%         t1 = tic;
%         int = zeros(n^2, 1);
%         for l = 1:n^2
%             int(l) = integral2(@(sx,sy) kern(sx,sy,xx(l),yy(l)).*f(sx,sy), 0, 1, 0, 1, 'AbsTol', 1e-14, 'RelTol', 1e-14);
%         end
%         closeMats1(:,n*j+k+1,1) = int;
%         fprintf('  j=%d, k=%d: %gs\n', j, k, toc(t1));
%         % Save everything
%         if ( doSave )
%             save(sprintf('laplace_%d_new.mat', n), 'naiveMats', 'closeMats', 'closeMats1'); %#ok<UNRCH>
%         end
%     end
% end
% fprintf('# Total time for self: %gs\n', toc(t));

opts = struct('balance', false, 'neighbors', false);
kern = @(sx,sy,tx,ty) -log(sqrt((tx-sx).^2 + (ty-sy).^2)) / (2*pi);
for i = 1:length(nbrs)
    t = tic;
    fprintf('# Precomputing neighbor quadrature: %d / %d\n', i, length(nbrs));
    tx = nbrs{i}(:,1);
    ty = nbrs{i}(:,2);
    for j = 0:n-1
        for k = 0:n-1
            f = @(x,y) mylegendre(j, 2*x-1) .* mylegendre(k, 2*y-1);
            t1 = tic;
            %int = arrayfun(@(tx1,ty1) integral2(treefun2(@(sx,sy) kern(sx,sy,tx1,ty1).*f(sx,sy), [0 1 0 1], 32, opts)), tx, ty);
            int = zeros(n^2, 1);
            for l = 1:n^2
                int(l) = integral2(treefun2(@(sx,sy) kern(sx,sy,tx(l),ty(l)).*f(sx,sy), [0 1 0 1], 32, opts));
            end
            closeMats1(:,n*j+k+1,i+1) = int;
            fprintf('  j=%d, k=%d: %gs\n', j, k, toc(t1));
            % Save everything
            if ( doSave )
                save(sprintf('laplace_%d_new.mat', n), 'naiveMats', 'closeMats', 'closeMats1'); %#ok<UNRCH>
            end
        end
    end
    fprintf('# Total time for neighbor %d: %gs\n', toc(t));
end

% Save everything
if ( doSave )
    save(sprintf('laplace_%d.mat', n), 'naiveMats', 'closeMats'); %#ok<UNRCH>
end

end

function D = laplace_kernel_self(x, y)
% Laplace kernel with no self interaction
    D = -log(sqrt((x(:)-x(:).').^2 + (y(:)-y(:).').^2)) / (2*pi);
    D(1:size(D,1)+1:end) = 0;
end

function D = laplace_kernel(sx, sy, tx, ty)
    D = -log(sqrt((tx(:)-sx(:).').^2 + (ty(:)-sy(:).').^2)) / (2*pi);
end
