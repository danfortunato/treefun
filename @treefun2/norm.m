function normF = norm(f, varargin)
%NORM   Norm of a TREEFUN2
%   For TREEFUN2 objects:
%       NORM(F) = sqrt(integral of abs(F)^2).
%       NORM(F, 2) is the same as NORM(F).
%       NORM(F, 'fro') is also the same as NORM(F).
%       NORM(F, 1) = integral of abs(F).
%       NORM(F, P) = (integral of abs(F)^P)^(1/P).
%       NORM(F, inf) = estimated global maximum in absolute value.
%       NORM(F, -inf) = estimated global minimum in absolute value.
%       NORM(F, 'max') is the same as NORM(F, inf).
%       NORM(F, 'min') is the same as NORM(F, -inf).
%       NORM(F, 'H1') = sqrt( ||u||_2^2 + ||grad(u)||_2^2 )
%       NORM(F, 'H2') = sqrt( ||u||_H1^2 + ||grad(diffx(u))||_2^2
%                                          + ||grad(diffy(u))||_2^2 )
%       NORM(F, 'lap') = sqrt( ||u||_2 + ||lap(u)||_2^2 )
%
%   NORM(F, 'all') and NORM(F, P, 'all') return an array of norms on each
%   box of F.

% Parse arguments.
p = 2;
reduce = true;
if ( nargin == 2 )
    if ( strcmp(varargin{1}, 'all') )
        reduce = false;
    else
        p = varargin{1};
    end
elseif ( nargin == 3 )
    p = varargin{1};
    if ( strcmp(varargin{2}, 'all') )
        reduce = false;
    end
end

% Empty TREEFUN2 has norm 0.
if ( isempty(f) )
    normF = 0;
    return
end

switch ( p )
    case 'fro'
        normF = norm(f, 2);
        reduceFun = @(x) sqrt( sum(x.^2) );

    case {inf, 'inf', 'max'}
        ids = leaves(f);
        normF = zeros(length(ids), 1);
        for k = 1:length(ids)
            id = ids(k);
            V = treefun2.coeffs2vals(f.coeffs{id});
            normF(k) = max(abs(V(:)));
        end
        reduceFun = @max;

    case {-inf, '-inf', 'min'}
        ids = leaves(f);
        normF = zeros(length(ids), 1);
        for k = 1:length(ids)
            id = ids(k);
            V = treefun2.coeffs2vals(f.coeffs{id});
            normF(k) = min(abs(V(:)));
        end
        reduceFun = @min;

    case 'H1'
        [fx, fy] = grad(f);
        norm_L2 = norm(f,   'all');
        norm_x  = norm(fx, 'all');
        norm_y  = norm(fy, 'all');
        normF = sqrt( norm_L2.^2 + norm_x.^2 + norm_y.^2 );
        reduceFun = @(x) sqrt( sum(x.^2) );

    case 'H2'
        [fx, fy]   = grad(f);
        [fxx, fxy] = grad(fx);
        [fyx, fyy] = grad(fy);
        norm_H1 = norm(f, 'H1', 'all');
        norm_xx = norm(fxx, 'all'); norm_xy  = norm(fxy, 'all');
        norm_yx = norm(fyx, 'all'); norm_yy  = norm(fyy, 'all');
        normF = sqrt( norm_H1.^2 + norm_xx.^2 + norm_xy.^2 + ...
                                   norm_yx.^2 + norm_yy.^2 );
        reduceFun = @(x) sqrt( sum(x.^2) );

    case 'lap'
        norm_L2  = norm(f,      'all');
        norm_lap = norm(lap(f), 'all');
        normF = sqrt( norm_L2.^2 + norm_lap.^2 );
        reduceFun = @(x) sqrt( sum(x.^2) );

    otherwise
        if ( isnumeric(p) && isreal(p) )
            if ( abs(round(p) - p) < eps )
                normF = integrate(f, p);
                reduceFun = @(x) ( sum(x.^p) ).^(1/p);
            else
                error('TREEFUN2:norm:norm', ...
                    'TREEFUN2 does not support this norm.');
            end
        else
            error('TREEFUN2:norm:unknown', 'Unknown norm.');
        end

end

% Combine norms on each element.
if ( reduce )
    normF = reduceFun(normF);
end

end

function int = integrate(f, p)
%INTEGRATE   Compute (integral of abs(F)^P)^(1/P)

% P should be an integer.
p = round(p);

% If a box uses an N x N discretization, then quadrature is performed on
% that box using N*P points.
nx = f.n;
ny = f.n;
qx = nx*p;
qy = ny*p;
wx = chebtech2.quadwts(qx); wx = wx(:);
wy = chebtech2.quadwts(qy); wy = wy(:);

ids = leaves(f);
int = zeros(length(ids), 1);
for k = 1:length(ids)
    id = ids(k);
    U = zeros(qy, qx);
    U(1:ny,1:nx) = f.coeffs{id};
    V = treefun2.coeffs2vals(U);
    int(k) = sum(sum(abs(V).^p .* wy .* wx.')).^(1/p);
end

end
