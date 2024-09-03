function [LLD, LRD, ULD, URD, LLT, LRT, ULT, URT] = coeffs2children3d(coeffs)
%COEFFS2CHILDREN   Convert 3D Chebyshev coefficients on a parent to 2D
% Chebyshev coefficients on its four children.
% follow Z curve, double check with ngrid

persistent EvalL EvalR Fp pstored
p = size(coeffs, 1);

if ( isempty(EvalL) || isempty(EvalR) || isempty(Fp) || p ~= pstored )
    pstored = p;
    x = -cos(pi*(2*(1:p)'-1)/(2*p));
    xL = 1/2*x-1/2;
    xR = 1/2*x+1/2;
    EvalL = ones(p, p); EvalR = ones(p, p); 
    EvalL(:,2) = xL; EvalR(:,2) = xR;
    for k=3:p
      EvalL(:,k) = 2*xL.*EvalL(:,k-1)-EvalL(:,k-2);
      EvalR(:,k) = 2*xR.*EvalR(:,k-1)-EvalR(:,k-2);
    end
    Fp = 2*cos(pi*((1:p)-1)'*(2*(p:-1:1)-1)/(2*p))/p;
    Fp(1,:) = 1/2*Fp(1,:);
end

% Lower left down
tmp1 = permute(tensorprod(EvalL,coeffs,2,1),[2 3 1 4]);
tmp2 = permute(tensorprod(EvalL,tmp1,2,1),[2 3 1 4]);
valsLLD = squeeze(permute(tensorprod(EvalL,tmp2,2,1),[2 3 1 4]));
tmp1hat = permute(tensorprod(Fp,valsLLD,2,1),[2 3 1 4]);
tmp2hat = permute(tensorprod(Fp,tmp1hat,2,1),[2 3 1 4]);
LLD     = permute(tensorprod(Fp,tmp2hat,2,1),[2 3 1 4]);
% Lower right down
% [xx0LRD, yy0LRD, zz0LRD] = ndgrid(xR,xL,xL);
tmp1 = permute(tensorprod(EvalR,coeffs,2,1),[2 3 1 4]);
tmp2 = permute(tensorprod(EvalL,tmp1,2,1),[2 3 1 4]);
valsLRD = squeeze(permute(tensorprod(EvalL,tmp2,2,1),[2 3 1 4]));
tmp1hat = permute(tensorprod(Fp,valsLRD,2,1),[2 3 1 4]);
tmp2hat = permute(tensorprod(Fp,tmp1hat,2,1),[2 3 1 4]);
LRD     = permute(tensorprod(Fp,tmp2hat,2,1),[2 3 1 4]);
% Upper left down
% [xx0ULD, yy0ULD, zz0ULD] = ndgrid(xL,xR,xL);
tmp1 = permute(tensorprod(EvalL,coeffs,2,1),[2 3 1 4]);
tmp2 = permute(tensorprod(EvalR,tmp1,2,1),[2 3 1 4]);
valsULD = squeeze(permute(tensorprod(EvalL,tmp2,2,1),[2 3 1 4]));
tmp1hat = permute(tensorprod(Fp,valsULD,2,1),[2 3 1 4]);
tmp2hat = permute(tensorprod(Fp,tmp1hat,2,1),[2 3 1 4]);
ULD     = permute(tensorprod(Fp,tmp2hat,2,1),[2 3 1 4]);
% Upper right down
% [xx0URD, yy0URD, zz0URD] = ndgrid(xR,xR,xL);
tmp1 = permute(tensorprod(EvalR,coeffs,2,1),[2 3 1 4]);
tmp2 = permute(tensorprod(EvalR,tmp1,2,1),[2 3 1 4]);
valsURD = squeeze(permute(tensorprod(EvalL,tmp2,2,1),[2 3 1 4]));
tmp1hat = permute(tensorprod(Fp,valsURD,2,1),[2 3 1 4]);
tmp2hat = permute(tensorprod(Fp,tmp1hat,2,1),[2 3 1 4]);
URD     = permute(tensorprod(Fp,tmp2hat,2,1),[2 3 1 4]);

% Lower left top
% [xx0LLT, yy0LLT, zz0LLT] = ndgrid(xL,xL,xR);
tmp1 = permute(tensorprod(EvalL,coeffs,2,1),[2 3 1 4]);
tmp2 = permute(tensorprod(EvalL,tmp1,2,1),[2 3 1 4]);
valsLLT = squeeze(permute(tensorprod(EvalR,tmp2,2,1),[2 3 1 4]));
tmp1hat = permute(tensorprod(Fp,valsLLT,2,1),[2 3 1 4]);
tmp2hat = permute(tensorprod(Fp,tmp1hat,2,1),[2 3 1 4]);
LLT     = permute(tensorprod(Fp,tmp2hat,2,1),[2 3 1 4]);
% Lower right top
% [xx0LRT, yy0LRT, zz0LRT] = ndgrid(xR,xL,xR);
tmp1 = permute(tensorprod(EvalR,coeffs,2,1),[2 3 1 4]);
tmp2 = permute(tensorprod(EvalL,tmp1,2,1),[2 3 1 4]);
valsLRT = squeeze(permute(tensorprod(EvalR,tmp2,2,1),[2 3 1 4]));
tmp1hat = permute(tensorprod(Fp,valsLRT,2,1),[2 3 1 4]);
tmp2hat = permute(tensorprod(Fp,tmp1hat,2,1),[2 3 1 4]);
LRT     = permute(tensorprod(Fp,tmp2hat,2,1),[2 3 1 4]);
% Upper left top
% [xx0ULT, yy0ULT, zz0ULT] = ndgrid(xL,xR,xR);
tmp1 = permute(tensorprod(EvalL,coeffs,2,1),[2 3 1 4]);
tmp2 = permute(tensorprod(EvalR,tmp1,2,1),[2 3 1 4]);
valsULT = squeeze(permute(tensorprod(EvalR,tmp2,2,1),[2 3 1 4]));
tmp1hat = permute(tensorprod(Fp,valsULT,2,1),[2 3 1 4]);
tmp2hat = permute(tensorprod(Fp,tmp1hat,2,1),[2 3 1 4]);
ULT     = permute(tensorprod(Fp,tmp2hat,2,1),[2 3 1 4]);
% Upper right top
% [xx0URT, yy0URT, zz0URT] = ndgrid(xR,xR,xR);
tmp1 = permute(tensorprod(EvalR,coeffs,2,1),[2 3 1 4]);
tmp2 = permute(tensorprod(EvalR,tmp1,2,1),[2 3 1 4]);
valsURT = squeeze(permute(tensorprod(EvalR,tmp2,2,1),[2 3 1 4]));
tmp1hat = permute(tensorprod(Fp,valsURT,2,1),[2 3 1 4]);
tmp2hat = permute(tensorprod(Fp,tmp1hat,2,1),[2 3 1 4]);
URT     = permute(tensorprod(Fp,tmp2hat,2,1),[2 3 1 4]);

end
