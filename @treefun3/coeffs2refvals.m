function vals = coeffs2refvals(coeffs)
%
%

persistent Eval pstored
p = size(coeffs, 1);
nrefpts = 2*p;

if ( isempty(Eval) || p ~= pstored )
    pstored = p;
    Eval = ones(nrefpts, p);
    x = linspace(-1, 1, nrefpts).';
    Eval(:,2) = x;
    for k=3:p
        Eval(:,k) = 2*x.*Eval(:,k-1)-Eval(:,k-2);
    end
end

tmp1 = permute(tensorprod(Eval,coeffs,2,1),[2 3 1]);
tmp2 = permute(tensorprod(Eval,tmp1,2,1),[2 3 1]);
vals = permute(tensorprod(Eval,tmp2,2,1),[2 3 1]);

end