function refvals = chebvals2refvals(chebvals)
%CHEBVALS2REFVALS   Convert values at tensor-product Chebyshev nodes to
% values at reference points.

persistent Eval pstored
p = size(chebvals, 1);
nrefpts = 2*p;
%nrefpts = p+2;

if ( isempty(Eval) || p ~= pstored )
    pstored = p;
    xcheb = chebpts(p, [-1 1]);
    xref = linspace(-1, 1, nrefpts).';
    Eval = barymat(xref, xcheb);
end

refvals = Eval * chebvals * Eval.';

end
