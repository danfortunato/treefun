function f = log10(f)
%LOG10   Base 10 logarithm of a TREEFUN2.
%   LOG10(F) returns the base 10 logarithm of the TREEFUN2 F. If F has any
%   roots in its domain, then the representation is likely to be
%   inaccurate.

% Empty check:
if ( isempty(f) )
    return
end

f = compose(@(x) real(log10(x)), f);

end
