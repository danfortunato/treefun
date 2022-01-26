function f = abs(f)
%ABS   Absolute value of a TREEFUN2.
%   ABS(F) returns the absolute value of the TREEFUN2 F.

% Empty check:
if ( isempty(f) )
    return
end

f = compose(@abs, f);

end
