function f = sqrt(f)
%SQRT   Square root of a TREEFUN2.
%   SQRT(F) returns the square root of the TREEFUN2 F.
%
%   See also POWER, COMPOSE.

% Empty check:
if ( isempty(f) )
    return
end

f = compose(@sqrt, f);

end
