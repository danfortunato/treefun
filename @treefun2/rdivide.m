function f = rdivide(f, g)
%./   Pointwise right divide for TREEFUN2.
%   F./G divides F by G, where F and G may be TREEFUN2 objects or scalars.
%
%   See also LDIVIDE, COMPOSE.

% If either F or G are empty then return an empty TREEFUN2 object.
if ( isempty(f) )
    return
elseif ( isempty(g) )
    f = g;
    return
end

if ( isa(f, 'treefun2') && isa(g, 'treefun2') )
    % Divide two TREEFUN2s:
    % TODO: Check that F and G have the same domain.
    [ss, wzero] = singleSignTest(g);
    if ( ss && ~wzero )
        f = compose(@rdivide, f, g);
    else
        error('TREEFUN2:rdivide:zero', ...
            'Attempting to invert a treefun2 with a root.');
    end
elseif ( isa(f, 'treefun2') && isnumeric(g) && isscalar(g) )
    % Divide TREEFUN2 F by scalar G:
    f = f .* (1/g);
elseif ( isnumeric(f) && isscalar(f) && isa(g, 'surfacefun') )
    % Divide scalar F by TREEFUN2 G:
    [ss, wzero] = singleSignTest(g);
    if ( ss && ~wzero )
        f = compose(@rdivide, f, g);
    else
        error('TREEFUN2:rdivide:zero', ...
            'Attempting to invert a treefun2 with a root.');
    end
else
    error('TREEFUN2:rdivide:invalid', ...
        'F and G must be scalars or treefun2 objects.')
end

end
