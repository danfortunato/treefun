function varargout = subsref(sol, index)
%SUBSREF   Subscripted reference for a TREEFUN2.
%   SOL(X, Y) returns the values of SOL evaluated at (X, Y). See FEVAL for
%   further details.
%
%   See also FEVAL.

idx = index(1).subs;
switch index(1).type
    case '()'
        % Where to evaluate:
        x = idx{1};
        y = idx{2};
        out = feval(sol, x, y);
        
    case '.'
        % Call GET() for .PROP access.
        out = sol.(idx);
        if ( numel(index) > 1 )
            % Recurse on SUBSREF():
            index(1) = [];
            out = subsref(out, index);
        end

    otherwise
        error('Unknown reference.');
end

% Convert to a cell:
varargout = {out};

end
