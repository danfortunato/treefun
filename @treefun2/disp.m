function disp(f)
%DISP   Display a TREEFUN2 to the command line.
% 
% See also DISPLAY.

loose = strcmp(get(0, 'FormatSpacing'), 'loose');

% Get display style and remove trivial empty TREEFUN2 case.
if ( isempty(f) )
    fprintf('    empty treefun2\n')
    if ( loose )
        fprintf('\n');
    end
    return
end

% Get information that we want to display:
dom = f.domain(:,1); % Domain
nleaves = length(f); % Number of leaves
nlevels = height(f); % Number of levels

% Display the information:
disp('   treefun2 object')
fprintf('       domain                order   length   levels\n');
fprintf('[%4.2g,%4.2g] x [%4.2g,%4.2g]    %5i   %6i   %6i\n', ...
    dom(1), dom(2), dom(3), dom(4), f.n, nleaves, nlevels);

if ( loose )
    fprintf('\n');
end

end
