function disp(f)
%DISP   Display a TREEFUN to the command line.
% 
% See also DISPLAY.

loose = strcmp(get(0, 'FormatSpacing'), 'loose');

% Get display style and remove trivial empty TREEFUN case.
if ( isempty(f) )
    fprintf('    empty treefun\n')
    if ( loose )
        fprintf('\n');
    end
    return
end

% Get information that we want to display:
dom = f.domain;      % Domain
nleaves = length(f); % Number of leaves
nlevels = height(f); % Number of levels

% Display the information:
disp('   treefun object')
fprintf('       domain          order   length   levels\n');
fprintf('[%8.2g,%8.2g]    %5i   %6i   %6i\n', ...
    dom(1), dom(2), f.n-1, nleaves, nlevels);

if ( loose )
    fprintf('\n');
end

end
