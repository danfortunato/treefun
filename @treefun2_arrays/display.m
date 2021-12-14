function display(X) %#ok<DISPLAY>
%DISPLAY   Display information about a TREEFUN2.
%   DISPLAY(F) outputs important information about the TREEFUN2 F to the
%   command window, including its domain of definition, the number of boxes
%   used to represent it, and a summary of its structure.
%
%   It is called automatically when the semicolon is not used at the end of
%   a statement that results in a TREEFUN2.
%
% See also DISP.

if ( isequal(get(0, 'FormatSpacing'), 'compact') )
	disp([inputname(1), ' =']);
	disp(X);
else
	disp(' ');
	disp([inputname(1), ' =']);
	disp(' ');
    disp(X);
end

end
