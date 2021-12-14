function out = treefunroot()
%TREEFUNROOT   Root directory of treefun installation.
%   S = TREEFUNROOT() returns a string that is the name of the directory
%   where treefun is installed.

out = fileparts(which('treefunroot'));

end
