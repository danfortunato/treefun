function [x, y] = morton2cartesian_32(morton)

x = Compact1By1_32(morton);
y = Compact1By1_32(bitshift(morton, -1));

end

% Inverse of Part1By1 - "delete" all odd-indexed bits
function x = Compact1By1_32(x)

x = bitand(x, 0x55555555);
x = bitand(bitxor(x, bitshift(x, -1)), 0x33333333);
x = bitand(bitxor(x, bitshift(x, -2)), 0x0f0f0f0f);
x = bitand(bitxor(x, bitshift(x, -4)), 0x00ff00ff);
x = bitand(bitxor(x, bitshift(x, -8)), 0x0000ffff);

end
