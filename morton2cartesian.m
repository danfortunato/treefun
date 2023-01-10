function [x, y] = morton2cartesian(morton)

x = Compact1By1_64(morton);
y = Compact1By1_64(bitshift(morton, -1));

end

function x = Compact1By1_64(x)
% Inverse of Part1By1 - "delete" all odd-indexed bits

x = bitand(x, 0x5555555555555555);
x = bitand(bitxor(x, bitshift(x, -1)),  0x3333333333333333);
x = bitand(bitxor(x, bitshift(x, -2)),  0x0f0f0f0f0f0f0f0f);
x = bitand(bitxor(x, bitshift(x, -4)),  0x00ff00ff00ff00ff);
x = bitand(bitxor(x, bitshift(x, -8)),  0x0000ffff0000ffff);
x = bitand(bitxor(x, bitshift(x, -16)), 0x00000000ffffffff);
    
end
