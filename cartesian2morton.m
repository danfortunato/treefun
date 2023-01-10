function morton = cartesian2morton(x, y)

morton = bitshift(Part1By1_64(y), 1) + Part1By1_64(x);

end

function x = Part1By1_64(x)
% "Insert" a 0 bit after each of the 16 low bits of x

x = bitand(x, 0x00000000ffffffff);
x = bitand(bitxor(x, bitshift(x, 16)), 0x0000ffff0000ffff);
x = bitand(bitxor(x, bitshift(x, 8)),  0x00ff00ff00ff00ff);
x = bitand(bitxor(x, bitshift(x, 4)),  0x0f0f0f0f0f0f0f0f);
x = bitand(bitxor(x, bitshift(x, 2)),  0x3333333333333333);
x = bitand(bitxor(x, bitshift(x, 1)),  0x5555555555555555);

end
