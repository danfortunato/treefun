function morton = cartesian2morton_32(x, y)

morton = bitshift(Part1By1_32(y), 1) + Part1By1_32(x);

end

function x = Part1By1_32(x)
% "Insert" a 0 bit after each of the 16 low bits of x

x = bitand(x, 0x0000ffff);
x = bitand(bitxor(x, bitshift(x, 8)), 0x00ff00ff);
x = bitand(bitxor(x, bitshift(x, 4)), 0x0f0f0f0f);
x = bitand(bitxor(x, bitshift(x, 2)), 0x33333333);
x = bitand(bitxor(x, bitshift(x, 1)), 0x55555555);

end
