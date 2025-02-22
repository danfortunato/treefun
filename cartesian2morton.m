function morton = cartesian2morton(x, y, z)

if nargin == 2
  morton = bitshift(Part1By1_64(y), 1) + Part1By1_64(x);
end

if nargin == 3
  morton = bitshift(Part1By2_64(z), 2) + bitshift(Part1By2_64(y), 1) + Part1By2_64(x);
end

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

function x = Part1By2_64(x)
% "Insert" two 0 bits after each of the 16 low bits of x
% https://github.com/trevorprater/pymorton/blob/f248683c19d90e904193c58bbd03e77bc2c43768/pymorton/pymorton.py#L20

x = bitand(x, uint64(0x1fffff));
x = bitand(bitxor(x, bitshift(x, 32)), uint64(0x1f00000000ffff));
x = bitand(bitxor(x, bitshift(x, 16)), uint64(0x1f0000ff0000ff));
x = bitand(bitxor(x, bitshift(x, 8)),  uint64(0x100f00f00f00f00f));
x = bitand(bitxor(x, bitshift(x, 4)),  uint64(0x10c30c30c30c30c3));
x = bitand(bitxor(x, bitshift(x, 2)),  uint64(0x1249249249249249));

end


