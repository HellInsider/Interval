# kinterval demo

x = kinterval(10,6)
y = kinterval(-2,1)

M = [x y; 2*x y/2]

z1 = x + y
z2 = x - y
z3 = x * y
z4 = y / x
z5 = x / y

M1 = M * M';
M2 = M .* M';
M3 = dual(M2);
