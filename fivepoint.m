clear;
clc;

sym E;
x = sym('x','real');
y = sym('y','real');
z = sym('z','real');

A = sym('A%d%d',[3 3]);
A = sym(A,'real');

B = sym('B%d%d',[3 3]);
B = sym(B,'real');

C = sym('C%d%d',[3 3]);
C = sym(C,'real');

D = sym('D%d%d',[3 3]);
D = sym(D,'real');

v = sym('v',[7 1]);
v = sym(v, 'real');

E = x*A+y*B+z*C+D;
EEE1 = expand(E*E'*E);
EEE2 = 0.5*trace(E*E')*eye(3);
EEE = expand(EEE1-EEE2);
EEE = EEE(:);
DE = det(E);
AAA = [EEE;DE];



AAA = subs(AAA,x^3,v(1));
AAA = subs(AAA,y^3,v(2));
AAA = subs(AAA,x^2*y,v(3));
AAA = subs(AAA,x*y^2,v(4));
AAA = subs(AAA,x^2,v(5));
AAA = subs(AAA,y^2,v(6));
AAA = subs(AAA,x*y,v(7));

coeffA = sym('coeffA',[10 10]);
coeffA = sym(coeffA, 'real');
for i = 1:10
    [coeffA(i,:),tx]= coeffs(AAA(i),[v(1:7)',x,y]);
end
cz = det(coeffA);
[coeffz tz] = coeffs(cz,z);
save AAA AAA;
save coeffA coeffA;
save tx tx;








