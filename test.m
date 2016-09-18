clear;
clc;
randn('seed',0);
R = [0.4120    0.1411    0.9002;
    0.7969   -0.5349   -0.2808;
    0.4419    0.8330   -0.3328];
tran = [1 2 3]';
Etruth = [0 -3 2;
          3 0 -1;
          -2 1 0]*R;
Etruth = Etruth/Etruth(3,3);


p1a = randn(3,1);
p1b = R*p1a+tran;
p1a = p1a/p1a(3);
p1b = p1b/p1b(3);

p2a = randn(3,1);
p2b = R*p2a+tran;
p2a = p2a/p2a(3);
p2b = p2b/p2b(3);

p3a = randn(3,1);
p3b = R*p3a+tran;
p3a = p3a/p3a(3);
p3b = p3b/p3b(3);

p4a = randn(3,1);
p4b = R*p4a+tran;
p4a = p4a/p4a(3);
p4b = p4b/p4b(3);

p5a = randn(3,1);
p5b = R*p5a+tran;
p5a = p5a/p5a(3);
p5b = p5b/p5b(3);

p6a = randn(3,1);
p6b = R*p6a+tran;
p6a = p6a/p6a(3);
p6b = p6b/p6b(3);


p7a = randn(3,1);
p7b = R*p7a+tran;
p7a = p7a/p7a(3);
p7b = p7b/p7b(3);

p8a = randn(3,1);
p8b = R*p8a+tran;
p8a = p8a/p8a(3);
p8b = p8b/p8b(3);
% p1b'*Etruth*p1a
% p2b'*Etruth*p2a
% p3b'*Etruth*p3a
% p4b'*Etruth*p4a
% p5b'*Etruth*p5a
A = zeros(5,9);
A(1,:) = reshape(p1a*p1b',[1 9]);
A(2,:) = reshape(p2a*p2b',[1 9]);
A(3,:) = reshape(p3a*p3b',[1 9]);
A(4,:) = reshape(p4a*p4b',[1 9]);
A(5,:) = reshape(p5a*p5b',[1 9]);
[U,S,V] = svd(A);
E1 = reshape(V(:,6),[3 3])';
E2 = reshape(V(:,7),[3 3])';
E3 = reshape(V(:,8),[3 3])';
E4 = reshape(V(:,9),[3 3])';

x = sym('x','real');
y = sym('y','real');
z = sym('z','real');

v = sym('v',[7 1]);
v = sym(v, 'real');

Esym = x*E1+y*E2+z*E3+E4;
EEE1 = 2*expand((Esym*Esym')*Esym);
EEE2 = trace(Esym*Esym')*Esym;
EEE = expand(EEE1-EEE2);
% EEE = EEE';
EEE = EEE(:);
DE = expand(det(Esym));
AAA = [EEE;DE];





coeffA = sym('coeffA',[10 10]);
coeffA = sym(coeffA, 'real');



for i = 1:10
    [coeffA(i,:),tx]= coeffs(AAA(i),[x,y]);
end
    
% for i = 1:10
%     for j = 1:10
%         [ccc,ddd]= coeffs(coeffA(i,j),z);
%         if size(ccc,2) == 1
%             C0(i,j) = ccc;
%         else if size(ccc,2) == 2
%                 C0(i,j) = ccc(2);
%                 C1(i,j) = ccc(1);
%             else if size(ccc,2) == 3
%                     C0(i,j) = ccc(3);
%                     C1(i,j) = ccc(2);
%                     C2(i,j) = ccc(1);
%                 else if size(ccc,2) == 4
%                         C0(i,j) = ccc(4);
%                         C1(i,j) = ccc(3);
%                         C2(i,j) = ccc(2);
%                         C3(i,j) = ccc(1);
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% [DDD FFF] = polyeig(C0,C1,C2,C3);
% index = find(FFF ~= Inf & imag(FFF) ==0);
% for i = 1:size(index,1)
%     Eestimate = DDD(6,index(i))/DDD(10,index(i))*E1+ DDD(9,index(i))/DDD(10,index(i))*E2+FFF(index(i))*E3+E4;
%     det(Eestimate)
%     Eestimate/Eestimate(3,3)-Etruth
%     
% end


cz = det(coeffA);
[coeffz tz] = coeffs(cz,z);
zroots = roots(coeffz);
xyz = zeros(4,3);

for i = 1:4
    Ai = subs(coeffA,zroots(i));
    [UU,SS,VV] = svd(Ai);
    xyz(i,:) = [VV(6,end)/VV(10,10),VV(9,end)/VV(10,10),zroots(i)];
end

for i = 1:4
    Eestimate = xyz(i,1)*E1+xyz(i,2)*E2+xyz(i,3)*E3+E4;
    display('determination:');
    det(Eestimate)
    display('sampson distance:');
    p1b'*Eestimate*p1a
    p2b'*Eestimate*p2a
    p3b'*Eestimate*p3a
    p4b'*Eestimate*p4a
    p5b'*Eestimate*p5a
    p6b'*Eestimate*p6a
    p7b'*Eestimate*p7a
    p8b'*Eestimate*p8a
    display('differrence');
    Eestimate/Eestimate(3,3)-Etruth
end



