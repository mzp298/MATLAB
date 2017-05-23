% Program to get the quadrature points
% and weight for Gauss-Legendre Quadrature
% Rule
clc
clear all
syms x
% Input n: Quad pt rule
n=30;
% Calculating the Pn(x)
% Legendre Polynomial
% Using recursive relationship
% P(order of polynomial, value of x)
% P(0,x)=1; P(1,x)=0;
% (i+1)*P(i+1,x)=(2*i+1)*x*P(i,x)-i*P(i-1,x)
m=n-1;
P0=1;
P1=x;
for i=1:1:m
    Pn=((2.0*i+1)*x*P1-i*P0)/(i+1.0);
    P0=P1;
    P1=Pn;
end
if n==1
    Pn=P1;
end
Pn=expand(Pn);
quadpts=solve(vpa(Pn,32));
quadpts=sort(quadpts);
% Finding the weights
% Formula for weights is given at
% http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
% Equation (13)
for k=1:1:n
    P0=1;
    P1=x;
    m=n;
    % Calculating P(n+1,x)
    for i=1:1:m
        Pn=((2.0*i+1)*x*P1-i*P0)/(i+1.0);
        P0=P1;
        P1=Pn;
    end
    Pn=P1;
    weights(k)=vpa(2*(1-quadpts(k)^2)/(n+1)^2/subs(Pn,x,quadpts(k))^2,32);
end
    fprintf('Quad point rule for n=%g \n',n)
disp('  ')
disp('Abscissas')
disp(quadpts)
disp('  ')
disp('Weights')
disp(weights')

quadpts=double(quadpts);
weights=double(weights');
txt={'x','Weight'};
xlswrite('Gauss-Legendre30.xlsx', txt(1),'Sheet1', 'A1') 
xlswrite('Gauss-Legendre30.xlsx', txt(2),'Sheet1', 'B1') 
xlswrite('Gauss-Legendre30.xlsx', quadpts,'Sheet1', 'A2:A31') 
xlswrite('Gauss-Legendre30.xlsx', weights, 'Sheet1', 'B2:B31') 
