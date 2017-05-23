clear;clc
%---------------------Verified parameters in random loading case-----------------------------
b1=1.5;                    %weakening scales distribution exponent
b2=5;                    %weakening scales distribution exponent
a=0.1;                %sensitivity of sequence effect(control alp>0)
WF=5.58e5;            %dissipated energy to failure per unit volume
X=1.065;                 %Major damage threshold (meaning S_{max} must be greater than X/2*yield to activate magnification)
pb=1.4;                %Major damage power
%---------------------Verified parameters in constant loading case-----------------------------
y=230e6;           %macroscopic yield stress
E=72e9;              %Young's modulus
k=6e8;                 %hardening parameter
nu=0.3;                     %poisson's ratio
sigu=320e6;             %ultimite stress
n0=1;                   %number of initial local defects
lam=0.3;               %hydrostatic pressure sensitivity
gam=0.1;             %material parameter from Chaboche law(Wohler curve exponent)
maxstress11=190e5;
stress11=0:1000:maxstress11;
Smax=sqrt(2/3).*stress11;
hydro=1/3.*stress11;
yield=y-lam.*hydro;

sequence=((Smax.*yield.^-1).*(X-Smax.*yield.^-1).^-1).^pb;
alp=1-a.*sequence;

NF1=((1+gam)*a*sequence).^-1*WF*E*(E+k*nu)*b1*(b1+1)*(n0*(4*(E-k)*(1+nu)*(b1-1)))^-1.*yield.^(b1-1).*Smax.^(-b1-1);
NF2=((1+gam)*a*sequence).^-1*WF*E*(E+k*nu)*b2*(b2+1)*(n0*(4*(E-k)*(1+nu)*(b2-1)))^-1.*yield.^(b2-1).*Smax.^(-b2-1);
hold on;
sn1=semilogx(NF1,Smax,'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 10, ...
            'MarkerEdgeColor', [238 99 99]/255, 'MarkerFaceColor',[238 99 99]/255);
sn2=semilogx(NF2,Smax,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 10, ...
            'MarkerEdgeColor', [153 50 204]/255, 'MarkerFaceColor',[153 50 204]/255);
