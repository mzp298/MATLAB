clear;clc;close all;
%---------------------Verified parameters in random loading case-----------------------------
b1=1.5;                    %weakening scales distribution exponent
b2=1.5;                    %weakening scales distribution exponent
a=0.1;                %sensitivity of sequence effect(control alp>0)
WF=5.58e8;            %dissipated energy to failure per unit volume
X1=1.065;                 %Major damage threshold (meaning S_{max} must be greater than X/2*yield to activate magnification)
X2=1.065;                 %Major damage threshold (meaning S_{max} must be greater than X/2*yield to activate magnification)
pb1=1.4;                %Major damage power
pb2=3;                %Major damage power
%---------------------Verified parameters in constant loading case-----------------------------
y=230e6;           %macroscopic yield stress
E=72e9;              %Young's modulus
k=6e8;                 %hardening parameter
nu=0.3;                     %poisson's ratio
sigu=320e6;             %ultimite stress
n0=1;                   %number of initial local defects
lam=0.1;               %hydrostatic pressure sensitivity
gam=0.1;             %material parameter from Chaboche law(Wohler curve exponent)
maxstress11=sigu;
stress11=1000:1000:maxstress11;
Smax=sqrt(2/3).*stress11;
hydro=1/3.*stress11;
yield=y-lam.*hydro;

sequence1=((Smax.*yield.^-1).*(X1-Smax.*yield.^-1).^-1).^pb1;
sequence2=((Smax.*yield.^-1).*(X2-Smax.*yield.^-1).^-1).^pb2;

NF1=((1+gam)*a*sequence1).^-1*WF*E*(E+k*nu)*b1*(b1+1)*(n0*(4*(E-k)*(1+nu)*(b1-1)))^-1.*yield.^(b1-1).*Smax.^(-b1-1);
NF2=((1+gam)*a*sequence2).^-1*WF*E*(E+k*nu)*b2*(b2+1)*(n0*(4*(E-k)*(1+nu)*(b2-1)))^-1.*yield.^(b2-1).*Smax.^(-b2-1);

sn1=semilogx(NF1,stress11,'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 10, ...
            'MarkerEdgeColor', [238 99 99]/255, 'MarkerFaceColor',[238 99 99]/255);
hold on
sn2=semilogx(NF2,stress11,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 10, ...
            'MarkerEdgeColor', [153 50 204]/255, 'MarkerFaceColor',[153 50 204]/255);
grid on;
grid minor;
axis([1 1e8 0 maxstress11]);
set(gca,'XTick',[1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9] );
% hTitle =title('S-N curve at constant amplitude stress with R=-1 in our model with different $f(\beta)$');
hXLabel = xlabel('N_F');
hYLabel =ylabel('Stress');

hLegend=legend([sn1,sn2], 'f(\beta)=1.4','f(\beta)=3','Location','northeast');
set([hLegend, gca], 'FontSize', 40)
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white 
% Adjust font
set(gca, 'FontName', 'Helvetica')
set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hXLabel, hYLabel], 'FontSize',40)

% Adjust axes properties
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
% Maximize print figure
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1200 1000]); %set(gcf,'PaperPosition',[left,bottom,width,height])
 saveas(gcf,'F:\Git\Anew\figures\SNpb.png');