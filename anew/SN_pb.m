clear;clc;
% close all;
%---------------------Verified parameters in random loading case-----------------------------
b1=1.1;                    %weakening scales distribution exponent
b2=1.1;                    %weakening scales distribution exponent
pb1=1.1;
pb2=2;
a1=0.1;                %sensitivity of sequence effect(control alp>0)
a2=0.1;                %sensitivity of sequence effect(control alp>0)
W0=3.27e8;            %dissipated energy to failure per unit volume
%---------------------Verified parameters in constant loading case-----------------------------
y=230e6;           %macroscopic yield stress
E=72e9;              %Young's modulus
k=6e8;                 %hardening parameter
nu=0.3;                     %poisson's ratio
sigu=320e6;             %ultimite stress
lam=0.1;               %hydrostatic pressure sensitivity

%---------------a b correspond to high and low bounds of cyclic load, 1 2
%correspond to  different parameters comparison------------------
stress11=1000:1000:y;
m=0.1*stress11;
Smax=sqrt(2/3).*stress11;
hydroa=1/3.*(stress11+m);
yielda=y-lam.*hydroa;
smina=yielda.*Smax.^-1;
hydrob=1/3.*(-stress11+m);
yieldb=y-lam.*hydrob;
sminb=yieldb.*Smax.^-1;
yieldm=(yielda+yieldb)/2;
alphamax1=(1-a1.*(smina-1).^-pb1);
alphamin1=(1-a1.*(sminb-1).^-pb1);
alphamax2=(1-a2.*(smina-1).^-pb2);
alphamin2=(1-a2.*(sminb-1).^-pb2);
alpm1=(alphamax1+alphamin1)/2;
alpm2=(alphamax2+alphamin2)/2;

NF1=(1-alpm1).^-1*W0*E*(E+k*nu)*b1*(b1+1)*((4*(E-k)*(1+nu)*(b1-1)))^-1.*yieldm.^(b1-1).*Smax.^(-b1-1);
NF2=(1-alpm2).^-1*W0*E*(E+k*nu)*b2*(b2+1)*((4*(E-k)*(1+nu)*(b2-1)))^-1.*yieldm.^(b2-1).*Smax.^(-b2-1);


sn1=semilogx(NF1,stress11,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 10, ...
            'MarkerEdgeColor', [238 99 99]/255, 'MarkerFaceColor',[238 99 99]/255);
hold on
sn2=semilogx(NF2,stress11,'LineStyle', 'none','LineWidth', 1,'Marker', 's', 'MarkerSize', 10, ...
            'MarkerEdgeColor', [153 50 204]/255, 'MarkerFaceColor',[153 50 204]/255);
grid on;
grid minor;
xlim([1e3 1e9]);
set(gca,'XTick',[1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9] );
% hTitle =title('S-N curve at constant amplitude stress with R=-1 in our model with different $\beta$','FontSize',20);
hXLabel = xlabel('N_F');
hYLabel =ylabel('Stress');

hLegend=legend([sn1,sn2], 'f(\beta)=1.1','f(\beta)=2','Location','northeast');
set([hLegend, gca], 'FontSize', 40);
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
%  saveas(gcf,'F:\Git\Anew\figures\SNpb.png');