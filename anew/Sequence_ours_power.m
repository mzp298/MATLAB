clear;clc
data=xlsread('F:\stress values 4points bending');
F1=1e6;
F2=0.8e6;
L=data(3,1);
tt=data(15,1);
ff=data(17,1);
I=data(11,1);
ac=(tt-ff/sqrt(3))/(ff/3);
ad=(3*tt)/(ff)-(1.5);
ap=3*(tt/ff-1/2);
bc=tt;
bd=tt;
bp=tt;
sigu=1.67e6;
r=0.1;
y=3 ;
t=pi ; 
sig1 = [0 0 0;0 0 0;0 0 -F1*L*y/I*sin(t+pi/2);] ;
sig2 = [0 0 0;0 0 0;0 0 -F2*L*y/I*sin(t+pi/2);] ;

pm1=1/3*sum(diag(sig1));
pm2=1/3*sum(diag(sig2));

S1=2*sig1-(1/3*sum(diag(sig1)))*diag([1,1,1]);
S2=2*sig2-(1/3*sum(diag(sig2)))*diag([1,1,1]);

sqrj1=1/2*sqrt(1/2*(S1(1,1)^2+S1(2,2)^2+S1(3,3)^2+2*(S1(1,2)^2)+...
2*(S1(1,3)^2)+2*(S1(2,3)^2)));
sqrj2=1/2*sqrt(1/2*(S2(1,1)^2+S2(2,2)^2+S2(3,3)^2+2*(S2(1,2)^2)+...
2*(S2(1,3)^2)+2*(S2(2,3)^2)));

b=10*ac/(3*ff);%we set parameter '3*b*ff' 10 times as 'ac' in Crossland
sigy=1.3e6;
lam=0.3;
yield1=sigy-lam*pm1;
yield2=sigy-lam*pm2;

%---------------High-Low sequence-------------
X=1.065;
fb=1.4;
eta1c=((sqrj1/yield1/(X-sqrj1/yield1))/(sqrj2/yield2/(X-sqrj2/yield2))).^fb
%---------------Low-High sequence-------------
eta2c=((sqrj2/yield2/(X-sqrj2/yield2))/(sqrj1/yield1/(X-sqrj1/yield1))).^fb

%---------------High-Low sequence-------------
fb=3;
eta1s=((sqrj1/yield1/(X-sqrj1/yield1))/(sqrj2/yield2/(X-sqrj2/yield2))).^fb
%---------------Low-High sequence-------------
eta2s=((sqrj2/yield2/(X-sqrj2/yield2))/(sqrj1/yield1/(X-sqrj1/yield1))).^fb

%---------------plot life ratio-------------
N1NF1=0:0.001:1;
N2NF2c=1-N1NF1.^eta1c;
N2NF22c=1-N1NF1.^eta2c;
N2NF2s=1-N1NF1.^eta1s;
N2NF22s=1-N1NF1.^eta2s;
N2NFm=1-N1NF1;
hold on
figure(1);
sequence=plot(N1NF1,N2NF22c,'-.c',N1NF1,N2NF2c,'--g',N1NF1,N2NF22s,'db',N1NF1,N2NF2s,'om',N1NF1,N2NFm,'k','LineWidth',3);
axis equal;
axis([0 1 0 1]);
grid on;
hTitle =title('Two-stress level loading lifetime proportion');
hXLabel = xlabel('T_{1}/T_{F1}');
hYLabel =ylabel('T_{2}/T_{F2}');
set(gca,'XTick',0:0.1:1);
set(gca,'YTick',0:0.1:1);
hLegend=legend( 'Low-High loading sequence(magnification power=1.4)','High-Low loading sequence(magnification power=1.4)',...
    'Low-High loading sequence(magnification power=3)','High-Low loading sequence(magnification power=3)','Miner law','Location','Bestoutside');
set([hLegend, gca], 'FontSize', 20)
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white 
% Adjust font
set(gca, 'FontName', 'Helvetica')
set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hXLabel, hYLabel], 'FontSize',30)
set(hTitle, 'FontSize', 30, 'FontWeight' , 'bold')
% Adjust axes properties
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
% Maximize print figure
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1080 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])
%  saveas(gcf,'F:\Git\Anew\figures\sequence_ours.png');

