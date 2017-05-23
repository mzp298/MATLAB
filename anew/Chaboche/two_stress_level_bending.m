clear;clc;close all;
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
%1stOpt when y=3 ; t=pi sqrj reaches sqrjmax(t)
y=3 ;
t=pi ; 
sig1 = [0 0 0;0 0 0;0 0 -F1*L*y/I*sin(t+pi/2);] 
sig2 = [0 0 0;0 0 0;0 0 -F2*L*y/I*sin(t+pi/2);]

p1=1/3*sum(diag(sig1));
p2=1/3*sum(diag(sig2));

S1=2*sig1-(1/3*sum(diag(sig1)))*diag([1,1,1]);
S2=2*sig2-(1/3*sum(diag(sig2)))*diag([1,1,1]);

sqrj1=1/2*sqrt(1/2*(S1(1,1)^2+S1(2,2)^2+S1(3,3)^2+2*(S1(1,2)^2)+...
2*(S1(1,3)^2)+2*(S1(2,3)^2)))
sqrj2=1/2*sqrt(1/2*(S2(1,1)^2+S2(2,2)^2+S2(3,3)^2+2*(S2(1,2)^2)+...
2*(S2(1,3)^2)+2*(S2(2,3)^2)))

pm1=p1
pm2=p2

b=10*ac/(3*ff);%we set parameter '3*b*ff' 10 times as 'ac' in Crossland

%---------------Modification with Crossland and DangVan-------------
cross1=sqrj1+ac*pm1-bc
cross2=sqrj2+ac*pm2-bc
tau1= 0.5*(sig1(3,3)-0);
dangvan1=tau1+ad*pm1-bd;
tau2= 0.5*(sig2(3,3)-0);
dangvan2=tau2+ad*pm2-bd;

%---------------High-Low sequence-------------
eta1s=(sqrt(4/3)*sqrj2-ff*(1-3*b*pm2))/(sqrt(4/3)*sqrj1-ff*(1-3*b*pm1))*...
(sigu-2*sqrj1)/(sigu-2*sqrj2)
eta1c=(cross2)/(cross1)*...
(sigu-2*sqrj1)/(sigu-2*sqrj2)
eta1d=(dangvan2)/(dangvan1)*...
(sigu-2*sqrj1)/(sigu-2*sqrj2)
%---------------Low-High sequence-------------
eta2s=(sqrt(4/3)*sqrj1-ff*(1-3*b*pm1))/(sqrt(4/3)*sqrj2-ff*(1-3*b*pm2))*...
(sigu-2*sqrj2)/(sigu-2*sqrj1)
eta2c=(cross1)/(cross2)*...
(sigu-2*sqrj2)/(sigu-2*sqrj1)
eta2d=(dangvan1)/(dangvan2)*...
(sigu-2*sqrj2)/(sigu-2*sqrj1)
%---------------plot life ratio-------------
N1NF1=0:0.01:1;
N2NF2s=1-N1NF1.^eta1s;
N2NF22s=1-N1NF1.^eta2s;
N2NF2c=1-N1NF1.^eta1c;
N2NF22c=1-N1NF1.^eta2c;
N2NF2d=1-N1NF1.^eta1d;
N2NF22d=1-N1NF1.^eta2d;

PAs2=plot(N1NF1,N2NF2s,'g','LineStyle', ':','LineWidth', 7);
hold on
PAs22=plot(N1NF1,N2NF22s,'g','LineStyle', ':','LineWidth', 7);
PAc2=plot(N1NF1,N2NF2c,'r','LineStyle', '-','LineWidth', 7);
PAc22=plot(N1NF1,N2NF22c,'r','LineStyle', '-','LineWidth', 7);
PAd2=plot(N1NF1,N2NF2d,'b','LineStyle', '--','LineWidth', 7);
PAd22=plot(N1NF1,N2NF22d,'b','LineStyle', '--','LineWidth', 7);

axis equal;
axis([0 1 0 1]);
grid on;
grid minor;
set(gca ,'FontSize',30);
hXLabel = xlabel('N_{1}/N_{F1}' ,'Fontsize' ,30);
hYLabel = ylabel('N_{2}/N_{F2}', 'Fontsize' ,30);
hTitle = title('Two-stress level loading in 4-point bending at edge of bar' ,'Fontsize' ,30);
set(hTitle, 'FontSize', 30, 'FontWeight' , 'bold')
hLegend=legend([PAs2,PAc2,PAd2],'Sines criteria','Crossland criteria','Dangvan criteria',...
    'Location','bestoutside');
set(hLegend, 'FontSize', 20)
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white
% Adjust font
set(gca, 'FontName', 'Helvetica')
set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')
% Adjust axes properties
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
% Maximize print figure
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1280 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])


saveas(gcf,'F:\Git\PhDreport\5thesis\figures\2stressB.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('Job done');
