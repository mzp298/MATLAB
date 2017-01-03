%*****************Fitting********************
clear;clc;
%***************** steel 30NCD16 data from Jabbado thesis*****************
%***************** bending*****************
stress11=[8.2E8 7.95E8 7.9E8 7.85E8 7.8E8 7.65E8 7.52E8 7.25E8 7.20E8 7.15E8 7.08E8];
stress12=zeros(1,length(stress11));
stress13=zeros(1,length(stress11));
stress21=zeros(1,length(stress11));
stress22=zeros(1,length(stress11));
stress23=zeros(1,length(stress11));
stress31=zeros(1,length(stress11));
stress32=zeros(1,length(stress11));
stress33=zeros(1,length(stress11));

for n=1:length(stress11)
hydro=1/3*sum(stress11(n)+stress22(n)+stress33(n));
% M=(1-3*hydro/sigu); %M function in Chaboche model
% M(M<0)=0;
% yield(n)=y*M; %macro yield strength considering mean stress effect
dev=[stress11(n) stress12(n) stress13(n);stress21(n) stress22(n) stress23(n);stress31(n) stress32(n) stress33(n)]-hydro*eye(3);
normdevben(n)=sqrt(sum(dev(:).^2));
n=n+1;
end
Smaxben=normdevben;
Nfben=[51000 80000 90000 95000 100000 120000 140000 200000 210000 230000 250000] ;


k=1e9;
e=191e9;
%**************************************************** fitted parameters**************************************
sigmay=1080e6;
bben=8.716;
Wfben=1.657e8;
btor=7.043;
Wftor=1.243e9;

Smax=4e8:1000:sigmay;
Nf=Wfben*(4*(191e9-1e9)*Smax.^(bben+1)/(191e9*(191e9+1e9*0.38)*bben*(bben+1)*1080e6^(bben-1))).^-1;

figure(1)
fittingben=plot(Nf,Smax,'-b','LineWidth',4);
axis([1e4 3e5 4e8 sigmay]);
hold on;
expben=plot(Nfben,Smaxben,'d','MarkerSize',15, 'MarkerFaceColor','b');
%***************** torsion*****************
stress12=[5.27e8 5.05e8 5e8 4.97e8 4.95e8 4.82e8 4.7e8 4.5e8 4.46e8 4.45e8 4.4e8];
stress21=[5.27e8 5.05e8 5e8 4.97e8 4.95e8 4.82e8 4.7e8 4.5e8 4.46e8 4.45e8 4.4e8];
stress11=zeros(1,length(stress12));
stress13=zeros(1,length(stress12));
stress22=zeros(1,length(stress12));
stress23=zeros(1,length(stress12));
stress31=zeros(1,length(stress12));
stress32=zeros(1,length(stress12));
stress33=zeros(1,length(stress12));

for n=1:length(stress12)
hydro=1/3*sum(stress11(n)+stress22(n)+stress33(n));
% M=(1-3*hydro/sigu); %M function in Chaboche model
% M(M<0)=0;
% yield(n)=y*M; %macro yield strength considering mean stress effect
dev=[stress11(n) stress12(n) stress13(n);stress21(n) stress22(n) stress23(n);stress31(n) stress32(n) stress33(n)]-hydro*eye(3);
normdevtor(n)=sqrt(sum(dev(:).^2));
n=n+1;
end
Smaxtor=normdevtor;
Nftor=[51000 80000 90000 95000 100000 120000 140000 200000 210000 230000 250000];

Nf=Wftor*(4*(191e9-1e9)*Smax.^(btor+1)/(191e9*(191e9+1e9*0.38)*btor*(btor+1)*1080e6^(btor-1))).^-1;
hold on;
fittingtor=plot(Nf,Smax,'--r','LineWidth',4);
exptor=plot(Nftor,Smaxtor,'^','MarkerSize',15, 'MarkerFaceColor','r');

grid on;
grid minor;
set(gca,'FontSize',30);
hLegend=legend('Bending fitted curve','Bending data','Torsion fitted curve','Torsion data');
set([hLegend, gca], 'FontSize', 30);
hTitle = title('Bending and torsion tests of 30NCD16 steel','Fontsize',30);
hXLabel = xlabel('N_f','Fontsize',30);
hYLabel = ylabel('S_{max}','Fontsize',30);
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
set(gcf, 'PaperPosition', [0 0 1920 1080]); %set(gcf,'PaperPosition',[left,bottom,width,height])
% saveas(gcf,'F:\Git\Anew\figures\30NCD60noH.png');


%*****************Fitting********************
cftool
Smaxben=normdevben;
Nfben=[51000 80000 90000 95000 100000 120000 140000 200000 210000 230000 250000] ;
cftool
Nf=Wfben*(4*(191e9-1e9)*Smax.^(bben+1)/(191e9*(191e9+1e9*0.38)*bben*(bben+1)*1080e6^(bben-1))).^-1;

clear;clc;
Smaxtor=normdevtor;
Nftor=[51000 80000 90000 95000 100000 120000 140000 200000 210000 230000 250000] ;
cftool
Nf=Wftor*(4*(191e9-1e9)*Smax.^(btor+1)/(191e9*(191e9+1e9*0.38)*btor*(btor+1)*1080e6^(btor-1))).^-1;

Wf=1e9~5e11
b=1~50

