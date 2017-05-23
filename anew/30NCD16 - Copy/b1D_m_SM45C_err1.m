clc;
clear;
close all;
% steel SM45C data from Jabbado thesis
%-------------------------------------------------bending--------------------------------
y=638e6;
k=1e9;
E=213e9;
nu=0.3;                     %poisson's ratio
a=0.1;
W0=1000e9;
lam=0.55;
b=3.8;
fb=b;
m=196e6;%mean tension

%---------------a b correspond to high and low bounds of cyclic load, 1 2
%correspond to  different parameters comparison------------------
Smax=3e8:1000:y;
NFben=[4.50E+04
5.40E+04
7.00E+04
7.10E+04
1.10E+05
1.50E+05
2.00E+05
2.10E+05
3.00E+05
4.30E+05
6.90E+05
5.80E+05] ;
stressben=1e6.*[540
515
520
485
485
475
460
455
435
415
410
390];%to get Smaxben

n=1;
while n<length(stressben)+1
hydro=1/3*sum(stressben(n)+m +0+0);
dev1=[stressben(n)+m 0 0 ;0 0 0 ;0 0 0 ]-hydro*eye(3);
Smaxben(n)=norm(dev1,'fro');
%--------to get average smin--------
hydromax=1/3.*(stressben(n)+m);
yieldmin=y-lam.*hydromax;
sminmin=yieldmin.*Smaxben(n).^-1;

hydromin=1/3.*(-stressben(n)+m);
yieldmax=y-lam.*hydromin;
sminmax=yieldmax.*Smaxben(n).^-1;

yieldm=(yieldmin+yieldmax)/2;
alphamax=(1-a.*(sminmax-1).^-fb);
alphamin=(1-a.*(sminmin-1).^-fb);
alpm_ben(n)=(alphamax+alphamin)/2;
%--------use average smin to get NF--------
NFben_num(n)=(1-alpm_ben(n)).^-1*W0*E*(E+k*nu)*b*(b+1)*((4*(E-k)*(1+nu)*(b-1)))^-1.*yieldm.^(b-1).*Smaxben(n).^(-b-1);
% NFben_num(n)=W0*E*(E+k*nu)*b*(b+1)*((4*(E-k)*(1+nu)*(b-1)))^-1.*yieldm.^(b-1).*Smaxben(n).^(-b-1)
n=n+1;
end

figure(1);
ben_exp=semilogx(NFben,Smaxben,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[208 32 144]/255, 'MarkerFaceColor','none');
NFben_num=(1-mean(alpm_ben)).^-1*W0*E*(E+k*nu)*b*(b+1)*((4*(E-k)*(1+nu)*(b-1)))^-1.*yieldm.^(b-1).*Smaxben.^(-b-1);
hold on;
ben_num=semilogx(NFben_num,Smaxben,'-','LineWidth', 3);
set(gca ,'FontSize',30);
hXLabel = xlabel('NF_{num}','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('S_{max}','Fontsize',30, 'FontWeight' , 'bold');
hLegend=legend([ben_exp,ben_num],...
    'Bending experimental result with \sigma_m=196 MPa','Bending numerical result with \sigma_m=196 MPa','location','bestoutside');
set(hLegend, 'FontSize', 18);
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white 
% Adjust font
set(gca, 'FontName', 'Helvetica')
% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid', 'on',...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
% Maximize print figure
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1800 1000]); %set(gcf,'PaperPosition',[left,bottom,width,height])

%% 
figure(2);
err_ben_m = loglog(NFben,NFben_num,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[208 32 144]/255, 'MarkerFaceColor','none');
hold on;
set(gca ,'FontSize',30);
hXLabel = xlabel('NF_{exp}','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('NF_{num}','Fontsize',30, 'FontWeight' , 'bold');
x=1e4:1e5:1e8;
y0=x;
hold on;
py0=loglog(x,y0,'k','LineWidth',3);
y1=2.*x;
y2=0.5.*x;
py1=loglog(x,y1, '--k','LineWidth',3);
py2=loglog(x,y2, '--k','LineWidth',3);
axis equal;
axis([1e4 1e6 1e4 1e6]);
set(gca,'xtick',[1e4 1e5 1e6]); 
set(gca,'ytick',[1e4 1e5 1e6]); 
hLegend=legend([err_ben_m],...
    'Bending test on SM45(R=-1,\sigma_m=196 MPa)','location','bestoutside');
set(hLegend, 'FontSize', 18);
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white 
% Adjust font
set(gca, 'FontName', 'Helvetica')
% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid', 'on',...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
% Maximize print figure
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1800 1000]); %set(gcf,'PaperPosition',[left,bottom,width,height])
figure(1);
saveas(gcf,'F:\Git\Anew\figures\b1D_m_SM45C_sn.png');
figure(2);
saveas(gcf,'F:\Git\Anew\figures\b1D_m_SM45C_err1.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('done done done');
