clc;
clear;
close all;
% steel 30 NCD 16 data from Jabbado thesis
%-------------------------------------------------bending--------------------------------
y=1080e6;
k=1e9;
E=191e9;
nu=0.3;                     %poisson's ratio
a=0.1;
W0=44e9;
lam=0.85;
b=1.85;
fb=b;
m=0;%mean tension
%---------------a b correspond to high and low bounds of cyclic load, 1 2
%correspond to  different parameters comparison------------------
NFben=[51000 80000 90000 95000 100000 120000 140000 200000 210000 230000 250000] ;
stressben=[8.2e8 7.95e8 7.9e8 7.85e8 7.8e8 7.65e8 7.52e8 7.25e8 7.20e8 7.15e8 7.08e8];%to get Smaxben
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

%-------------------------------------------------torsion--------------------------------
NFtor=NFben;
stresstor=[5.27e8 5.05e8 5e8 4.97e8 4.95e8 4.82e8 4.7e8 4.5e8 4.46e8 4.45e8 4.4e8];%to get Smaxtor
n=1;
while n<length(stresstor)+1
hydro=1/3*sum(m +0+0);
dev1=[0 stresstor(n) 0 ;0 stresstor(n) 0 ;0 0 0 ]-hydro*eye(3);
Smaxtor(n)=norm(dev1,'fro'); % deviatoric stress
Smaxtor_VM(n)=VM_Stress(0, 0, 0, stresstor(n), 0, 0);
Smaxtor_tresca(n)=stresstor(n);
%--------to get average smin--------
hydromax=1/3.*(m);
yieldmin=y-lam.*hydromax;
sminmin=yieldmin.*Smaxtor(n).^-1;

hydromin=1/3.*(m);
yieldmax=y-lam.*hydromin;
sminmax=yieldmax.*Smaxtor(n).^-1;

yieldm=(yieldmin+yieldmax)/2;%-------------------------yield line should be lower in torsion to fit---------------
alphamax=(1-a.*(sminmax-1).^-fb);
alphamin=(1-a.*(sminmin-1).^-fb);
alpm_tor(n)=(alphamax+alphamin)/2;
%--------use average smin to get NF--------
NFtor_num(n)=(1-alpm_tor(n)).^-1*W0*E*(E+k*nu)*b*(b+1)*((4*(E-k)*(1+nu)*(b-1)))^-1.*yieldm.^(b-1).*Smaxtor(n).^(-b-1);
NFtor_num_VM(n)=(1-alpm_tor(n)).^-1*W0*E*(E+k*nu)*b*(b+1)*((4*(E-k)*(1+nu)*(b-1)))^-1.*yieldm.^(b-1).*Smaxtor_VM(n).^(-b-1);
NFtor_num_tresca(n)=(1-alpm_tor(n)).^-1*W0*E*(E+k*nu)*b*(b+1)*((4*(E-k)*(1+nu)*(b-1)))^-1.*yieldm.^(b-1).*Smaxtor_tresca(n).^(-b-1);
n=n+1;
end

figure(1);
ben_exp=semilogx(NFben,Smaxben,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[208 32 144]/255, 'MarkerFaceColor','none');
hold on;
tor_exp=semilogx(NFtor,Smaxtor,'v','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor','b', 'MarkerFaceColor','none');
NFben_num=(1-alpm_ben).^-1*W0*E*(E+k*nu)*b*(b+1)*((4*(E-k)*(1+nu)*(b-1)))^-1.*yieldm.^(b-1).*Smaxben.^(-b-1);
NFtor_num=(1-alpm_tor).^-1*W0*E*(E+k*nu)*b*(b+1)*((4*(E-k)*(1+nu)*(b-1)))^-1.*yieldm.^(b-1).*Smaxtor.^(-b-1);
ben_num=semilogx(NFben_num,Smaxben,'-r','LineWidth', 3);
tor_num=semilogx(NFtor_num,Smaxtor,'--b','LineWidth', 3);
set(gca ,'FontSize',30);
hXLabel = xlabel('NF_{num}','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('S_{max}','Fontsize',30, 'FontWeight' , 'bold');
hLegend=legend([ben_exp,ben_num,tor_exp,tor_num],...
    'Bending experimental result','Bending numerical result',...
    'Torsion experimental result','Torsion numerical result','location','bestoutside');
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
saveas(gcf,'F:\Git\Anew\figures\bt1D_30NCD16_sn.png');
%% 
figure(2);
err_ben = loglog(NFben,NFben_num,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[208 32 144]/255, 'MarkerFaceColor','none');
hold on;
err_tor = loglog(NFtor,NFtor_num,'v','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor','b', 'MarkerFaceColor','none');
err_tor_VM = loglog(NFtor,NFtor_num_VM,'^','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor','c', 'MarkerFaceColor','none');
err_tor_tresca = loglog(NFtor,NFtor_num_tresca,'s','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor','g', 'MarkerFaceColor','none');
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
hLegend=legend([err_ben,err_tor,err_tor_VM,err_tor_tresca],...
    'Bending test on 30NCD16(R=-1)','Torsion test on 30NCD16(R=-1)',...
    'Torsion test on 30NCD16 with Von Mises stress(R=-1)','Torsion test on 30NCD16 with Tresca stress(R=-1)','location','bestoutside');
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
saveas(gcf,'F:\Git\Anew\figures\bt1D_30NCD16_err1.png');
saveas(gcf,'F:\Git\PhDreport\3Plastisity\figures\bt1D_30NCD16_err1.png');
sp=actxserver('SAPI.SpVoice');
sp.Speak('done done done');
figure(2);
plot(alpm_ben,'b')
hold on 
plot(alpm_tor,'r')
% figure(3);
% plot(Smaxben,'b')
% hold on 
% plot(Smaxtor,'r')
