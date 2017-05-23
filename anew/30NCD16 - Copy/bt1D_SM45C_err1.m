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
m=0;%mean tension

%---------------a b correspond to high and low bounds of cyclic load, 1 2
%correspond to  different parameters comparison------------------
NFben=[1.9E4 3.6E4 5.8E4 9E4 1.7E5 2.2E5 4.4E5 8E5 1.4E6 1.5E6 2.4E6 3.3E6 6E6 6.6E6 7E6 9E6] ;
stressben=[6.2E8 5.9E8 5.52E8 5.35E8 5.05E8 4.9E8 4.7E8 4.65E8 4.62E8 4.65E8 4.6E8 4.6E8 4.4E8 4.58E8 4.3E8 4.48E8];%to get Smaxben
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
n=n+1;
end

%-------------------------------------------------torsion--------------------------------
stresstor=[4.05E8 3.99E8 3.80E8 3.65E8 3.55E8 3.49E8 3.40E8 3.35E8 3.33E8 3.25E8];%to get Smaxtor
NFtor=[2.9E4 5E4 7.8E4 1E5 1.8E5 1.9E5 3E5 5.8E5 8.2E5 2.25E6];
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
axis([1e4 1e7 1e4 1e7]);
set(gca,'xtick',[1e4 1e5 1e6 1e7]); 
set(gca,'ytick',[1e4 1e5 1e6 1e7]); 
hLegend=legend([err_ben,err_tor,err_tor_VM,err_tor_tresca],...
    'Bending test on SM45(R=-1)','Torsion test on SM45(R=-1)',...
    'Torsion test on SM45 with Von Mises stress(R=-1)','Torsion test on SM45 with Tresca stress(R=-1)','location','bestoutside');
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
saveas(gcf,'F:\Git\Anew\figures\bt1D_SM45C_sn.png');
figure(2);
saveas(gcf,'F:\Git\Anew\figures\bt1D_SM45C_err1.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('done done done');

% figure(3);
% plot(Smaxben,'b')
% hold on 
% plot(Smaxtor,'r')
