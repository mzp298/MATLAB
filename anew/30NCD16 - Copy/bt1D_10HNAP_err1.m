clc;
clear;
close all;
% steel 10HNAP data from Jabbado thesis
%-------------------------------------------------bending--------------------------------
y=1080e6;
k=800e6;
E=215e9;
nu=0.3;                     %poisson's ratio
a=0.1;
W0=20e7;
lam=5.7;
b=2.8;
fb=b;
m=0;                     %mean tension
%---------------a b correspond to high and low bounds of cyclic load, 1 2
%correspond to  different parameters comparison------------------
NFben=[1.00E+05
2.00E+05
3.00E+05
4.00E+05
5.00E+05
6.00E+05
7.00E+05
8.00E+05
9.00E+05
1.00E+06] ;
stressben=1e6.*[326.69
304.42
292.11
283.68
277.30
272.20
267.96
264.34
261.19
258.41];%to get Smaxben
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
stresstor=1e6.*[260.70
242.70
232.17
224.70
218.90
214.16
210.16
206.69
203.63
200.90];%to get Smaxtor
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
    'Bending test on 10HNAP(R=-1)','Torsion test on 10HNAP(R=-1)',...
    'Torsion test on 10HNAP with Von Mises stress(R=-1)','Torsion test on 10HNAP with Tresca stress(R=-1)','location','bestoutside');
set(hLegend, 'FontSize', 18);
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white 
set(gca, 'FontName', 'Helvetica')
set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid', 'on',...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1800 1000]); %set(gcf,'PaperPosition',[left,bottom,width,height])
% figure(1);
% saveas(gcf,'F:\Git\Anew\figures\bt1D_10HNAP_sn.png');
% figure(2);
% saveas(gcf,'F:\Git\Anew\figures\bt1D_10HNAP_err1.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('done done done');
% figure(2);
% plot(alpm_ben,'b')
% hold on 
% plot(alpm_tor,'r')
% figure(3);
% plot(Smaxben,'b')
% hold on 
% plot(Smaxtor,'r')
