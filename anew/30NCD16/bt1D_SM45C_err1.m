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
W0=45e9;
lam=0.1;
b=1.2;
fb=b;
m=0;%mean tension

%---------------a b correspond to high and low bounds of cyclic load, 1 2
%correspond to  different parameters comparison------------------
NFben=[1.9E4 3.6E4 5.8E4 9E4 1.7E5 2.2E5 4.4E5 8E5 1.4E6 1.5E6 2.4E6 3.3E6 6E6 6.6E6 7E6 9E6] ;
stressben=[6.2E8 5.9E8 5.52E8 5.35E8 5.05E8 4.9E8 4.7E8 4.65E8 4.62E8 4.65E8 4.6E8 4.6E8 4.4E8 4.58E8 4.3E8 4.48E8];%to get Smaxben
n=1;
while n<length(stressben)+1
hydromax=1/3*sum(stressben(n)+m +0+0);
hydromin=1/3*sum(-stressben(n)+m +0+0);
hydromean=1/2*sqrt(hydromax^2+hydromin^2);
devmax=[stressben(n)+m 0 0 ;0 0  0 ;0 0 0 ]-hydromax*eye(3);
devmin=[-stressben(n)+m 0 0 ;0 0 0 ;0 0 0 ]-hydromin*eye(3);
devmean=1/2*(norm(devmax,'fro')+norm(devmin,'fro'));
Smaxben_m(n)=norm(devmean,'fro');
%--------to get average smin--------
yield_min=y-lam.*hydromax;
Smaxbenmax(n)=norm(devmax,'fro');
sminmin(n)=yield_min.*Smaxbenmax(n).^-1;

yield_max=y-lam.*hydromin;
Smaxbenmin(n)=norm(devmin,'fro');
sminmax(n)=yield_max.*Smaxbenmin(n).^-1;

yield_m(n)=(yield_min+yield_max)/2;
alphamax(n)=(1-a.*(sminmax(n)-1).^-fb);
alphamin(n)=(1-a.*(sminmin(n)-1).^-fb);
alpben_m(n)=(alphamax(n)+alphamin(n))/2;

n=n+1;
end
%control lam to make min(sminmin)>1
min(sminmin)
%--------use average smin to get NF--------
NFben_num=(1-alpben_m).^-1*W0*E*(E+k*nu)*b*(b+1)*((4*(E-k)*(1+nu)*(b-1)))^-1.*yield_m.^(b-1).*Smaxben_m.^(-b-1);

%-------------------------------------------------torsion--------------------------------
stresstor=[4.05E8 3.99E8 3.80E8 3.65E8 3.55E8 3.49E8 3.40E8 3.35E8 3.33E8 3.25E8];%to get Smaxtor
NFtor=[2.9E4 5E4 7.8E4 1E5 1.8E5 1.9E5 3E5 5.8E5 8.2E5 2.25E6];
n=1;
while n<length(stresstor)+1
hydromax=1/3*sum(m +0+0);
hydromin=1/3*sum(m +0+0);
hydromean=1/2*sqrt(hydromax^2+hydromin^2);
devmax=[m stresstor(n) 0 ;stresstor(n) 0  0 ;0 0 0 ]-hydromax*eye(3);
devmin=[m stresstor(n) 0 ;stresstor(n) 0 0 ;0 0 0 ]-hydromin*eye(3);
devmean=1/2*(norm(devmax,'fro')+norm(devmin,'fro'));
Smaxtor_m(n)=norm(devmean,'fro');
%--------to get average smin--------
yield_min=y-lam.*hydromax;
Smaxtormax(n)=norm(devmax,'fro');
sminmin(n)=yield_min.*Smaxtormax(n).^-1;

yield_max=y-lam.*hydromin;
Smaxtormin(n)=norm(devmin,'fro');
sminmax(n)=yield_max.*Smaxtormin(n).^-1;

yieldtor_m(n)=(yield_min+yield_max)/2;
alphamax(n)=(1-a.*(sminmax(n)-1).^-fb);
alphamin(n)=(1-a.*(sminmin(n)-1).^-fb);
alptor_m(n)=(alphamax(n)+alphamin(n))/2;

n=n+1;
end
%control lam to make min(sminmin)>1
min(sminmin)
%--------use average smin to get NF--------
NFtor_num=(1-alptor_m).^-1*W0*E*(E+k*nu)*b*(b+1)*((4*(E-k)*(1+nu)*(b-1)))^-1.*yieldtor_m.^(b-1).*Smaxtor_m.^(-b-1);


figure(1);
ben_exp=semilogx(NFben,Smaxben_m,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[208 32 144]/255, 'MarkerFaceColor','none');
hold on;
tor_exp=semilogx(NFtor,Smaxtor_m,'v','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor','b', 'MarkerFaceColor','none');
ben_num=semilogx(NFben_num,Smaxben_m,'-','color',[208 32 144]/255,'LineWidth', 3);
tor_num=semilogx(NFtor_num,Smaxtor_m,'--b','LineWidth', 3);
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
set(gca ,'FontSize',30);
hXLabel = xlabel('NF_{exp}','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('NF_{num}','Fontsize',30, 'FontWeight' , 'bold');
x=1e4:1000:1e8;
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
hLegend=legend([err_ben,err_tor],...
    'Bending test on SM45(R=-1)','Torsion test on SM45(R=-1)',...
   'location','bestoutside');
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
