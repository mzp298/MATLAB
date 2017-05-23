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
m=196e6;%mean tension
%---------------a b correspond to high and low bounds of cyclic load, 1 2
%correspond to  different parameters comparison------------------
NF=[53E3 59.2E3 70.1E3 86.3E3 89.9E3 92.1E3 102E3 135E3 351E3 394E3] ;
stressben=[441E6 286E6 464E6 473E6 173E6 403E6 437E6 167E6 357E6 182E6];%to get Smaxben
stresstor=[215E6 309E6 155E6 136E6 334E6 209E6 177E6 321E6 179E6 274E6];%to get Smaxtor
n=1;
while n<length(stressben)+1
hydromax=1/3*sum(stressben(n)+m +0+0);
hydromin=1/3*sum(-stressben(n)+m +0+0);
hydromean=1/2*sqrt(hydromax^2+hydromin^2);
devmax=[stressben(n)+m -stresstor(n) 0 ;-stresstor(n) 0  0 ;0 0 0 ]-hydromax*eye(3);
devmin=[-stressben(n)+m stresstor(n) 0 ;stresstor(n) 0 0 ;0 0 0 ]-hydromin*eye(3);
devmean=1/2*(norm(devmax,'fro')+norm(devmin,'fro'));
Smaxbt_m(n)=norm(devmean,'fro');
%--------to get average smin--------
yield_min=y-lam.*hydromax;
Smaxbtmax(n)=norm(devmax,'fro');
sminmin(n)=yield_min.*Smaxbtmax(n).^-1;

yield_max=y-lam.*hydromin;
Smaxbtmin(n)=norm(devmin,'fro');
sminmax(n)=yield_max.*Smaxbtmin(n).^-1;

yield_m(n)=(yield_min+yield_max)/2;
alphamax(n)=(1-a.*(sminmax(n)-1).^-fb);
alphamin(n)=(1-a.*(sminmin(n)-1).^-fb);
alpbt_m(n)=(alphamax(n)+alphamin(n))/2;

n=n+1;
end
%control lam to make min(sminmin)>1
min(sminmin)
%--------use average smin to get NF--------
NFbt_num=(1-alpbt_m).^-1*W0*E*(E+k*nu)*b*(b+1)*((4*(E-k)*(1+nu)*(b-1)))^-1.*yield_m.^(b-1).*Smaxbt_m.^(-b-1);

figure(1);
smax_exp=semilogx(NF,Smaxbt_m,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[208 32 144]/255, 'MarkerFaceColor','none');
hold on;
smax_num=semilogx(NFbt_num,Smaxbt_m,'-r','LineWidth', 3);

set(gca ,'FontSize',30);
hXLabel = xlabel('NF_{num}','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('S_{max}','Fontsize',30, 'FontWeight' , 'bold');
hLegend=legend([smax_exp,smax_num],...
    'Bending-torsion experimental result','Bending-torsion numerical result','location','bestoutside');
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
err_bt_m = loglog(NF,NFbt_num,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[139 69 19]/255, 'MarkerFaceColor','none');
set(gca ,'FontSize',30);
hXLabel = xlabel('NF_{exp}','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('NF_{num}','Fontsize',30, 'FontWeight' , 'bold');
x=1e4:1000:1e6;
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
hLegend=legend([err_bt_m],...
    'Bending-torsion 90 degree out-of-phase test with mean stress on SM45C','location','bestoutside');
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
saveas(gcf,'F:\Git\Anew\figures\bt2D_m_SM45C_sn.png');
figure(2);
saveas(gcf,'F:\Git\Anew\figures\bt2D_m_SM45C_err1.png');

figure(3);
plot(alpbt_m,'b','LineWidth', 3);
sp=actxserver('SAPI.SpVoice');
sp.Speak('done done done');

