clear;clc

alpha=0.8;

gam1=0.1;
gam2=0.5;
gam3=1.5;
ratio=0:0.1:1;

D1=(1-(1-ratio.^(1/(1-alpha)))).^(1/(gam1+1));
D2=(1-(1-ratio.^(1/(1-alpha)))).^(1/(gam2+1));
D3=(1-(1-ratio.^(1/(1-alpha)))).^(1/(gam3+1));

hold on
pd1= plot(ratio,D1,':r','LineWidth', 10);
pd2= plot(ratio,D2,'--b','LineWidth', 10);
pd3= plot(ratio,D3,'-.g','LineWidth', 10);


% title('Damage versus fatigue cycle ratio with different \gamma with \alpha=0.8');
xlabel('N/N_{F}');
ylabel('D');
grid on;
grid minor;
hXLabel = xlabel('N/N_{F}');
hYLabel =ylabel('D');

hLegend=legend([pd1,pd2,pd3], '\gamma=0.1','\gamma=0.5','\gamma=1.5','Location','southeast');
set([hLegend, gca], 'FontSize', 40)
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
saveas(gcf,'F:\Git\PhDreport\5thesis\figures\Dratio2.png');