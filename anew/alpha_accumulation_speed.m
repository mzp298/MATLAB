clear;clc
alpha1=0.9;
alpha2=0.6;
alpha3=0.3;
alpha4=0;
alpha5=-1;
alpha6=-5;
alpha7=-80;
gam=0.5;

ratio=0:0.05:1;

D1=ratio.^(1/(1-alpha1));
D2=ratio.^(1/(1-alpha2));
D3=ratio.^(1/(1-alpha3));
D4=ratio.^(1/(1-alpha4));
% D1=(1-(1-ratio.^(1/(1-alpha1)))).^(1/(gam+1));
% D2=(1-(1-ratio.^(1/(1-alpha2)))).^(1/(gam+1));
% D3=(1-(1-ratio.^(1/(1-alpha3)))).^(1/(gam+1));
hold on;
PD1 = plot(ratio,D1,'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 10, ...
     'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'm');
PD2 =  plot(ratio,D2,'LineWidth', 2, 'Marker', 'd', 'MarkerSize', 10, ...
     'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'g');
PD3 =  plot(ratio,D3,'LineWidth', 2, 'Marker', '^', 'MarkerSize', 10, ...
     'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'r');
PD4 =  plot(ratio,D4,'LineWidth', 2, 'Marker', 'h', 'MarkerSize', 10, ...
     'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'k');
% PD5 =  plot(ratio,D5,'LineWidth', 2, 'Marker', 's', 'MarkerSize', 10, ...
%      'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'b');
% PD6 =  plot(ratio,D6,'LineWidth', 2, 'Marker', 'd', 'MarkerSize', 10, ...
%      'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'g');
% PD7 =  plot(ratio,D7,'LineWidth', 2, 'Marker', '^', 'MarkerSize', 10, ...
%      'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'r');
xlabel('T/T_{F}');
ylabel('D');
% % ---------------------plot settings-----------------------------
grid on;
grid minor;
set(gca ,'FontSize',30);
hXLabel = xlabel('T/T_F' ,'Fontsize' ,30);
hYLabel =ylabel('D', 'Fontsize' ,30);
hLegend=legend([PD1,PD2,PD3,PD4],'\alpha=0.9','\alpha=0.6',...
   '\alpha=0.3','\alpha=0','Location','Best');
set([hLegend, gca], 'FontSize', 30)
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white 
% 
 % Adjust font
set(gca, 'FontName', 'Helvetica')
set([ hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hXLabel, hYLabel], 'FontSize', 30)
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
set(gcf, 'PaperPosition', [0 0 1000 800]);
% saveas(gcf,'F:\Git\Anew\figures\alpha_accumulation_speed.png');

