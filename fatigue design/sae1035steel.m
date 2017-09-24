clear;clc;
tt=172.4;
ff=234;
a=(tt-ff/sqrt(3))/(ff/3);
b=tt;
ap=3*(tt/ff-1/2);
bp=tt;
ad=3*(tt/ff-1/2);
bd=tt;
lg=0.2861;
lgp=0.2861;
lgd=0.2861;
y=1:0.01:30;

sigmax  =b*(1/sqrt(3)+a/3-lg*((sqrt(3)*y).^-1+a/3*y.^-1)).^-1;
sigmaxp=b*(1/2+ap/3-lgp*((2*y).^-1+ap/3*y.^-1)).^-1;
sigmaxd=b*(1/2+ad/3-lgd*((2*y).^-1+ad/3*y.^-1)).^-1;

figure(1);
hold on;
fr=plot(y,sigmax, '-b','LineWidth',2,'Marker', 'o', 'MarkerSize', 12,'MarkerEdgeColor',  'b', 'MarkerFaceColor' , 'none');
frp=plot(y,sigmaxp,'-.r','LineWidth',2,'Marker', '^', 'MarkerSize', 8,'MarkerEdgeColor',  'r', 'MarkerFaceColor' , 'none');
frd=plot(y,sigmaxd, ':g','LineWidth',2,'Marker', 'v', 'MarkerSize', 4,'MarkerEdgeColor',  'g', 'MarkerFaceColor' , 'none');
plot(2,269,'d','MarkerSize',20, 'MarkerFaceColor','k');
plot(3,268,'d','MarkerSize',20, 'MarkerFaceColor','k');
plot(6.3,241,'d','MarkerSize',20, 'MarkerFaceColor','k');
plot(12.7,237,'d','MarkerSize',20, 'MarkerFaceColor','k');
plot(24,234,'d','MarkerSize',20, 'MarkerFaceColor','k');
class=plot([2 30],[ff ff],'k','LineWidth',5);
hold off;
grid on;
% grid minor;
set(gca ,'FontSize',30);
hXLabel = xlabel('radius','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('\sigma_{max}','Fontsize',30, 'FontWeight' , 'bold');
text(0.2,ff,'ref','Fontsize',30);

hTitle=title({['SAE 1035 steel'],[  'Rotating Cantilever Bending'],[ '(data from Moore 1944)']; },'Fontsize',30);
set(hTitle, 'FontSize', 30, 'FontWeight' , 'bold')
hLegend = legend([class,fr,frp,frd],'Crossland\_Classical','Crossland\_Gradient','Papadopoulos\_Gradient','DangVan\_Gradient');
set([hLegend, gca], 'FontSize', 25)
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
set(gcf, 'PaperPosition', [0 0 800 600]); %set(gcf,'PaperPosition',[left,bottom,width,height])
saveas(gcf,'1035steel.png');
saveas(gcf,'1035steel.fig');
%*****************Fitting********************
R=[2 3 6.3 12.7 24];
f=[269 268 241 237 234];
createFit1035(R, f)
% cftool
% t/(1/sqrt(3)+(t-234/sqrt(3))/234-lg*(1/(sqrt(3)*x)+(t-234/sqrt(3))/(234*x)))
% lg=0.001~0.5
% t=150~200