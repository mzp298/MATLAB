clear;clc;
x=1e4:1e5:1e8;
y0=x;
py0=loglog(x,y0,'k','LineWidth',3);
hold on;

cyclic_exp=[9.98000E+04
4.14000E+05];
experiments_a=[
2.59000E+06
4.22000E+06
3.90000E+06
2.44000E+06
5.23000E+06
5.32000E+06
1.39000E+07
1.04000E+07
1.14000E+07
1.10000E+07
];
experiments_b=[3.18421E+06
1.06842E+07
1.19474E+07
1.36316E+07
1.45526E+07
1.61053E+07
6.65789E+06
5.15789E+06
4.68421E+06
3.23684E+06
];
cyclic_num=[1.00958E+05
3.98367E+05];
numeric_a=[
2.65672E+06
3.06662E+06
3.07835E+06
3.15143E+06
3.03257E+06
2.94109E+06
1.11144E+07
1.10223E+07
1.09542E+07
1.09782E+07
];
numeric_b=[3.18533E+06
1.15166E+07
1.15101E+07
1.14230E+07
1.13941E+07
1.15575E+07
3.35360E+06
3.49238E+06
3.34020E+06
3.37473E+06
];
err_cyc = loglog(cyclic_exp,cyclic_num,'o','MarkerSize',15,'LineWidth', 3,'MarkerEdgeColor',[208 32 144]/255, 'MarkerFaceColor','none');
err_a = loglog(experiments_a,numeric_a,'^','MarkerSize',15,'LineWidth', 3,'MarkerEdgeColor',[238 118 0]/255, 'MarkerFaceColor','none');
err_b = loglog(experiments_b,numeric_b,'s','MarkerSize',15,'LineWidth', 3,'MarkerEdgeColor',[106 90 205]/255, 'MarkerFaceColor','none');
set(gca ,'FontSize',35);
hXLabel = xlabel('n_{exp}','Fontsize',40, 'FontWeight' , 'bold');
hYLabel = ylabel('n_{num}','Fontsize',40, 'FontWeight' , 'bold');

y1=exp(1).*x;
y2=exp(-1).*x;
py1=loglog(x,y1, '--k','LineWidth',3);
py2=loglog(x,y2, '--k','LineWidth',3);
axis equal;
axis([1e4 1e8 1e4 1e8]);

set(gca,'xtick',[1e4 1e5 1e6 1e7 1e8]); 
set(gca,'ytick',[1e4 1e5 1e6 1e7 1e8]); 
%add legend
hLegend=legend([err_cyc,err_a,err_b],'Constant amplitude test of batch a',...
    'Random amplitude test of batch a','Random amplitude test of batch b','location','bestoutside');
set([hLegend, gca], 'FontSize', 20);
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
set(gcf, 'PaperPosition', [0 0 1000 700]); %set(gcf,'PaperPosition',[left,bottom,width,height])
saveas(gcf,'F:\Git\Anew\figures\Cetim_err.png');
