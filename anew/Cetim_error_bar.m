clear;clc;close all;
x=1e4:1e5:1e8;
y0=x;
py0=loglog(x,y0,'k','LineWidth',3);
hold on;

cyclic_exp=[99892
414298];
experiments_a=[2500000
4105263
3815789
2368421
5105263
5184211
13552632
10131579
11157895
10763158];
experiments_b=[3184211
10684211
11947368
13631579
14552632
16105263
6657895
5157895
4684211
3236842];
cyclic_num=[1.42499E+05
4.90448E+05
];
numeric_a=[3.48329E+06
6.16976E+06
4.48783E+06
3.99816E+06
4.96007E+06
6.33447E+06
8.38421E+06
7.11392E+06
7.20350E+06
8.09421E+06
];
numeric_b=[4.90717E+06
7.72917E+06
8.69919E+06
8.23941E+06
7.02594E+06
6.89510E+06
5.54561E+06
5.90673E+06
4.54425E+06
4.06431E+06
];
err_cyc = loglog(cyclic_exp,cyclic_num,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[208 32 144]/255, 'MarkerFaceColor','none');
err_a = loglog(experiments_a,numeric_a,'^','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[238 118 0]/255, 'MarkerFaceColor','none');
err_b = loglog(experiments_b,numeric_b,'s','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[106 90 205]/255, 'MarkerFaceColor','none');
set(gca ,'FontSize',30);
hXLabel = xlabel('n_{exp}','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('n_{num}','Fontsize',30, 'FontWeight' , 'bold');

y1=2.*x;
y2=0.5.*x;
py1=loglog(x,y1, '--k','LineWidth',3);
py2=loglog(x,y2, '--k','LineWidth',3);
axis equal;
axis([1e4 1e8 1e4 1e8]);

set(gca,'xtick',[1e4 1e5 1e6 1e7 1e8]); 
set(gca,'ytick',[1e4 1e5 1e6 1e7 1e8]); 
%add legend
hLegend=legend([err_cyc,err_a,err_b],'Constant amplitude test of batch a',...
    'Random amplitude test of batch a','Random amplitude test of batch b','location','bestoutside');
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
set(gcf, 'PaperPosition', [0 0 1200 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])
% saveas(gcf,'F:\Git\Anew\figures\Cetim_err.png');
sp=actxserver('SAPI.SpVoice');
sp.Speak('done done done');
