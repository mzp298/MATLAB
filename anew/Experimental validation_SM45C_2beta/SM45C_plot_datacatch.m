clear all; clc; close all;
NFben=[17520.33619
33991.10739
52427.26419
91077.18396
156881.5804
222261.1047
446114.6594
822487.2063
1279413.916
1453321.25
2440359.869
3428114.749
6880791.155
6213809.484
9342857.067
7240667.386]' ;

stressben=1e6.*[632.13256
590.05764
552.01729
529.5389
506.48415
489.76945
466.7147
463.83285
459.2219
463.83285
454.03458
455.18732
450
437.31988
441.93084
424.0634]';%to get Smaxben

stresstor=1e6.*[404.11255
394.87734
375.25253
363.13131
354.4733
345.8153
338.31169
331.38528
329.07648
322.15007
];%to get Smaxtor
NFtor=[27957.25096
47749.27586
76193.65712
100000
162305.0764
182806.9151
296704.9032
575635.5425
822487.2063
2203806.358
];

m=196e6;%mean tension
NFbenm=[43043.40644
55522.84924
74724.73326
75361.55953
110407.2374
146089.5545
194951.3693
212217.5648
297989.8518
440285.7643
678727.3216
597603.0102
] ;
stressbenm=1e6.*[541.48629
511.47186
514.35786
493.57864
490.69264
471.06782
455.48341
452.0202
430.08658
409.30736
407.57576
386.79654
];

figure(1)
grid on;
grid minor;
axis([1e4 1e7 250*1e6 650*1e6]);
benplot=semilogx(NFben,stressben,'ks','MarkerSize',15,'LineWidth', 3,'MarkerFaceColor','none');
hold on
torplot=semilogx(NFtor,stresstor,'k^','MarkerSize',15,'LineWidth', 3,'MarkerFaceColor','none');
benmplot=semilogx(NFbenm,stressbenm,'ks','MarkerSize',15,'LineWidth', 3,'MarkerFaceColor','k');
set(gca ,'FontSize',30);
hXLabel = xlabel('Cycles to failure of SM45C','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('Stress amplitude (Pa)','Fontsize',30, 'FontWeight' , 'bold');
hLegend=legend([benplot,benmplot,torplot],...
    'Bending experiments',['Bending experiments with ',sprintf('\n'),'196 MPa mean stress'],...
    'Torsion experiments','location','best');
set(hLegend, 'FontSize', 28);
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
set(gcf, 'PaperPosition', [0 0 1080 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])
figure(1);
saveas(gcf,'F:\Git\Doctor_thesis_Zepeng\figures\SM45C_SN.png');