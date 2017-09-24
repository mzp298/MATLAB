% E=72e9;               %Young's modulus
% nu=0.3;                 %poisson's ratio
% k=6e8;                  %hardening parameter
% b=1.1;                      %weakening scales distribution exponent (between 1 and 2)
% y=230e6;            %macroscopic yield stress
% a=0.1;
% W0=3.22e7;             %dissipated energy to failure per unit volume
% lam=0.1;               %hydrostatic pressure sensitivity
% m=0;
% load=2.25e8;            %cyclic load
% loadtensor= [load 0 0;0 0 0;0 0 0];
% stepnumber=10;        %devide one cycle in X parts
stepnumber=[2 10 1e2 1e3 1e4 1e5 ] ;%steps in unit cycle

NF_direct =[
     6.717222067146224e+02
     1.717092975449184e+03
     1.733711135857042e+03
     1.733556280409072e+03
     1.733660257156788e+03
     1.733566592854177e+03];
NF_iter =[     
     6.760000000000000e+02
     1.716800000000000e+03
     1.731750000000000e+03
     1.731450000000000e+03
     1.731537700000000e+03
     1.731443450000000e+03];
NF_num_fix =[
     6.755000000000000e+02
     1.717000000000000e+03
     1.731370000000000e+03
     1.730952000000000e+03
     1.731204500000000e+03
     1.730940250000000e+03];
NF_num_change =[ 
     6.605000000000000e+02
     1.153900000000000e+03
     1.606400000000000e+03
     1.667738000000000e+03
     1.673879300000000e+03
     1.674469860000000e+03];
figure(9)

pNF_direct=semilogx(stepnumber,NF_direct,'LineStyle', '-','Color', [160 32 240]/255,'LineWidth', 4, 'Marker', 'd', 'MarkerSize',20, ...
    'MarkerEdgeColor',  [160 32 240]/255, 'MarkerFaceColor' ,'none');
hold on
pNF_iter=semilogx(stepnumber,NF_iter,'LineStyle', '-','Color', [102 205 0]/255,'LineWidth', 4, 'Marker', 'o', 'MarkerSize',20, ...
    'MarkerEdgeColor',  [102 205 0]/255, 'MarkerFaceColor' ,'none');
pNF_num_fix=semilogx(stepnumber,NF_num_fix ,'LineStyle', '-','Color', [72 61 139]/255,'LineWidth', 4, 'Marker', 's', 'MarkerSize',20, ...
    'MarkerEdgeColor',  [72 61 139]/255, 'MarkerFaceColor' , 'none');
pNF_num_change=semilogx(stepnumber,NF_num_change,'LineStyle', '-','Color', 'r','LineWidth', 4, 'Marker', 'v', 'MarkerSize',20, ...
    'MarkerEdgeColor',  'r', 'MarkerFaceColor' , 'none');
grid on;
grid minor;
set(gca ,'FontSize',30);
hXLabel = xlabel('Time steps in unit cycle' ,'Fontsize' ,30);
hYLabel = ylabel('Number of cycles to failure', 'Fontsize' ,30);
hLegend=legend([pNF_direct,pNF_iter,pNF_num_fix,pNF_num_change],...
    '$N_F$ with integration of D($N_F=\frac{W_0}{(1-\alpha_m)W_{cyc}}$)',...
    '$N_F$ with mean $\alpha$ iteration($\delta D=D^{\alpha_m}\frac{W_{cyc}}{W_0}\delta N$)',...
    '$N_F$ with numerical method($\delta D=D^{\alpha_m}\frac{\dot{W}}{W_0}\delta t$)',...
    '$N_F$ with numerical method($\delta D=D^\alpha\frac{\dot{W}}{W_0}\delta t$)','Location','best');
set(hLegend,'Interpreter','latex');
set(hLegend, 'FontSize', 25)
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white
% Adjust font
set(gca, 'FontName', 'Helvetica')
set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
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
set(gcf, 'PaperPosition', [0 0 1280 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])
% figure(9)
% saveas(gcf,'F:\Git\Anew\figures\NF_4methods_steps_unitcycle.png');


     