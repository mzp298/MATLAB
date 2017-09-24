% Program to get the Gauss-Legendre Quadrature results (Vectorized)
clear;clc;close all;
format long e
load('gaussian.mat');

E=215e9;               %Young's modulus
nu=0.38;                 %poisson's ratio
k=1e9;                  %hardening parameter
b=5;                      %weakening scales distribution exponent (between 1 and 2)
fb=1.1;
a=0.5;
delta_alp=1e-4;
y=418e6;            %macroscopic yield stress
stressben=326e6;            %cyclic stressben
stressbentensor= [stressben 0 0;0 0 0;0 0 0];
stepnumber=100;        %devide one cycle in 200 parts
W0=3e4;             %dissipated energy to failure per unit volume
cycles=2;
lamratio=1e-3;
lamplus=0;
lamminus=lamratio.*lamplus;
% m=0;                   % mean stress
m=50e6;                   % mean stress
sigm=m;

%---------------------2 Numerical method with alp(Smin)-----------------------------
D=1e-16;
n=1;
tensor = [m+stressben*sind(n*360/stepnumber) 0 0 ;...
    0 0 0 ;...
    0 0 0 ];
run('Damiter1.m')
g=1; %first reference index when creating
D=D+D^alp*W/W0;
while D<1
    tensor = [m+stressben*sind(n*360/stepnumber) 0 0 ;...
        0 0 0 ;...
        0 0 0 ];
    hydro(n)=1/3*trace(tensor);
    dev1=tensor-hydro(n)*eye(3)-scentre;
    dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
    dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
    dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
    tensor = [m+stressben*sind((n+1)*360/stepnumber) 0 0 ;...
        0 0 0 ;...
        0 0 0 ];
    run('Damiter2.m')
    figure(1);
    hold on;
    W_change_alp=plot ((n+1),W,'LineStyle', 'none','LineWidth', 1, 'Marker', 's', 'MarkerSize',8, ...
        'MarkerEdgeColor',  [255 193 37]/255, 'MarkerFaceColor' , [255 193 37]/255);
    sig_scale=plot (n+1,tensor(1,1)/4.2e5,'LineStyle', 'none','LineWidth', 1,'Marker', 's', 'MarkerSize', 11, ...
        'MarkerEdgeColor','none', 'MarkerFaceColor',[238 18 137]/255);

    D(n+1)=D(n)+D(n)^alp(n+1)*W/W0;
    n=n+1;
 end
hydroplus=1/sqrt(2)*(max(hydro)-mean(hydro))+mean(hydro);
hydrominus=1/sqrt(2)*(min(hydro)-mean(hydro))+mean(hydro);
disp(['num_original= ' num2str(n) ' .']);

% %---------------------2 Numerical opt-----------------------------
% n=1;
% tensor = [m+stressben*sind(n*360/stepnumber) 0 0 ;...
%     0 0 0 ;...
%     0 0 0 ];
% run('Damiter1.m')
% g=1; %first reference index when creating
% while n<cycles*stepnumber
%     tensor = [m+stressben*sind(n*360/stepnumber) 0 0 ;...
%         0 0 0 ;...
%         0 0 0 ];
%     hydro(n)=1/3*trace(tensor);
%     dev1=tensor-hydro(n)*eye(3)-scentre;
%     dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
%     dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
%     dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
%     tensor = [m+stressben*sind((n+1)*360/stepnumber) 0 0 ;...
%         0 0 0 ;...
%         0 0 0 ];
%     run('Damiter2.m')
%     n=n+1;
% end
%     e=1; %first D index
%     j=1; %first reference alp, W, n index when iterate(after adaptation)
%     D= 1e-16;
%     while D<1 %-----------the optimal time steps can be iterated with scalar
%         D=D+D^alp_ref(j)*W_ref(j)/W0;
%         j=j+1;
%         if j>=length(alp_ref)
%             j=1;
%         end
%         e=e+1;
%     end
% % plot(n_ref(1:100),W_ref(1:100),'b^')
% disp(['num_opt= ' num2str(e) ' .']);
%%
% 
%---------------------3 analytical_beforeintgration-----------------------------
alp_cyc_v2= mean(alp);
yield_plus=y-lamplus*hydroplus;     %mean positive
yield_minus=y-lamminus*hydrominus;     %mean negative
D_cyc_v2(1)=1e-16;
n_cyc_v2=1;
Smax_cyc_v2=max(Smax);
W_cyc_v2=2*(E-k)*(1+nu)*(b-1)/(E*(E+k*nu)*b*(b+1))*(Smax_cyc_v2.^(b+1)*yield_plus.^(1-b)+Smax_cyc_v2.^(b+1)*yield_minus.^(1-b)) ;
while D_cyc_v2(n_cyc_v2)<1
    D_cyc_v2(n_cyc_v2+1) = D_cyc_v2(n_cyc_v2)+ D_cyc_v2(n_cyc_v2)^alp_cyc_v2*W_cyc_v2/stepnumber/W0;
    figure(1)
    hold on;
    plot_W_cyc_v2=plot (n_cyc_v2+1,W_cyc_v2/stepnumber,'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize',10, ...
        'MarkerEdgeColor',  [102 205 0]/255, 'MarkerFaceColor' ,[102 205 0]/255);
    n_cyc_v2=n_cyc_v2+1;
end
disp(['analytical_beforeintgration= ' num2str(n_cyc_v2) ' .']);
% 
% % ---------------------3 analytical-----------------------------
% NFben_ana=W0.*E.*(E+k.*nu).*b.*(b+1).*...
%         (2.*(1-[alp_cyc_v2]).*(E-k).*(1+nu).*(b-1)).^-1.*...
%         ([Smax_cyc_v2].^(b+1).*(yield_plus).^(1-b)...
%         +[Smax_cyc_v2].^(b+1).*(yield_minus).^(1-b)).^-1;
% n=stepnumber.*NFben_ana;
% disp(['analytical= ' num2str(n) ' .']);
% 
% error=(n-e)/n;
% disp(['Error (analytical-num_opt)/analytical = ' num2str(error) '.']);

% %---------------------in loop 2 scales plot settings-----------------------------
figure(1)
grid on;
grid minor;
set(gca ,'FontSize',30);
hXLabel = xlabel('Time step(100 steps in one cycle)' ,'Fontsize' ,30);
hYLabel = ylabel('Dissipated energy per timestep', 'Fontsize' ,30);
hTitle = title('Dissipated energy comparison' ,'Fontsize' ,30);
set(hTitle, 'FontSize', 30, 'FontWeight' , 'bold')
hLegend=legend([sig_scale,plot_W_cyc_v2,W_change_alp],'Scaled \Sigma=\Sigma_{bending}/4.2e5','W with analytical method(\alpha_m)',...
    'W with numerical method(varying \alpha)','Location','southeast');
set(hLegend, 'FontSize', 25)
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
set(gcf, 'PaperPosition', [0 0 1280 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])

figure(2)
pD_cyc=plot (D_cyc_v2,'LineStyle', 'none','LineWidth', 2, 'Marker', 'o', 'MarkerSize',10, ...
    'MarkerEdgeColor',  [102 205 0]/255, 'MarkerFaceColor' ,'none');
hold on;
pD_change_alp=plot (D,'LineStyle', 'none','LineWidth', 2, 'Marker', 'v', 'MarkerSize',10, ...
    'MarkerEdgeColor',  'r', 'MarkerFaceColor' , 'none');
figure(2);
set(gca ,'FontSize',30);
hXLabel = xlabel('Time step' ,'Fontsize' ,30);
hYLabel = ylabel('Damage', 'Fontsize' ,30);
hTitle = title('Damage accumulation with time' ,'Fontsize' ,30);
set(hTitle, 'FontSize', 30, 'FontWeight' , 'bold')
hLegend=legend([pD_cyc,pD_change_alp],...
    'D with analytical method($\delta D=D^{\alpha_m}\frac{W_{cyc}}{W_0}\delta N$)',...
    'D with numerical method($\delta D=D^\alpha\frac{\dot{W}}{W_0}\delta t$)','Location','best');
set(hLegend,'Interpreter','latex');
set(hLegend, 'FontSize', 30)
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
set(gcf, 'PaperPosition', [0 0 1280 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])


% figure(1)
% saveas(gcf,'F:\Git\Anew\figures\W_3methods.png');
% % saveas(gcf,'F:\Git\Anew\figures\W_3methods_enlarge.png');
% figure(2)
% saveas(gcf,'F:\Git\Anew\figures\D_3methods2_100steps.png');

figure(1)
saveas(gcf,'F:\Git\Anew\figures\W3methods_bigbeta.png');
figure(2)
saveas(gcf,'F:\Git\Anew\figures\damsin_bigbeta.png');



sp=actxserver('SAPI.SpVoice');
sp.Speak('stressben done');
%mail2me('job finished',['Elapsed time is ' num2str(toc) ' seconds. Real test time is ' testtime ' seconds. Number of points to failure is ' NF ' points.']);
