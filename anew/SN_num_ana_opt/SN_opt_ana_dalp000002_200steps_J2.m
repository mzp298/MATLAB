clear;clc;close all;format long;
load('gaussian.mat');
%---------------------Verified parameters in random loading case-----------------------------
b=1.1;                    %weakening scales distribution exponent
W0=3.27e8;            %dissipated energy to failure per unit volume
%---------------------Verified parameters in constant loading case-----------------------------
y=1080e6;           %macroscopic yield stress
E=191e9;              %Young's modulus
k=1e9;                 %hardening parameter
nu=0.38;                     %poisson's ratio
lamplus=0.3;               %hydrostatic pressure sensitivity
lamminus=0.3;
fb=1.1;                %Major damage power
%---------------------Verified parameters in constant loading case-----------------------------
stepnumber=200;        %devide one cycle in 200 parts
delta_alp=2e-5;       %optimal alpha filter threshold
cycles=2;          %numerical cycles to adaptation
m=0;                   % mean stress
sigm=m;
load = [y/1.2 y/1.3 y/1.5 y/2 y/2.5 y/3 y/3.5 y/4 y/4.5 y/5 y/5.5 y/6 y/7];
a=0.01

%---------------------1 Numerical method with optimal time steps-----------------------------
for  i=1:length(load) %experimental points index
    n=1;
    tensor = [load(i)*sind(n*360/stepnumber)+m 0 0 ;...
        0 0 0 ;...
        0 0 0 ];
    run('Damiter1.m')
    
    while n<cycles*stepnumber
        tensor = [load(i)*sind(n*360/stepnumber)+m, 0, 0 ;...
            0, 0, 0 ;...
            0, 0, 0 ; ];
        hydro=1/3*trace(tensor);
        dev1=tensor-hydro*eye(3);
        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
        dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
        dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
        tensor = [load(i)*sind((n+1)*360/stepnumber)+m, 0, 0 ;...
            0, 0,  0 ;...
            0, 0,  0 ; ];
        run('Damiter2.m')
        n=n+1;
    end
    alp_m(i)=mean(alp);
    
    %----------scalar iteration
    j=ceil(length(alp_ref)/cycles); %first reference alp, W, n index when iterate(after adaptation)
    D=1e-16;
    n=1; %first D index
    while D<1 %-----------the optimal time steps can be iterated with scalar
        D=D+D^alp_ref(j)*W_ref(j)/W0; % explicit(require D \neq 0, normally D=1e-16)
        j=j+1;
        if j+1>=length(alp_ref)
            j=ceil(length(alp_ref)/cycles);
        end
        n=n+1;
    end
    NF_num_opt(i)=(n-1)/(length(n_ref));
end
reduction_factor=length(alp_ref)/stepnumber  % to see how much optimal timesteps comparing to the original


%---------------------3 Analytical calculation after integrated of D with mean alp from numerical original-----------------------------
m=0;
D0=1e-16;
for i=1:length(load)
    tensor= [load(i) 0 0;0 0 0;0 0 0];
    hydro=1/3*trace(tensor);
    yield=y-lamplus*hydro; %macro yield strength considering mean stress effect
    dev1=tensor-hydro*eye(3);
    Smax1(i)=1/sqrt(2).*norm(dev1,'fro');
    Wcyc1(i)=4*(E-k)*(1+nu)*(b-1)/(E*(E+k*nu)*b*(b+1)).*Smax1(i).^(b+1).*yield.^(1-b) ;
    NF_direct(i)=W0.*(Wcyc1(i).*(1-alp_m(i))).^-1*(1-D0^(1-alp_m(i)));
end
meanalp=alp_m'
optimalNF=NF_num_opt'
analyticalNF=NF_direct'
error=((NF_num_opt(1:length(load))-NF_direct(1:length(load))).*NF_num_opt(1:length(load)).^-1)'
save('SN_opt_ana_dalp000002_200steps_J2.mat','y','load','NF_num_opt','NF_direct','error');
%% 




%-------------plot-----------------------------
clear;clc;close all;
load('SN_opt_ana_dalp000002_200steps_J2.mat');
figure(1)
numerical_opt=semilogx(NF_num_opt,load, 'o','MarkerSize',35, 'MarkerFaceColor','none', 'LineWidth', 5,'MarkerEdgeColor','g');
hold on;
analytical_integration=semilogx(NF_direct,load,  '^','MarkerSize',35, 'MarkerFaceColor','none', 'LineWidth', 5,'MarkerEdgeColor','r');
 hLegend=legend([...
      numerical_opt,...
  analytical_integration,...
    ],...
    '1 Numerical results with optimal time steps ($\delta D=D^{\alpha_{ref}}\frac{W_{ref}}{W_0}\delta t_{ref}$)',...
    '2 Analytical results with $\alpha_m$ from numerical calculation ($N_F=\frac{W_0}{( 1-\alpha_m)W_{cyc}}$)',...
          'Location','best');
%             '1 Numerical results with original time steps ($\delta D=D^\alpha\frac{\dot{W}}{W_0}\delta t$)',...
%         '3 Analytical results with $\alpha_m$ from numerical calculation ($N_F=\frac{W_0}{( 1-\alpha_m)W_{cyc}}$)',...set(hLegend,'Interpreter','latex');
set([hLegend, gca], 'FontSize', 35)
set(hLegend,'Interpreter','latex');
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white
grid on;
grid minor;
ylim([0 y]);
set(gca,'XTick',[1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9],'FontSize',40 );
hXLabel = xlabel('N_F');
hYLabel =ylabel('Stress');
set(gca, 'FontName', 'Helvetica')
set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hXLabel, hYLabel], 'FontSize',40)
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1600 1000]); %set(gcf,'PaperPosition',[left,bottom,width,height])

figure(2)
errorbar_plot=semilogx(NF_direct,error,'k','Marker','*','MarkerSize',35,'MarkerFaceColor','none','LineStyle','none', 'LineWidth', 5);
hLegend=legend([errorbar_plot],...
    ['Relative error between analytical and numerical results',sprintf('\n'),'($\lambda_{+-}=0.3$,timesteps=200)'],...
    'Location','northeast');
set(hLegend,'Interpreter','latex');
set([hLegend, gca], 'FontSize', 35)
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white
grid on;
grid minor;
ylim([-0.04 0.04]);
set(gca,'XTick',[1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9],'FontSize',40 );
hXLabel = xlabel('N_F');
hYLabel =ylabel('Relative error');
set(gca, 'FontName', 'Helvetica')
set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hXLabel, hYLabel], 'FontSize',40)
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1600 1000]); %set(gcf,'PaperPosition',[left,bottom,width,height])

figure(1)
saveas(gcf,'F:\Git\Anew\figures\SN_opt_ana_200_delta_alp=0.00002.png');
figure(2)
saveas(gcf,'F:\Git\Anew\figures\SN_opt_ana_200_delta_alp=0.00002_err.png');
