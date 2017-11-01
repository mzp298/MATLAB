% Program to get the Gauss-Legendre Quadrature results (Vectorized)
clear;clc;close all;

dbstop if error
format long e
load('gaussian.mat');
% vpa(x,289);
% vpa(weight,289);

E=191e9;               %Young's modulus
nu=0.38;                 %poisson's ratio
k=1e9;                  %hardening parameter
b=1.1;                      %weakening scales distribution exponent (between 1 and 2)
pb=1.1;
y=1080e6;            %macroscopic yield stress
a=0.5;
W0=7e5;             %dissipated energy to failure per unit volume
lam=0.3;               %hydrostatic pressure sensitivity
m=0;
load=8.5e8;            %cyclic load
loadtensor= [load 0 0;0 0 0;0 0 0];
stepnumber=100;        %devide one cycle in X parts


%---------------------3 Numerical method with alp(Smin)-----------------------------
D_change_alp= 1e-16; 
n=1;       
%---------------------to get the the first Sb-----------------------------
stress11=load*sind(n*360/stepnumber);
hydro=1/3*sum(stress11+0+0);

dev1=[stress11 0 0;0 0 0;0 0 0]-hydro*eye(3);
dev11=dev1(1,1); dev12=dev1(1,3); dev13=dev1(1,3);
dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
Smax=1/sqrt(2).*norm(dev1,'fro');

trial11=dev11; trial12=dev12; trial13=dev13;
trial21=dev21; trial22=dev22; trial23=dev23;
trial31=dev31; trial32=dev32; trial33=dev33;
Smaxtrial(1)=1/sqrt(2).*norm([trial11, trial12, trial13; trial21, trial22, trial23;trial31, trial32, trial33],'fro');

s= (x/2+1/2).^(1/(1-b)); 
yield(1)=y-lam*hydro;
eta=bsxfun(@minus,bsxfun(@times,Smaxtrial(1)/yield,s),1); 
eta(eta<0)=0; 

Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));

Ws=(bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield(1),s))<=0).*...
    (0)+...
    (bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield(1),s))>0).*...
    ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield,s)),yield),s)));

W= sum(Ws);
Wchange=W;
sequence=((Smax*yield(1)^-1)*(1-Smax*yield(1)^-1)^-1)^pb;
sequence(sequence<0)=0;
alp(1)=1-a*sequence;
D_change_alp(2)=D_change_alp(1)+D_change_alp(1)^alp(1)*W/W0;

while D_change_alp(n)<1
    stress11=load*sind(n*360/stepnumber);
    hydro=1/3*sum(stress11+0+0);
    dev1=[stress11 0 0;0 0 0;0 0 0]-hydro*eye(3);
    dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
    dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
    dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
    
    stress11=load*sind((n+1)*360/stepnumber);
    hydro=1/3*sum(stress11+0+0);
    devn=[stress11 0 0;0 0 0;0 0 0]-hydro*eye(3);
    dev11g=devn(1,1); dev12g=devn(1,2); dev13g=devn(1,3);
    dev21g=devn(2,1); dev22g=devn(2,2); dev23g=devn(2,3);
    dev31g=devn(3,1); dev32g=devn(3,2); dev33g=devn(3,3);
    Smax=1/sqrt(2).*norm(devn,'fro');
    
    trial11=bsxfun(@plus,Sb11,(dev11g-dev11)); trial12=bsxfun(@plus,Sb12,(dev12g-dev12));trial13=bsxfun(@plus,Sb13,(dev13g-dev13));
    trial21=bsxfun(@plus,Sb21,(dev21g-dev21)); trial22=bsxfun(@plus,Sb22,(dev22g-dev22));trial23=bsxfun(@plus,Sb23,(dev23g-dev23));
    trial31=bsxfun(@plus,Sb31,(dev31g-dev31)); trial32=bsxfun(@plus,Sb32,(dev32g-dev32));trial33=bsxfun(@plus,Sb33,(dev33g-dev33));
    trialtensor=[trial11; trial12; trial13; trial21; trial22; trial23;trial31; trial32; trial33];
    Smaxtrial=1/sqrt(2).*sqrt(sum(trialtensor.^2));
    yield(n+1)=y-lam*hydro;
    eta=bsxfun(@minus,bsxfun(@times,Smaxtrial/yield(n+1),s),1); 
    eta(eta<0)=0;
    
    Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
    Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
    Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
    
    Ws=(bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield(n+1),s))<=0).*...
        (0)+...
        (bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield(n+1),s))>0).*...
        ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*...
        bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield(n+1),s)),yield(n+1)),s)));
    
    W= sum(Ws);
    Wchange(n+1)=Wchange(n)+W;
    sequence=((Smax*yield(n+1)^-1)*(1-Smax*yield(n+1)^-1)^-1)^pb;
    sequence(sequence<0)=0;
    alp(n+1)=1-a*sequence;
    D_change_alp(n+1)=D_change_alp(n)+D_change_alp(n)^alp(n+1)*W/W0;
    n=n+1;
end
NF_num_change=n/stepnumber
disp(['Mean alpha to give fix alpha value is ' num2str(mean(alp)) ' .']);
figure(1)
hold on;
W_change_alp=plot (1:n,Wchange(1:n),'LineStyle', 'none','LineWidth', 2, 'Marker', 'v', 'MarkerSize',10, ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor' , 'none');
figure(3)
hold on;
pD_change_alp=plot (D_change_alp,'LineStyle', 'none','LineWidth', 2, 'Marker', 'v', 'MarkerSize',10, ...
    'MarkerEdgeColor',  'r', 'MarkerFaceColor' , 'none');

%% 
alp0=mean(alp);
hydro=1/3*trace(loadtensor);
yield_m=y-lam*hydro;
%---------------------2 Direct without iteration-----------------------------
NF_direct=(1-alp0).^-1*W0*E*(E+k*nu)*b*(b+1)*((4*(E-k)*(1+nu)*(b-1)))^-1.*yield_m.^(b-1).*(1/sqrt(2).*norm(loadtensor-hydro*eye(3),'fro')).^(-b-1)*(1-1e-16^(1-alp0));
%---------------------2 Cyclic load calculation-----------------------------
D_cyc(1)=1e-16;
n=1;
Wcyc=4*(E-k)*(1+nu)*(b-1)/(E*(E+k*nu)*b*(b+1))*(1/sqrt(2).*norm(loadtensor-hydro*eye(3),'fro')).^(b+1)*yield_m.^(1-b) ;
Wcyc_step=Wcyc/stepnumber;
Wcycsum(1)=Wcyc_step;
while D_cyc(n)< 1
    D_cyc(n+1) = D_cyc(n)+ D_cyc(n)^alp0*Wcyc_step/W0;
    Wcycsum(n+1)=Wcycsum(n)+Wcyc_step;
    n=n+1;
end
NF_iter=n/stepnumber

figure(1)
hold on;
W_cyc=plot (1:n,Wcycsum(1:n),'LineStyle', 'none','LineWidth', 2, 'Marker', 'o', 'MarkerSize',10, ...
    'MarkerEdgeColor',  [102 205 0]/255, 'MarkerFaceColor' ,'none');
figure(3)
hold on;
pD_cyc=plot (D_cyc,'LineStyle', 'none','LineWidth', 2, 'Marker', 'o', 'MarkerSize',10, ...
    'MarkerEdgeColor',  [102 205 0]/255, 'MarkerFaceColor' ,'none');



%% 

%---------------------3 Numerical method with fixed alp-----------------------------
clear Smaxtrial;
D_fix_alp= 1e-16; 
n=1;       
%---------------------to get the the first Sb-----------------------------
stress11=load*sind(n*360/stepnumber); 
hydro=1/3*sum(stress11+0+0);

dev1=[stress11 0 0;0 0 0;0 0 0]-hydro*eye(3);
dev11=dev1(1,1); dev12=dev1(1,3); dev13=dev1(1,3);
dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
Smax=1/sqrt(2).*norm(dev1,'fro');

trial11=dev11; trial12=dev12; trial13=dev13;
trial21=dev21; trial22=dev22; trial23=dev23;
trial31=dev31; trial32=dev32; trial33=dev33;
Smaxtrial(1)=1/sqrt(2).*norm([trial11, trial12, trial13; trial21, trial22, trial23;trial31, trial32, trial33],'fro');

s= (x/2+1/2).^(1/(1-b)); 
yield=y-lam*hydro;
eta=bsxfun(@minus,bsxfun(@times,Smaxtrial(1)/yield,s),1); 
eta(eta<0)=0;

Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));

Ws=(bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield,s))<=0).*...
    (0)+...
    (bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield,s))>0).*...
    ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield,s)),yield),s)));

W= sum(Ws);
Wfix=W; 
D_fix_alp(2)=D_fix_alp(1)+D_fix_alp(1)^alp0*W/W0;

while D_fix_alp(n)<1
    stress11=load*sind(n*360/stepnumber);
    hydro=1/3*sum(stress11+0+0);
    dev1=[stress11 0 0;0 0 0;0 0 0]-hydro*eye(3);
    dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
    dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
    dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
    
    stress11=load*sind((n+1)*360/stepnumber);
    hydro=1/3*sum(stress11+0+0);
    devn=[stress11 0 0;0 0 0;0 0 0]-hydro*eye(3);
    dev11g=devn(1,1); dev12g=devn(1,2); dev13g=devn(1,3);
    dev21g=devn(2,1); dev22g=devn(2,2); dev23g=devn(2,3);
    dev31g=devn(3,1); dev32g=devn(3,2); dev33g=devn(3,3);
    Smax=1/sqrt(2).*norm(devn,'fro');
    
    trial11=bsxfun(@plus,Sb11,(dev11g-dev11)); trial12=bsxfun(@plus,Sb12,(dev12g-dev12));trial13=bsxfun(@plus,Sb13,(dev13g-dev13));
    trial21=bsxfun(@plus,Sb21,(dev21g-dev21)); trial22=bsxfun(@plus,Sb22,(dev22g-dev22));trial23=bsxfun(@plus,Sb23,(dev23g-dev23));
    trial31=bsxfun(@plus,Sb31,(dev31g-dev31)); trial32=bsxfun(@plus,Sb32,(dev32g-dev32));trial33=bsxfun(@plus,Sb33,(dev33g-dev33));
    trialtensor=[trial11; trial12; trial13; trial21; trial22; trial23;trial31; trial32; trial33];
    Smaxtrial=1/sqrt(2).*sqrt(sum(trialtensor.^2));
    yield=y-lam*hydro;
    eta=bsxfun(@minus,bsxfun(@times,Smaxtrial/yield,s),1); 
    eta(eta<0)=0;
    
    Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
    Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
    Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));

   
    Ws=(bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield,s))<=0).*...
        (0)+...
        (bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield,s))>0).*...
        ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield,s)),yield),s)));
    
    W= sum(Ws);
    Wfix(n+1)=Wfix(n)+W;
        
    D_fix_alp(n+1)=D_fix_alp(n)+D_fix_alp(n)^alp0*W/W0;
    
    n=n+1;
end
NF_num_fix=n/stepnumber

%% 

figure(1)
hold on;
W_fix_alp=plot (1:n,Wfix(1:n),'LineStyle', 'none','LineWidth', 2, 'Marker', 's', 'MarkerSize',10, ...
    'MarkerEdgeColor',  [72 61 139]/255, 'MarkerFaceColor' , 'none');
figure(3)
hold on;
pD_fix_alp=plot (D_fix_alp,'LineStyle', 'none','LineWidth', 2, 'Marker', 's', 'MarkerSize',10, ...
    'MarkerEdgeColor',  [72 61 139]/255, 'MarkerFaceColor' , 'none');

%---------------------in loop 2 scales plot settings-----------------------------
figure(1)
grid on;
grid minor;
set(gca ,'FontSize',30);
hXLabel = xlabel('Time step' ,'Fontsize' ,30);
hYLabel = ylabel('Total dissipated energy', 'Fontsize' ,30);
hTitle = title('Dissipated energy accumulation speed' ,'Fontsize' ,30);
set(hTitle, 'FontSize', 30, 'FontWeight' , 'bold')
hLegend=legend([W_cyc,W_fix_alp,W_change_alp],...
    'W with analytical method($\delta D=D^{\alpha_m}\frac{W_{cyc}}{W_0}\delta N$)',...
    'W with numerical method($\delta D=D^{\alpha_m}\frac{\dot{W}}{W_0}\delta t$)',...
    'W with numerical method($\delta D=D^\alpha\frac{\dot{W}}{W_0}\delta t$)','Location','best');
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


figure(3);
set(gca ,'FontSize',30);
hXLabel = xlabel('Time step' ,'Fontsize' ,30);
hYLabel = ylabel('Damage', 'Fontsize' ,30);
hTitle = title('Damage accumulation with time' ,'Fontsize' ,30);
set(hTitle, 'FontSize', 30, 'FontWeight' , 'bold')
hLegend=legend([pD_cyc,pD_fix_alp,pD_change_alp],...
    'D with analytical method($\delta D=D^{\alpha_m}\frac{W_{cyc}}{W_0}\delta N$)',...
    'D with numerical method($\delta D=D^{\alpha_m}\frac{\dot{W}}{W_0}\delta t$)',...
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


error_W=(Wcycsum(1:round(n*0.90))-Wchange(1:round(n*0.90))).*(Wcycsum(1:round(n*0.90)).^-1);
error_D=(D_cyc(1:round(n*0.90))-D_change_alp(1:round(n*0.90))).*(D_cyc(1:round(n*0.90)).^-1);
figure(5);
W_diff_fix_change_alp=plot (1:round(n*0.90),error_W,...
    'LineStyle', 'none','LineWidth', 1, 'Marker', 's', 'MarkerSize',8, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor' , 'k');
grid on;
grid minor;
ylim([min(error_W) max(error_W)]);
set(gca ,'FontSize',30);
hXLabel = xlabel('Time step' ,'Fontsize' ,30);
hYLabel = ylabel('Relative difference ', 'Fontsize' ,30);
hLegend=legend(['Relative difference of accumulation between ',sprintf('\n'),'$\frac{W_{cyc}}{stepnumber}$ and $\dot{W}$'],...
    'Location','best');
set(hLegend,'Interpreter','latex');
set(hLegend, 'FontSize', 30)
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); 
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


figure(7);
D_diff_fix_change_alp=plot (1:round(n*0.90),error_D,...
    'LineStyle', 'none','LineWidth', 1, 'Marker', 's', 'MarkerSize',8, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor' , 'k');
grid on;
grid minor;
ylim([min(error_D) max(error_D)]);
set(gca ,'FontSize',30);
hXLabel = xlabel('Time step' ,'Fontsize' ,30);
hYLabel = ylabel('Relative difference ', 'Fontsize' ,30);
hLegend=legend(['Relative difference of accumulation between ',sprintf('\n'),'$\delta D=D^{\alpha_m}\frac{W_{cyc}}{W_0}\delta N$ and $\delta D=D^\alpha\frac{\dot{W}}{W_0}\delta t$'],...
    'Location','best');
set(hLegend,'Interpreter','latex');
set(hLegend, 'FontSize', 30)
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); 
set(gca, 'FontName', 'Helvetica')
set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1280 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])

% figure(1)
% saveas(gcf,'F:\Git\Doctor_thesis_Zepeng\figures\W_3methods2_100steps.png');
% figure(3)
% saveas(gcf,'F:\Git\Doctor_thesis_Zepeng\figures\D_3methods2_100steps.png');
% figure(5)
% saveas(gcf,'F:\Git\Doctor_thesis_Zepeng\figures\W_3methods_diff_100steps.png');
% figure(7)
% saveas(gcf,'F:\Git\Doctor_thesis_Zepeng\figures\D_3methods_diff_100steps.png');
% 
% figure(1) %after partial enlargement
% saveas(gcf,'F:\Git\Doctor_thesis_Zepeng\figures\W_3methods_100steps_enlarge.png');
% figure(3) %after partial enlargement
% saveas(gcf,'F:\Git\Doctor_thesis_Zepeng\figures\D_3methods_100steps_enlarge.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('energy damage');
%mail2me('job finished',['Elapsed time is ' num2str(toc) ' seconds. Real test time is ' testtime ' seconds. Number of points to failure is ' NF ' points.']);
