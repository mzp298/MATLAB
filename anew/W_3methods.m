
% Program to get the Gauss-Legendre Quadrature results (Vectorized)
clear;clc;close all;

dbstop if error
format long e

load('gaussian.mat');

E=191e9;               %Young's modulus
nu=0.38;                 %poisson's ratio
k=1e9;                  %hardening parameter
b=1.1;                      %weakening scales distribution exponent (between 1 and 2)
y=1080e6;            %macroscopic yield stress
sigu=1200e6;             %ultimite stress
ff=690e6;              %bending fatigue limit
tt=428e6;                  %torsion fatigue limit

load=6e8;            %cyclic load
loadtensor= [load 0 0;0 0 0;0 0 0];
stepnumber=200;        %devide one cycle in 200 parts
f=50;                            %frequency of load
steptime=1/f/stepnumber;
delta=(b+1)/(b-1);
a=0.5;
alp0=0.78165;
W0=1.5e5;             %dissipated energy to failure per unit volume
n0=3;                   %number of initial local defects
lam=0.3;               %hydrostatic pressure sensitivity
hydro=1/3*sum(diag(loadtensor));
yield=y-lam*hydro;                   %our mean stress

%---------------------2 Cyclic load calculation-----------------------------
% alp0=mean(alp)
Dcyc(1)=1e-16;
n=1;
Wcyc=4*(E-k)*(1+nu)*(b-1)/(E*(E+k*nu)*b*(b+1))*norm(loadtensor-hydro*eye(3),'fro').^(b+1)*yield.^(1-b) ;
while Dcyc(n)< 1
    Dcyc(n+1) = Dcyc(n)+ Dcyc(n)^alp0*Wcyc/stepnumber/W0;
    hold on;
    W_cyc=plot (n+1,Wcyc/stepnumber,'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize',10, ...
        'MarkerEdgeColor',  [	102 205 0]/255, 'MarkerFaceColor' ,[	102 205 0]/255);
    n=n+1;
end
n
%---------------------3 Numerical method with fixed alp-----------------------------
D= ones(1,1e6)*1e-16; %Pre-allocate memory for vectors
n=1;       %initial recording point
%---------------------to get the the first Sb-----------------------------
stress11=load*sind(n*360/stepnumber);
hydro=1/3*sum(stress11+0+0);

dev1=[stress11 0 0;0 0 0;0 0 0]-hydro*eye(3);
dev11=dev1(1,1); dev12=dev1(1,3); dev13=dev1(1,3);
dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
Smax=norm(dev1,'fro');

trial11=dev11; trial12=dev12; trial13=dev13;
trial21=dev21; trial22=dev22; trial23=dev23;
trial31=dev31; trial32=dev32; trial33=dev33;
normtrial(1)=norm([trial11, trial12, trial13; trial21, trial22, trial23;trial31, trial32, trial33],'fro');

s= (x/2+1/2).^(1/(1-b)); %1*25
yield=y-lam*hydro;
eta=bsxfun(@minus,bsxfun(@times,normtrial(1)/yield,s),1); %compare normtrial with yield/s
eta(eta<0)=0; %only keep normtrials which are larger than yield/s

Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
%1*25 for each Sb element
Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
normSb=sqrt(sum(Sbtensor.^2)); %sum(a) sums all the colume

% existsOnGPU(normSb)
Ws=(bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))<=0).*...
    (0)+...
    (bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))>0).*...
    ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s)),yield),s)));

W= sum(Ws);
D(2)=D(1)+D(1)^alp0*W/W0;
%%

tic;
while D(n)<1
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
    Smax=norm(devn,'fro');
    
    trial11=bsxfun(@plus,Sb11,(dev11g-dev11)); trial12=bsxfun(@plus,Sb12,(dev12g-dev12));trial13=bsxfun(@plus,Sb13,(dev13g-dev13));
    trial21=bsxfun(@plus,Sb21,(dev21g-dev21)); trial22=bsxfun(@plus,Sb22,(dev22g-dev22));trial23=bsxfun(@plus,Sb23,(dev23g-dev23));
    trial31=bsxfun(@plus,Sb31,(dev31g-dev31)); trial32=bsxfun(@plus,Sb32,(dev32g-dev32));trial33=bsxfun(@plus,Sb33,(dev33g-dev33));
    trialtensor=[trial11; trial12; trial13; trial21; trial22; trial23;trial31; trial32; trial33];
    normtrial=sqrt(sum(trialtensor.^2));
    yield=y-lam*hydro;
    eta=bsxfun(@minus,bsxfun(@times,normtrial/yield,s),1); %1*25
    eta(eta<0)=0;
    
    Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
    Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
    Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
    %1*25 for each Sb element
    Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
    normSb=sqrt(sum(Sbtensor.^2)); %sum(a) sums all the colume
    
    Ws=(bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))<=0).*...
        (0)+...
        (bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))>0).*...
        ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s)),yield),s)));
    
    W= sum(Ws);
    hold on;
    W_fix_alp=plot ((n+1),W,'LineStyle', 'none','LineWidth', 2, 'Marker', 's', 'MarkerSize',14, ...
        'MarkerEdgeColor',  [72 61 139]/255, 'MarkerFaceColor' , 'none');
    
    D(n+1)=D(n)+D(n)^alp0*W/W0;
    
    n=n+1;
end
toc;
n
t=n/stepnumber*1/f;
disp(['Time to failure is ' num2str(t) ' s.']);
Nf=n*stepnumber^-1;
disp(['Cycles to failure is ' num2str(Nf) ' cycles.']);


%---------------------3 Numerical method with alp(Smin)-----------------------------
D2= ones(1,1e6)*1e-16; %Pre-allocate memory for vectors
n=1;       %initial recording point
%---------------------to get the the first Sb-----------------------------
stress11=load*sind(n*360/stepnumber);
hydro=1/3*sum(stress11+0+0);

dev1=[stress11 0 0;0 0 0;0 0 0]-hydro*eye(3);
dev11=dev1(1,1); dev12=dev1(1,3); dev13=dev1(1,3);
dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
Smax=norm(dev1,'fro');

trial11=dev11; trial12=dev12; trial13=dev13;
trial21=dev21; trial22=dev22; trial23=dev23;
trial31=dev31; trial32=dev32; trial33=dev33;
normtrial(1)=norm([trial11, trial12, trial13; trial21, trial22, trial23;trial31, trial32, trial33],'fro');

s= (x/2+1/2).^(1/(1-b)); %1*25
yield=y-lam*hydro;
eta=bsxfun(@minus,bsxfun(@times,normtrial(1)/yield,s),1); %compare normtrial with yield/s
eta(eta<0)=0; %only keep normtrials which are larger than yield/s

Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
%1*25 for each Sb element
Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
normSb=sqrt(sum(Sbtensor.^2)); %sum(a) sums all the colume

% existsOnGPU(normSb)
Ws=(bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))<=0).*...
    (0)+...
    (bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))>0).*...
    ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s)),yield),s)));

W= sum(Ws);
sequence=((Smax*yield^-1)*(1-Smax*yield^-1)^-1)^b;
sequence(sequence<0)=0;
alp(1)=1-a*sequence;
D2(2)=D2(1)+D2(1)^alp(1)*W/W0;
%%

tic;
while D2(n)<1
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
    Smax=norm(devn,'fro');
    
    trial11=bsxfun(@plus,Sb11,(dev11g-dev11)); trial12=bsxfun(@plus,Sb12,(dev12g-dev12));trial13=bsxfun(@plus,Sb13,(dev13g-dev13));
    trial21=bsxfun(@plus,Sb21,(dev21g-dev21)); trial22=bsxfun(@plus,Sb22,(dev22g-dev22));trial23=bsxfun(@plus,Sb23,(dev23g-dev23));
    trial31=bsxfun(@plus,Sb31,(dev31g-dev31)); trial32=bsxfun(@plus,Sb32,(dev32g-dev32));trial33=bsxfun(@plus,Sb33,(dev33g-dev33));
    trialtensor=[trial11; trial12; trial13; trial21; trial22; trial23;trial31; trial32; trial33];
    normtrial=sqrt(sum(trialtensor.^2));
    yield=y-lam*hydro;
    eta=bsxfun(@minus,bsxfun(@times,normtrial/yield,s),1); %1*25
    eta(eta<0)=0;
    
    Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
    Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
    Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
    %1*25 for each Sb element
    Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
    normSb=sqrt(sum(Sbtensor.^2)); %sum(a) sums all the colume
    
    Ws=(bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))<=0).*...
        (0)+...
        (bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))>0).*...
        ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s)),yield),s)));
    
    W= sum(Ws);
    
    hold on;
    W_change_alp=plot ((n+1),W,'LineStyle', 'none','LineWidth', 1, 'Marker', 's', 'MarkerSize',8, ...
        'MarkerEdgeColor',  [255 193 37]/255, 'MarkerFaceColor' , [255 193 37]/255);
    Smax1=plot ((n+1),sign(stress11)*Smax/4.2e5,'LineStyle', 'none','LineWidth', 1, 'Marker', '^', 'MarkerSize',14, ...
        'MarkerEdgeColor',  'm', 'MarkerFaceColor' , 'none');
    sequence=((Smax*yield^-1)*(1-Smax*yield^-1)^-1)^b;
    sequence(sequence<0)=0;
    alp(n+1)=1-a*sequence;
    D2(n+1)=D2(n)+D2(n)^alp(n+1)*W/W0;
    n=n+1;
end
toc;
n
t=n/stepnumber*1/f;
disp(['Time to failure is ' num2str(t) ' s.']);
Nf=n*stepnumber^-1;
disp(['Cycles to failure is ' num2str(Nf) ' cycles.']);
disp(['Mean alpha to give fix alpha value is ' num2str(mean(alp)) ' .']);




% %---------------------in loop 2 scales plot settings-----------------------------
grid on;
grid minor;
set(gca ,'FontSize',30);
hXLabel = xlabel('Time step(200 steps in one cycle)' ,'Fontsize' ,30);
hYLabel = ylabel('Dissipated energy per timestep', 'Fontsize' ,30);
hTitle = title('Dissipated energy comparison of 3 methods' ,'Fontsize' ,30);
set(hTitle, 'FontSize', 30, 'FontWeight' , 'bold')
hLegend=legend([Smax1,W_cyc,W_fix_alp,W_change_alp],'S_{max}/4.2e5(scaled to plot)','W with analytical method(fixed \alpha=0.78)','W with numerical method(fixed \alpha=0.78)',...
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


saveas(gcf,'F:\Git\Doctor_thesis_Zepeng\figuresW_3methods.png');

% saveas(gcf,'F:\Git\Doctor_thesis_Zepeng\figuresW_3methods_enlarge.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('I have finished');
%mail2me('job finished',['Elapsed time is ' num2str(toc) ' seconds. Real test time is ' testtime ' seconds. Number of points to failure is ' NF ' points.']);
