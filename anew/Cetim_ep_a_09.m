clear;clc;
close all
format long e

fid = fopen('F:\Git\Cetim\ep_a_09\Acqui_CV.txt');
% fid = fopen('/home/ma/MATLAB/Cetim/ep_a_09/Acqui_CV.txt');
[force]=textscan(fid,'%*s%*s%s%*s%*s','headerlines',5);
area=2.29*10*1e-6; %meter square
stress11=1000*str2double(strrep(force{1,1},',','.')).*area^-1; %Pa
repetition=xlsread('F:\Git\MATLAB\anew\tuning.xlsx',1,'G11');
nbrep=xlsread('F:\Git\MATLAB\anew\tuning.xlsx',1,'F11');
copy=nbrep+300;
stress11=repmat(stress11(1:repetition),copy,1);
clear force;


y=230e6;           %macroscopic yield stress
E=72e9;              %Young's modulus
k=6e8;                 %hardening parameter
nu=0.3;                     %poisson's ratio
sigu=320e6;             %ultimite stress
%---------------------Verified parameters in constant loading case-----------------------------
b=4.425;                    %weakening scales distribution exponent
WF=6.14e8;             %dissipated energy to failure per unit volume
%---------------------Verified parameters in random loading case-----------------------------
n0=1;                   %number of initial local defects
lam=0.3;             %hydrostatic pressure sensitivity
a=0.7;                  %sensitivity of sequence effect
gam=6;    %material parameter from Chaboche law(Wohler curve exponent)

%---------------------Vecterization-----------------------------
D=0;                    %initial damage
n=1;                      %initial recording point
G = 0;

x= [0.999305042	0.996340117	0.991013371	0.983336254	0.973326828	0.9610088	0.946411375	0.929569172	0.910522137...
    0.889315446	0.865999398	0.840629296	0.813265315	0.783972359	0.752819907	0.71988185	0.685236313	0.648965471...
    0.611155355	0.571895646	0.531279464	0.489403146	0.446366017	0.402270158	0.357220158	0.311322872	0.264687162...
    0.217423644	0.16964442	0.121462819	0.072993122	0.024350293	-0.024350293	-0.072993122	-0.121462819	-0.16964442...
    -0.217423644	-0.264687162	-0.311322872	-0.357220158	-0.402270158	-0.446366017	-0.489403146	-0.531279464...
    -0.571895646	-0.611155355	-0.648965471	-0.685236313	-0.71988185	-0.752819907	-0.783972359	-0.813265315...
    -0.840629296	-0.865999398	-0.889315446	-0.910522137	-0.929569172	-0.946411375	-0.9610088	-0.973326828...
    -0.983336254	-0.991013371	-0.996340117	-0.999305042];
weight=[0.001783281	0.004147033	0.006504458	0.00884676	0.011168139	0.013463048	0.01572603	0.017951716	0.020134823...
    0.022270174	0.024352703	0.02637747	0.028339673	0.030234657	0.032057928	0.033805162	0.035472213	0.037055129	0.038550153...
    0.039953741	0.041262563	0.042473515	0.043583725	0.044590558	0.045491628	0.046284797	0.046968183	0.047540166	0.047999389...
    0.048344762	0.048575467	0.048690957	0.048690957	0.048575467	0.048344762	0.047999389	0.047540166	0.046968183	0.046284797...
    0.045491628	0.044590558	0.043583725	0.042473515	0.041262563	0.039953741	0.038550153	0.037055129	0.035472213	0.033805162...
    0.032057928	0.030234657	0.028339673	0.02637747	0.024352703	0.022270174	0.020134823	0.017951716	0.01572603	0.013463048...
    0.011168139	0.00884676	0.006504458	0.004147033	0.001783281];
% x=xlsread('Gauss-Legendre Quadrature','Sheet1','b1:z1');
% weight=xlsread('Gauss-Legendre Quadrature','Sheet1','b2:z2');

%---------------------to get the the first Sb-----------------------------
hydro=1/3*sum(stress11(1)+0+0);
yield=y-lam*hydro; %macro yield strength considering mean stress effect
dev1=[stress11(1) 0 0;0 0 0;0 0 0]-hydro*eye(3);
dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
Smax=norm(dev1,'fro');

trial11=dev11; trial12=dev12; trial13=dev13;
trial21=dev21; trial22=dev22; trial23=dev23;
trial31=dev31; trial32=dev32; trial33=dev33;
trialtensor=[trial11; trial12; trial13; trial21; trial22; trial23;trial31; trial32; trial33];
normtrial=sqrt(sum(trialtensor.^2));
s= (x/2+1/2).^(1/(1-b)); %1*64 weak scale
eta=bsxfun(@minus,bsxfun(@times,normtrial(1)/yield,s),1); %compare normtrial with yield/s 
eta(eta<0)=0;

Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
%1*64 for each Sb element
Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
normSb=sqrt(sum(Sbtensor.^2));
Ws=(bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))<=0).*...
    (0)+...
    (bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))>0).*...
    ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s)),yield),s)));
W= sum(Ws);
sequence=1-Smax/yield;
sequence(sequence<0)=0;
alp=1-a*sequence;
G = G+n0*(1-alp)*(gam+1)*W/WF;
D(1)=1-(1-G.^(1/(1-alp))).^(1/(gam + 1));
%%

tic;
while G<1
    hydro=1/3*sum(stress11(n)+0+0);
    dev1=[stress11(n) 0 0;0 0 0;0 0 0]-hydro*eye(3);
    dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
    dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
    dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
    
    hydro=1/3*sum(stress11(n+1)+0+0);
    yield=y-lam*hydro; %macro yield strength considering mean stress effect
    devn=[stress11(n+1) 0 0;0 0 0;0 0 0]-hydro*eye(3);
    dev11g=devn(1,1); dev12g=devn(1,2); dev13g=devn(1,3);
    dev21g=devn(2,1); dev22g=devn(2,2); dev23g=devn(2,3);
    dev31g=devn(3,1); dev32g=devn(3,2); dev33g=devn(3,3);
    Smax=norm(devn,'fro');
    
    trial11=bsxfun(@plus,Sb11,(dev11g-dev11)); trial12=bsxfun(@plus,Sb12,(dev12g-dev12));trial13=bsxfun(@plus,Sb13,(dev13g-dev13));
    trial21=bsxfun(@plus,Sb21,(dev21g-dev21)); trial22=bsxfun(@plus,Sb22,(dev22g-dev22));trial23=bsxfun(@plus,Sb23,(dev23g-dev23));
    trial31=bsxfun(@plus,Sb31,(dev31g-dev31)); trial32=bsxfun(@plus,Sb32,(dev32g-dev32));trial33=bsxfun(@plus,Sb33,(dev33g-dev33));
    trialtensor=[trial11; trial12; trial13; trial21; trial22; trial23;trial31; trial32; trial33];
    normtrial=sqrt(sum(trialtensor.^2));
    
    eta=bsxfun(@minus,bsxfun(@times,normtrial/yield,s),1); %1*64
    eta(eta<0)=0;
    
    Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
    Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
    Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
    %1*64 for each Sb element
    Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
    normSb=sqrt(sum((Sbtensor.^2)));
   
    Ws=(bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))<=0).*...
    (0)+...
    (bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))>0).*...
    ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s)),yield),s)));
    W= sum(Ws);
    sequence=1-Smax/yield;
    sequence(sequence<0)=0;
    alp=1-a*sequence;

    G = G+n0*(1-alp)*(gam+1)*W/WF;
    D(n+1)=1-(1-G.^(1/(1-alp))).^(1/(gam + 1));

   n=n+1;
end
toc;
disp(['Number of test points is ' num2str(n) ' points.']);
% 
% % %---------------------Plot stress-----------------------------
% figure(1);
% stress=plot(1:10000,stress11(1:10000),'LineWidth', 2);
% grid on;
% grid minor;
% set(gca ,'FontSize',25);
% hXLabel = xlabel('t(s)' ,'Fontsize' ,25);
% hTitle =title('Stress evlolution of EP\_09 random load' ,'Fontsize' ,25);
% hYLabel =ylabel('Stress', 'Fontsize' ,25);
% % Adjust font
% set(gca, 'FontName', 'Helvetica')
% set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')
% set([hXLabel, hYLabel], 'FontSize', 25)
% set(hTitle, 'FontSize', 25, 'FontWeight' , 'bold')
% % Adjust axes properties
% set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
%     'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
%     'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
%     'LineWidth', 1)
% set(gcf,'color','w'); %set figure background transparent
% set(gca,'color','w'); %set axis transparent
% % Maximize print figure
% set(gcf,'outerposition',get(0,'screensize'));
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
% set(gcf, 'PaperPosition', [0 0 1080 608]); %set(gcf,'PaperPosition',[left,bottom,width,height])
% saveas(gcf,'F:\Git\Anew\figures\ep_09_stress.png');

% % %---------------------Plot Damage evolution-----------------------------
% clear stress11;
% figure(2);
% DamageN=plot ((1:n),D(1:n),'LineStyle', 'none','LineWidth', 2, 'Marker', 'o', 'MarkerSize', 10, ...
%     'MarkerEdgeColor',  'r' , 'MarkerFaceColor' ,'none');
% % % ---------------------plot settings-----------------------------
% grid on;
% grid minor;
% set(gca ,'FontSize',25);
% hXLabel = xlabel('t(s)' ,'Fontsize' ,25);
% hTitle =title('Damage evolution of EP\_09 random load' ,'Fontsize' ,25);
% hYLabel =ylabel('D', 'Fontsize' ,25);
% % Adjust font
% set(gca, 'FontName', 'Helvetica')
% set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')
% set([hXLabel, hYLabel], 'FontSize', 25)
% set(hTitle, 'FontSize', 25, 'FontWeight' , 'bold')
% % Adjust axes properties
% set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
%     'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
%     'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
%     'LineWidth', 1)
% set(gcf,'color','w'); %set figure background transparent
% set(gca,'color','w'); %set axis transparent
% % Maximize print figure
% set(gcf,'outerposition',get(0,'screensize'));
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
% set(gcf, 'PaperPosition', [0 0 1080 608]); %set(gcf,'PaperPosition',[left,bottom,width,height])
% saveas(gcf,'F:\Git\Anew\figures\ep_09_damage.png');
% 
xlswrite('F:\Git\MATLAB\anew\tuning.xlsx',a,1,'C11');
xlswrite('F:\Git\MATLAB\anew\tuning.xlsx',n,1,'D11');
% xlswrite('/home/ma/MATLAB/tuning.xlsx',n0,1,'C11');
% xlswrite('/home/ma/MATLAB/tuning.xlsx',n0,1,'D11');