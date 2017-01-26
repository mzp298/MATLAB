profile on
clear;clc;
dbstop if error
format long e

load('FX_RAVG.mat');
signal.data=double(signal.data);
forcex= transpose(signal.data);
load('FY_RAVG.mat');
signal.data=double(signal.data);
forcey= transpose(signal.data);
load('FZ_RAVG.mat');
signal.data=double(signal.data);
forcez= transpose(signal.data);
copy=1.2;
forcex=repmat(forcex,copy,1);
forcey=repmat(forcey,copy,1);
forcez=repmat(forcez,copy,1);

% %---------------------Arithmetic sequence between every recorded points---------------------
% ari=5; %insert (ari-1) points between the two limits
% for i=2:(copy*802805)
%     %force(1+ari*(i-1):1+ari*i)=linspace(forceorigin(i),forceorigin(i+1),ari+1);
%     forcelx(1+ari*(i-2):1+ari*(i-1))=linspace(forcex(i-1),forcex(i),ari+1);
%     forcely(1+ari*(i-2):1+ari*(i-1))=linspace(forcey(i-1),forcey(i),ari+1);
%     forcelz(1+ari*(i-2):1+ari*(i-1))=linspace(forcez(i-1),forcez(i),ari+1);
% end;
% %  ari*(i-1)+1 ;%the number of points
ari=1;
forcelx=forcex;
forcely=forcey;
forcelz=forcez;

%------------------------build the stress tensor---------------------
A=2.25e-4; %size of the gauge in square meter
% A=1/6e4;
cx=10;
cy=60;
thetax=0.5;
thetay=0.6;
phix=0.3;
phiy=0.4;
stress11=1/A*(forcelz+cx*forcelx*cos(thetax)^2+cy*forcely*cos(thetay)^2);
stress12=1/A*(cx*forcelx*cos(thetax)*sin(thetax)*cos(phix)+cy*forcely*cos(thetay)*sin(thetay)*cos(phiy));
stress13=1/A*(cx*forcelx*cos(thetax)*sin(thetax)*sin(phix)+cy*forcely*cos(thetay)*sin(thetay)*sin(phiy));
stress21=stress12;
stress22=1/A*(cx*forcelx*sin(thetax)^2*cos(phix)^2+cy*forcely*sin(thetay)^2*cos(phiy)^2);
stress23=1/A*(cx*forcelx*sin(thetax)^2*cos(phix)*sin(phix)+cy*forcely*sin(thetay)^2*cos(phiy)*sin(phiy));
stress31=stress13;
stress32=stress23;
stress33=1/A*(cx*forcelx*sin(thetax)^2*sin(phix)^2+cy*forcely*sin(thetay)^2*sin(phiy)^2);
% [max(stress11) max(stress12) max(stress13);
% max(stress12) max(stress22) max(stress23);
% max(stress23) max(stress13) max(stress33);]
% [mean(stress11) mean(stress12) mean(stress13);
% mean(stress12) mean(stress22) mean(stress23);
% mean(stress23) mean(stress13) mean(stress33);]

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
y=6.38e8;           %macroscopic yield stress
lam=0.5;             %hydrostatic pressure sensitivity
E=2e11;              %Young¡¯s modulus
k=6e8;                 %hardening parameter
b=3;                    %weakening scales distribution exponent
nu=0.3;                     %poisson's ratio
tt=2e8;                 %torsion fatigue limit
ff=2.5e8;              %bending fatigue limit
ac=(tt-ff/sqrt(3))/(ff/3); %crossland criterial constant
bc=tt;                     %crossland criterial constant
sigu=8e8;             %ultimite stress
gam=b+1;    %material parameter from Chaboche law(Wohler curve exponent)
samplerate=256;   %recorded samples per second
delta=(b+1)/(b-1); %yield stress degradation sensitivity with D
n0=3;                   %number of initial local defects
%---------------------Vecterization-----------------------------

WF=3e8;             %dissipated energy to failure per unit volume
D=0;                    %initial damage
n=1;                      %initial recording point
step=1/samplerate/ari;
t=n*step;
G = 0;
yield = zeros(size(forcelx)); %Pre-allocate memory for vectors
D = yield;
normSb  = zeros(length(forcelx),length(x)); %Pre-allocate memory for tensors
normtrial = zeros(length(forcelx),length(x)); %Pre-allocate memory for tensors

%---------------------to get the the first Sb-----------------------------
hydro=1/3*sum(stress11(1)+stress22(1)+stress33(1));
M=(1-3*hydro/sigu); %M function in Chaboche model
M(M<0)=0;
yield(1)=y*(1-D(1))^delta*M; %macro yield strength considering mean stress effect
dev1=[stress11(1) stress12(1) stress13(1);stress21(1) stress22(1) stress23(1);stress31(1) stress32(1) stress33(1)]-hydro*eye(3);
dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
Smax=norm(dev1,'fro');

trial11=dev11; trial12=dev12; trial13=dev13;
trial21=dev21; trial22=dev22; trial23=dev23;
trial31=dev31; trial32=dev32; trial33=dev33;
trialtensor=[trial11; trial12; trial13; trial21; trial22; trial23;trial31; trial32; trial33];
normtrial(1,:)=sqrt(sum(trialtensor.^2));
s= (x/2+1/2).^(1/(1-b)); %1*25 weak scale
eta=bsxfun(@minus,bsxfun(@times,normtrial(1,1:length(x))/yield(1),s),1); %1*25
eta(eta<0)=0;

Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
%1*25 for each Sb element
Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
normSb(1,:)=sqrt(sum(Sbtensor.^2));
Ws=(bsxfun(@minus,normtrial(1,1:length(x)),bsxfun(@rdivide, yield(1),s))<=0).*...
    (0)+...
    (bsxfun(@minus,normtrial(1,1:length(x)),bsxfun(@rdivide, yield(1),s))>0).*...
    ((E-k)*(1+nu)/(2*E*(E+k*nu))*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial(1,1:length(x)),bsxfun(@rdivide, yield(1),s)),yield(1)),s)));
W= sum(Ws);
alp=1-W*((1-Smax/y)*sigu)^-1;
G = G+n0*(1-alp)*(gam+1)*W/WF; %1.322163316411401e-03
D(1)=1-(1-G.^(1/(1-alp))).^(1/(gam + 1));

tic;
while G<1
    hydro=1/3*sum(stress11(n)+stress22(n)+stress33(n));
    dev1=[stress11(n) stress12(n) stress13(n);stress21(n) stress22(n) stress23(n);stress31(n) stress32(n) stress33(n)]-hydro*eye(3);
    dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
    dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
    dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
    
    hydro=1/3*sum(stress11(n+1)+stress22(n+1)+stress33(n+1));
    M=(1-hydro/sigu); %M function in Chaboche model
    M(M<0)=0;
    yield(n+1)=y*(1-D(n))^delta*M; %macro yield strength considering mean stress effect
    
    devn=[stress11(n+1) stress12(n+1) stress13(n+1);stress21(n+1) stress22(n+1) stress23(n+1);stress31(n+1) stress32(n+1) stress33(n+1)]-hydro*eye(3);
    dev11g=devn(1,1); dev12g=devn(1,2); dev13g=devn(1,3);
    dev21g=devn(2,1); dev22g=devn(2,2); dev23g=devn(2,3);
    dev31g=devn(3,1); dev32g=devn(3,2); dev33g=devn(3,3);
    Smax=norm(devn,'fro');
    trial11=bsxfun(@plus,Sb11,(dev11g-dev11)); trial12=bsxfun(@plus,Sb12,(dev12g-dev12));trial13=bsxfun(@plus,Sb13,(dev13g-dev13));
    trial21=bsxfun(@plus,Sb21,(dev21g-dev21)); trial22=bsxfun(@plus,Sb22,(dev22g-dev22));trial23=bsxfun(@plus,Sb23,(dev23g-dev23));
    trial31=bsxfun(@plus,Sb31,(dev31g-dev31)); trial32=bsxfun(@plus,Sb32,(dev32g-dev32));trial33=bsxfun(@plus,Sb33,(dev33g-dev33));
    trialtensor=[trial11; trial12; trial13; trial21; trial22; trial23;trial31; trial32; trial33];
    normtrial(n+1,:)=sqrt(sum(trialtensor.^2));
    
    eta=bsxfun(@minus,bsxfun(@times,normtrial(n+1,:)/yield(n+1),s),1); %1*25
    eta(eta<0)=0;
    
    Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
    Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
    Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
    %1*25 for each Sb element
    Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
    
    normSb(n+1,:)=sqrt(sum((Sbtensor.^2)));
    
    Ws=(bsxfun(@minus,normtrial(n+1,:),bsxfun(@rdivide, yield(n+1),s))<=0).*...
        (0)+...
        (bsxfun(@minus,normtrial(n+1,:),bsxfun(@rdivide, yield(n+1),s))>0).*...
        ((E-k)*(1+nu)/(2*E*(E+k*nu))*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial(n+1,:),bsxfun(@rdivide, yield(n+1),s)),yield(n+1)),s)));
    W= sum(Ws);
    alp=1-W*((1-Smax/y)*sigu)^-1;
    G = G+n0*(1-alp)*(gam+1)*W/WF
    D(n+1)=1-(1-G.^(1/(1-alp))).^(1/(gam + 1));
    t=n*step;
    %            hold on;
    %         yield1=plot (n,yield(n)*s(1).^-1, 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 10, ...
    %             'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'c');
    %         Trial1=plot (n,normtrial(n,1),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 10, ...
    %             'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    %         Sb1=plot (n,normSb(n,1),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 10, ...
    %             'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    %         yield61=plot (n,yield(n)*s(61).^-1,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 10, ...
    %             'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b');
    %         Trial61=plot (n,normtrial(n,61),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize',10, ...
    %             'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor',[1 0.5 0]);
    %         Sb61=plot (n,normSb(n,61),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 10, ...
    %             'MarkerEdgeColor','k', 'MarkerFaceColor','k');
    
    %     hold on;
    %     DamageN=plot (t,D,'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize',10, ...
    %        'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'r');
    
    
    n=n+1;
end
toc;
disp(['Number of test points is ' num2str(n/ari+1) ' points.']);
disp(['Number of test time is ' num2str(t) ' seconds.']);
testtime=num2str(t);

%---------------------Plot Trial and Sb evolution-----------------------------
figure(1);
 hold on;
  Trial1=plot ((1:n)*step,normtrial(1:n,1),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize',12, ...
    'MarkerEdgeColor','r', 'MarkerFaceColor','none');
  Trial61=plot ((1:n)*step,normtrial(1:n,61),'LineStyle', 'none','LineWidth', 1,'Marker', 's', 'MarkerSize',12, ...
    'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor','none');
   Sb1=plot ((1:n)*step,normSb(1:n,1),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 8, ...
    'MarkerEdgeColor','none', 'MarkerFaceColor',[96 96 96]/255);
   Sb61=plot ((1:n)*step,normSb(1:n,61),'LineStyle', 'none','LineWidth', 1,'Marker', 's', 'MarkerSize',8, ...
    'MarkerEdgeColor','none', 'MarkerFaceColor','g');
  yield1=plot ((1:n)*step,yield(1:n)*s(1).^-1,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 8, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceColor','b');
  yield61=plot ((1:n)*step,yield(1:n)*s(61).^-1, 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
    'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'c');
%---------------------plot settings-----------------------------
grid on;
grid minor;
set(gca ,'FontSize',25);
hXLabel = xlabel('t(s)' ,'Fontsize' ,25);
hTitle = title('Microscopic stress evolution at 2 scales' ,'Fontsize' ,25);
hYLabel = ylabel('Stress(Pa)', 'Fontsize' ,25);
hLegend=legend([yield1,Sb1,Trial1,yield61,Sb61,Trial61],'(\sigma_y-\lambda\Sigma_H)/s_1     at scale s_1','||S-b||              at scale s_1',...
    '||S-b||_{trial}         at scale s_1', '(\sigma_y-\lambda\Sigma_H)/s_61     at scale s_{61}','||S-b||              at scale s_{61}','||S-b||_{trial}         at scale s_{61}');
set([hLegend, gca], 'FontSize', 25)
% Adjust font
set(gca, 'FontName', 'Helvetica')
set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hXLabel, hYLabel], 'FontSize', 25)
set(hTitle, 'FontSize', 25, 'FontWeight' , 'bold')
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
set(gcf, 'PaperPosition', [0 0 1920 1080]); %set(gcf,'PaperPosition',[left,bottom,width,height])

% saveas(gcf,'trialreal3d yield evolve with D.png');
%
% %---------------------Plot Damage evolution-----------------------------
figure(2);
DamageN=plot ((1:n)*step,D(1:n),'LineStyle', 'none','LineWidth', 2, 'Marker', 'o', 'MarkerSize', 10, ...
   'MarkerEdgeColor',  'r' , 'MarkerFaceColor' ,'none');
% % ---------------------plot settings-----------------------------
grid on;
grid minor;
set(gca ,'FontSize',25);
hXLabel = xlabel('t(s)' ,'Fontsize' ,25);
 hTitle =title('Damage evolution under multidimensional stress' ,'Fontsize' ,25);
 hYLabel =ylabel('D', 'Fontsize' ,25);
 % Adjust font
set(gca, 'FontName', 'Helvetica')
set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hXLabel, hYLabel], 'FontSize', 25)
set(hTitle, 'FontSize', 25, 'FontWeight' , 'bold')
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
set(gcf, 'PaperPosition', [0 0 1920 1080]); %set(gcf,'PaperPosition',[left,bottom,width,height])
% saveas(gcf,'damage3d yield evolve with D.png');

% sp=actxserver('SAPI.SpVoice');
% sp.Speak('I finished all the work finally. oh la la');
%  mail2me('job finished',['Elapsed time is ' num2str(toc) ' seconds. Real test time is ' testtime ' seconds. Number of test points is ' num2str(n/ari+1) ' points.']);
%
profile off
profile viewer
