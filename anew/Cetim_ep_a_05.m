% profile on
% profile off
% profile viewer
clear;clc;
close all
format long e
dbstop if error

[time,cycles,force,forcecommande,displacement]=textread('F:\Git\Cetim\ep_a_05\Acqui_CV.txt','%s%s%s%s%s','headerlines',5);
force=str2double(strrep(force,',','.')); %in unit (KN)
area=25.9695e-6; %sample section area (square milimeter)
stress11=1000*force.*area^-1; %in unit (MPa)
stress12= zeros(size(force));
stress13= zeros(size(force));
stress21=stress12;
stress22= zeros(size(force));
stress23= zeros(size(force));
stress31=stress13;
stress32=stress23;
stress33= zeros(size(force));
time=str2double(strrep(time,',','.')); 
step=time(length(time))/(length(time)-5);

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
%%
y=230e6;           %macroscopic yield stress
E=72e9;              %Young's modulus
k=6e8;                 %hardening parameter
nu=0.3;                     %poisson's ratio
sigu=320e6;             %ultimite stress
%---------------------Verified parameters in constant loading case-----------------------------
b=3.928;                    %weakening scales distribution exponent
WF=3.219e9;             %dissipated energy to failure per unit volume
%---------------------Verified parameters in random loading case-----------------------------
n0=13;                   %number of initial local defects
lam=0.3;             %hydrostatic pressure sensitivity
a=2;                  %sensitivity of sequence effect
gam=6;    %material parameter from Chaboche law(Wohler curve exponent)

%---------------------Vecterization-----------------------------
D=0;                    %initial damage
n=1;                      %initial recording point
t=n*step;
G = 0;
yield = zeros(size(force)); %Pre-allocate memory for vectors
D = zeros(size(force));
Smax = zeros(size(force));
normSb  = zeros(length(force),length(x)); %Pre-allocate memory for tensors
normtrial = zeros(length(force),length(x)); %Pre-allocate memory for tensors

%---------------------to get the the first Sb-----------------------------
hydro=1/3*sum(stress11(1)+stress22(1)+stress33(1));
yield(1)=y-lam*hydro; %macro yield strength considering mean stress effect
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
s= (x/2+1/2).^(1/(1-b)); %1*64 weak scale
eta=bsxfun(@minus,bsxfun(@times,normtrial(1,1:length(x))/yield(1),s),1); %1*64
eta(eta<0)=0;

Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
%1*64 for each Sb element
Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
normSb(1,:)=sqrt(sum(Sbtensor.^2));
Ws=(bsxfun(@minus,normtrial(1,1:length(x)),bsxfun(@rdivide, yield(1),s))<=0).*...
    (0)+...
    (bsxfun(@minus,normtrial(1,1:length(x)),bsxfun(@rdivide, yield(1),s))>0).*...
    ((E-k)*(1+nu)/(2*E*(E+k*nu))*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial(1,1:length(x)),bsxfun(@rdivide, yield(1),s)),yield(1)),s)));
W= sum(Ws);
alp=1-(1-Smax(1)/yield(1))^a;
G = G+n0*(1-alp)*(gam+1)*W/WF;
D(1)=1-(1-G.^(1/(1-alp))).^(1/(gam + 1));
%%

tic;
while G<1
    hydro=1/3*sum(stress11(n)+stress22(n)+stress33(n));
    dev1=[stress11(n) stress12(n) stress13(n);stress21(n) stress22(n) stress23(n);stress31(n) stress32(n) stress33(n)]-hydro*eye(3);
    dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
    dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
    dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
    
    hydro=1/3*sum(stress11(n+1)+stress22(n+1)+stress33(n+1));
    yield(n+1)=y-lam*hydro; %macro yield strength considering mean stress effect
    
    devn=[stress11(n+1) stress12(n+1) stress13(n+1);stress21(n+1) stress22(n+1) stress23(n+1);stress31(n+1) stress32(n+1) stress33(n+1)]-hydro*eye(3);
    dev11g=devn(1,1); dev12g=devn(1,2); dev13g=devn(1,3);
    dev21g=devn(2,1); dev22g=devn(2,2); dev23g=devn(2,3);
    dev31g=devn(3,1); dev32g=devn(3,2); dev33g=devn(3,3);
    Smax(n+1)=norm(devn,'fro');
    trial11=bsxfun(@plus,Sb11,(dev11g-dev11)); trial12=bsxfun(@plus,Sb12,(dev12g-dev12));trial13=bsxfun(@plus,Sb13,(dev13g-dev13));
    trial21=bsxfun(@plus,Sb21,(dev21g-dev21)); trial22=bsxfun(@plus,Sb22,(dev22g-dev22));trial23=bsxfun(@plus,Sb23,(dev23g-dev23));
    trial31=bsxfun(@plus,Sb31,(dev31g-dev31)); trial32=bsxfun(@plus,Sb32,(dev32g-dev32));trial33=bsxfun(@plus,Sb33,(dev33g-dev33));
    trialtensor=[trial11; trial12; trial13; trial21; trial22; trial23;trial31; trial32; trial33];
    normtrial(n+1,:)=sqrt(sum(trialtensor.^2));
    
    eta=bsxfun(@minus,bsxfun(@times,normtrial(n+1,:)/yield(n+1),s),1); %1*64
    eta(eta<0)=0;
    
    Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
    Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
    Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
    %1*64 for each Sb element
    Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
    
    normSb(n+1,:)=sqrt(sum((Sbtensor.^2)));
    
    Ws=(bsxfun(@minus,normtrial(n+1,:),bsxfun(@rdivide, yield(n+1),s))<=0).*...
        (0)+...
        (bsxfun(@minus,normtrial(n+1,:),bsxfun(@rdivide, yield(n+1),s))>0).*...
        ((E-k)*(1+nu)/(2*E*(E+k*nu))*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial(n+1,:),bsxfun(@rdivide, yield(n+1),s)),yield(n+1)),s)));
    W= sum(Ws);
    alp=1-(1-Smax(n+1)/yield(n+1))^a;
    G = G+n0*(1-alp)*(gam+1)*W/WF;
    D(n+1)=1-(1-G.^(1/(1-alp))).^(1/(gam + 1));
    t=n*step;
    
    %         hold on;
    %         yield1=plot (n+1,yield(n+1)*s(60).^-1, 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
    %             'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'c');
    %         Trial1=plot (n+1,normtrial(n+1,60),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 10, ...
    %             'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    %         Sb1=plot (n+1,normSb(n+1,60),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 10, ...
    %             'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    
    %         yield8=plot (n+1,yield(n+1)*s(64).^-1,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 8, ...
    %             'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b');
    %         Trial8=plot (n+1,normtrial(n+1,64),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 10, ...
    %             'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor',[1 0.5 0]);
    %         Sb8=plot (n+1,normSb(n+1,64),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 10, ...
    %             'MarkerEdgeColor','k', 'MarkerFaceColor','k');
    
   n=n+1;
end
toc;
disp(['Number of test points is ' num2str(n) ' points.']);
disp(['Numerical test time is ' num2str(toc) ' seconds.']);
disp(['Reallife test time is ' num2str(t) ' seconds.']);
% %---------------------Plot stress-----------------------------
figure(1);
stress=plot((1:10000)*step,stress11(1:10000),'LineWidth', 2);
grid on;
grid minor;
set(gca ,'FontSize',25);
hXLabel = xlabel('t(s)' ,'Fontsize' ,25);
hTitle =title('Stress evlolution of Ep\_05 random load' ,'Fontsize' ,25);
hYLabel =ylabel('Stress', 'Fontsize' ,25);
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
set(gcf, 'PaperPosition', [0 0 1080 608]); %set(gcf,'PaperPosition',[left,bottom,width,height])
saveas(gcf,'F:\Git\Anew\figures\ep05_stress.png');

% %---------------------Plot Damage evolution-----------------------------
figure(2);
DamageN=plot ((1:n)*step,D(1:n),'LineStyle', 'none','LineWidth', 2, 'Marker', 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor',  'r' , 'MarkerFaceColor' ,'none');
% % ---------------------plot settings-----------------------------
grid on;
grid minor;
set(gca ,'FontSize',25);
hXLabel = xlabel('t(s)' ,'Fontsize' ,25);
hTitle =title('Damage evolution of Ep\_05 random load' ,'Fontsize' ,25);
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
set(gcf, 'PaperPosition', [0 0 1080 608]); %set(gcf,'PaperPosition',[left,bottom,width,height])
 saveas(gcf,'F:\Git\Anew\figures\ep05_damage.png');

% sp=actxserver('SAPI.SpVoice');
% sp.Speak('I finished all the work finally. oh la la');
%
%

