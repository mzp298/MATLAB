
% Program to get the Gauss-Legendre Quadrature results (Vectorized)
clear;clc;close all;

dbstop if error
format long e

[x]= [0.999305042	0.996340117	0.991013371	0.983336254	0.973326828	0.9610088	0.946411375	0.929569172	0.910522137...
    0.889315446	0.865999398	0.840629296	0.813265315	0.783972359	0.752819907	0.71988185	0.685236313	0.648965471...
    0.611155355	0.571895646	0.531279464	0.489403146	0.446366017	0.402270158	0.357220158	0.311322872	0.264687162...
    0.217423644	0.16964442	0.121462819	0.072993122	0.024350293	-0.024350293	-0.072993122	-0.121462819	-0.16964442...
    -0.217423644	-0.264687162	-0.311322872	-0.357220158	-0.402270158	-0.446366017	-0.489403146	-0.531279464...
    -0.571895646	-0.611155355	-0.648965471	-0.685236313	-0.71988185	-0.752819907	-0.783972359	-0.813265315...
    -0.840629296	-0.865999398	-0.889315446	-0.910522137	-0.929569172	-0.946411375	-0.9610088	-0.973326828...
    -0.983336254	-0.991013371	-0.996340117	-0.999305042];
[weight]=[0.001783281	0.004147033	0.006504458	0.00884676	0.011168139	0.013463048	0.01572603	0.017951716	0.020134823...
    0.022270174	0.024352703	0.02637747	0.028339673	0.030234657	0.032057928	0.033805162	0.035472213	0.037055129	0.038550153...
    0.039953741	0.041262563	0.042473515	0.043583725	0.044590558	0.045491628	0.046284797	0.046968183	0.047540166	0.047999389...
    0.048344762	0.048575467	0.048690957	0.048690957	0.048575467	0.048344762	0.047999389	0.047540166	0.046968183	0.046284797...
    0.045491628	0.044590558	0.043583725	0.042473515	0.041262563	0.039953741	0.038550153	0.037055129	0.035472213	0.033805162...
    0.032057928	0.030234657	0.028339673	0.02637747	0.024352703	0.022270174	0.020134823	0.017951716	0.01572603	0.013463048...
    0.011168139	0.00884676	0.006504458	0.004147033	0.001783281];
% [x]=xlsread('Gauss-Legendre Quadrature','Sheet1','b1:z1');
% [weight]=xlsread('Gauss-Legendre Quadrature','Sheet1','b2:z2');

E=191e9;               %Young's modulus
nu=0.38;                 %poisson's ratio
k=1e9;                  %hardening parameter
b=1.5;                      %weakening scales distribution exponent (between 1 and 2)
a=0.1;
y=1080e6;            %macroscopic yield stress
sigu=1200e6;             %ultimite stress
ff=690e6;              %bending fatigue limit
tt=428e6;                  %torsion fatigue limit

load=6e8;            %cyclic load
loadtensor= [load 0 0;0 0 0;0 0 0];
stepnumber=100;        %devide one cycle in 200 parts
f=50;                            %frequency of load
steptime=1/f/stepnumber;
W0=5e8;             %dissipated energy to failure per unit volume
lam=0.7;               %hydrostatic pressure sensitivity
% m=3e8;                   % mean stress
m=0;                   % mean stress
hydrofix=1/3*(sum(diag(loadtensor))); 
dev=loadtensor-hydrofix*eye(3); %mean stress does not change deviatoric stress!!!!!!!!


%---------------------3 Numerical method-----------------------------
D=1e-16;
n=1;       %initial recording point
%---------------------to get the the first Sb-----------------------------
stress11=m+load*sind(n*360/stepnumber);
hydro=1/3*sum(stress11(1)+0+0);
yield=y-lam*hydro; %macro yield strength considering mean stress effect
dev1=[stress11 0 0 ;0 0 0 ;0 0 0 ]-hydro*eye(3);
dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
Smax=norm(dev1,'fro');
s= (x/2+1/2).^(1/(1-b)); %1*64

trial11=dev11; trial12=dev12; trial13=dev13;
trial21=dev21; trial22=dev22; trial23=dev23;
trial31=dev31; trial32=dev32; trial33=dev33;

normtrial(1)=norm([trial11, trial12, trial13; trial21, trial22, trial23;trial31, trial32, trial33],'fro');

eta=bsxfun(@minus,bsxfun(@times,normtrial(1)/yield,s),1); %compare normtrial with yield/s
eta(eta<0)=0; %only keep normtrials which are larger than yield/s

Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
%1*64 for each Sb element
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
alp=1-a*sequence;
D=D+D^alp*W/W0;
%%

tic;
while n<2*stepnumber
    stress11=m+load*sind(n*360/stepnumber);
    hydro=1/3*sum(stress11+0+0);
    
    dev1=[stress11 0 0 ;0 0 0 ;0 0 0 ]-hydro*eye(3);
    dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
    dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
    dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
    
    stress11=m+load*sind((n+1)*360/stepnumber);
    hydro=1/3*sum(stress11+0+0);
	yield=y-lam*hydro;
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

    eta=bsxfun(@minus,bsxfun(@times,normtrial/yield,s),1); %1*64
    eta(eta<0)=0;
    
    Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
    Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
    Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
    %1*64 for each Sb element
    Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
    normSb=sqrt(sum(Sbtensor.^2)); %sum(a) sums all the colume
    
    Ws=(bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))<=0).*...
        (0)+...
        (bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))>0).*...
        ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s)),yield),s)));
    
    W= sum(Ws);
sequence=((Smax*yield^-1)*(1-Smax*yield^-1)^-1)^b;
sequence(sequence<0)=0;
alp=1-a*sequence;
D=D+D^alp*W/W0;
    
        figure(2);
            hold on;
            yield1=plot (n+1,yield*s(33).^-1, 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
                'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'c');
            yield1n= plot (n+1,-yield*s(33).^-1, 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
                'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'c');
            Trial1=plot (n+1,sign(trial11(33))*normtrial(33),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 10, ...
                'MarkerEdgeColor','r', 'MarkerFaceColor','r');
            Sb1=plot (n+1,sign(Sb11(33))*normSb(33),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 10, ...
                'MarkerEdgeColor','g', 'MarkerFaceColor','g');
            dev=plot (n+1,sign(stress11)*Smax,'LineStyle', 'none','LineWidth', 1,'Marker', 's', 'MarkerSize', 11, ...
                    'MarkerEdgeColor','none', 'MarkerFaceColor',[238 18 137]/255);
            yield8=plot (n+1,yield*s(40).^-1,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 8, ...
                'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b');
            yield8n=plot (n+1,-yield*s(40).^-1,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 8, ...
                'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b');
            Trial8=plot (n+1,sign(trial11(40))*normtrial(40),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 10, ...
                'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor',[1 0.5 0]);
            Sb8=plot (n+1,sign(Sb11(40))*normSb(40),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 10, ...
                'MarkerEdgeColor','k', 'MarkerFaceColor','k');
    
    n=n+1;
end
toc;
n
t=n/stepnumber*1/f;
disp(['Time to failure is ' num2str(t) ' s.']);
Nf=n*stepnumber^-1;
disp(['Cycles to failure is ' num2str(Nf) ' cycles.']);



% %---------------------in loop 2 scales plot settings-----------------------------
figure(2);
grid on;
grid minor;
axis([0 n -8e8 8e8]);
set(gca ,'FontSize',30);
hXLabel = xlabel('Time step(1/100 s)' ,'Fontsize' ,30);
hYLabel = ylabel('Stress(Pa)', 'Fontsize' ,30);
hTitle = title('Microscopic stress evolution at 2 scales' ,'Fontsize' ,30);
set(hTitle, 'FontSize', 30, 'FontWeight' , 'bold')
hLegend=legend([dev,yield1,Sb1,Trial1,yield8,Sb8,Trial8],'S_{max}=dev\Sigma','(\sigma_y-\lambda\Sigma_H)/s_{33}     at scale s_{33}','||S-b||                at scale s_{33}',...
    '||S-b||_{trial}           at scale s_{33}', '(\sigma_y-\lambda\Sigma_H)/s_{40}     at scale s_{40}',...
    '||S-b||                at scale s_{40}','||S-b||_{trial}           at scale s_{40}','Location','northeast');
set(hLegend, 'FontSize', 15)
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
set(gcf, 'PaperPosition', [0 0 1080 720]); %set(gcf,'PaperPosition',[left,bottom,width,height])
%saveas(gcf,'F:\Git\Anew\figures\trialsin_m.png');
%  saveas(gcf,'F:\Git\Anew\figures\trialsin_0.png');
%
sp=actxserver('SAPI.SpVoice');
sp.Speak('Fuck that I finished all this shit finally');


