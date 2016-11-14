
% Program to get the Gauss-Legendre Quadrature results (Vectorized)
clear;clc;

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

E=2e11;               %Young��s modulus
k=6e8;                  %hardening parameter
b=1.5;                      %weakening scales distribution exponent
nu=0.3;                 %poisson's ratio
tt=2e8;                  %torsion fatigue limit
ff=2.5e8;              %bending fatigue limit
ac=(tt-ff/sqrt(3))/(ff/3); %crossland criterial constant
bc=tt;                    %crossland criterial constant
sigu=8e8;             %ultimite stress
gam=b+1;              %material parameter from Chaboche law(Wohler curve exponent)
y=6.38e8;            %macroscopic yield stress
load=5e8;            %cyclic load
loadtensor= [load 0 0;0 0 0;0 0 0];
stepnumber=200;        %devide one cycle in 200 parts
f=50;                            %frequency of load
delta=(b+1)/(b-1);
alp=0.5;
WF=5e8;             %dissipated energy to failure per unit volume

% %---------------------Plot 3 methods-----------------------------
% %---------------------1 Chaboche method-----------------------------
% Dcha(1)=0;             %initial damage
% n=1;       %initial recording point
%   G = (1 - (1 - Dcha(1)).^(gam + 1)).^(1-alp);
%   m=1/3*sum(diag(loadtensor));
%   S1=loadtensor-m*diag([1,1,1]);
%   sqrj1=1/2*sqrt(1/2)*norm(S1,'fro');
%   M=12.4262*ff*(1-3*m/sigu);
%   while G<1
%   NF=1/((gam+1)*(1-alp))*(sqrj1/M)^(-gam);
%   G = G+(1-alp)*(gam + 1)/stepnumber/NF;
%   Dcha(n+1)=1-(1-G.^(1/(1-alp))).^(1/(gam + 1));
%   t=n/stepnumber*1/f;
%   n=n+1;
%   end
%   n
%     figure(2);
%     hold on;
% DamageCha=plot ((1:n),Dcha(1:n), 'LineStyle', 'none','LineWidth', 1.2, 'Marker', 'o', 'MarkerSize', 10, ...
%    'MarkerEdgeColor',  'g', 'MarkerFaceColor','none');
% %---------------------2 Cyclic load calculation-----------------------------
% Dcyc(1)=0;
% n=1;
% Gcyc = (1 - (1 - Dcyc(1)).^(gam + 1)).^(1-alp);
% Wcyc=4*(E-k)*(1+nu)*(b-1)/(E*(E+k*nu)*b*(b+1))*norm(loadtensor-(1/3*sum(diag(loadtensor)))*diag([1,1,1]),'fro').^(b+1)*y.^(1-b) ;
% while Gcyc< 1
%  Gcyc = Gcyc+(1-alp)*(gam + 1)*Wcyc/stepnumber/WF;
%  Dcyc(n+1)=1-(1-Gcyc^(1/(1-alp)))^(1/(gam + 1));
%  n=n+1;
% end
% n
%   hold on;
%  Damagecyc=plot ((1:n),Dcyc(1:n),'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
%    'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'b');
% % %---------------------Cyclic load calculation before integration-----------------------------
% % Dstar(1)=1e-15;
% % n=1;
% % while Dstar(n)<1
% %  Dstar(n+1)=Dstar(n)+(1-(1-Dstar(n))^(gam+1))^alp*(1-Dstar(n))^-(b+1)*Wcyc/stepnumber/WF;
% %  n=n+1;
% % end
% % n
% %  hold on
% %  Damagestar=plot ((1:n),Dstar(1:n),'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 10, ...
% %    'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'r');

%---------------------3 Numerical method before integration-----------------------------
D=1e-16;       %initial damage
n=1;       %initial recording point
%---------------------to get the the first Sb-----------------------------
stress11=load*sin(2*pi/stepnumber);
m=1/3*sum(stress11+0+0);
dev1=[stress11 0 0 ;0 0 0 ;0 0 0 ]-m*diag([1,1,1]);
dev11=dev1(1,1); dev12=dev1(1,3); dev13=dev1(1,3);
dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
s= (x/2+1/2).^(1/(1-b)); %1*25

trial11=dev11; trial12=dev12; trial13=dev13;
trial21=dev21; trial22=dev22; trial23=dev23;
trial31=dev31; trial32=dev32; trial33=dev33;

normtrial(1)=norm([trial11, trial12, trial13; trial21, trial22, trial23;trial31, trial32, trial33],'fro');
yield=y*(1-D(n))^delta;
% yield=y;
eta=bsxfun(@minus,bsxfun(@times,normtrial(1)/yield,s),1); %1*25
eta(eta<0)=0;

Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
%1*25 for each Sb element
Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
normSb=sqrt(sum(Sbtensor.^2));
% existsOnGPU(normSb)
    Ws=(bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,[s]))<=0).*...
        (0)+...
        (bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,[s]))>0).*...
        ((E-k)*(1+nu)/(2*E*(E+k*nu))*bsxfun(@times,[weight],bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,[s])),yield),[s])));
    
    W= sum(Ws);
    D(1)=D+(1-(1-D)^(gam+1))^alp*W/WF;

    tic;

while D(n)<0.9
    stress11=load*sin((n)*2*pi/stepnumber);
    m=1/3*sum(stress11+0+0);
    dev1=[stress11 0 0 ;0 0 0 ;0 0 0 ]-m*diag([1,1,1]);
    dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
    dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
    dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
    
    stress11=load*sin((n+1)*2*pi/stepnumber);
    m=1/3*sum(stress11+0+0);
    devn=[stress11 0 0;0 0 0;0 0 0]-m*diag([1,1,1]);
    dev11g=devn(1,1); dev12g=devn(1,2); dev13g=devn(1,3);
    dev21g=devn(2,1); dev22g=devn(2,2); dev23g=devn(2,3);
    dev31g=devn(3,1); dev32g=devn(3,2); dev33g=devn(3,3);
    
    trial11=bsxfun(@plus,Sb11,(dev11g-dev11)); trial12=bsxfun(@plus,Sb12,(dev12g-dev12));trial13=bsxfun(@plus,Sb13,(dev13g-dev13));
    trial21=bsxfun(@plus,Sb21,(dev21g-dev21)); trial22=bsxfun(@plus,Sb22,(dev22g-dev22));trial23=bsxfun(@plus,Sb23,(dev23g-dev23));
    trial31=bsxfun(@plus,Sb31,(dev31g-dev31)); trial32=bsxfun(@plus,Sb32,(dev32g-dev32));trial33=bsxfun(@plus,Sb33,(dev33g-dev33));
    trialtensor=[trial11; trial12; trial13; trial21; trial22; trial23;trial31; trial32; trial33];
    normtrial=sqrt(sum(trialtensor.^2));
yield=y*(1-D(n))^delta;% yield function evolve with damage D
%     yield=y;
    eta=bsxfun(@minus,bsxfun(@times,normtrial/yield,s),1); %1*25
    eta(eta<0)=0;
    
    Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
    Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
    Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
    %1*25 for each Sb element
    Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
    normSb=sqrt(sum((Sbtensor.^2)));
    
    Ws=(bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,[s]))<=0).*...
        (0)+...
        (bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,[s]))>0).*...
        ((E-k)*(1+nu)/(2*E*(E+k*nu))*bsxfun(@times,[weight],bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,[s])),yield),[s])));
    
    W= sum(Ws);
    D(n+1)=D(n)+(1-(1-D(n))^(gam+1))^alp*W/WF;
    
%     hold on;
%     yield1=plot (n,yield*s(45).^-1, 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 10, ...
%         'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'c');
%     Trial1=plot (n,sign(trial11(45))*normtrial(45),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 10, ...
%         'MarkerEdgeColor','r', 'MarkerFaceColor','r');
%     Sb1=plot (n,sign(Sb11(45))*normSb(45),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 10, ...
%         'MarkerEdgeColor','g', 'MarkerFaceColor','g');
%     yield8=plot (n,yield*s(61).^-1,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 10, ...
%         'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b');
%     Trial8=plot (n,sign(trial11(61))*normtrial(61),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 10, ...
%         'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor',[1 0.5 0]);
%     Sb8=plot (n,sign(Sb11(61))*normSb(61),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 10, ...
%         'MarkerEdgeColor','k', 'MarkerFaceColor','k');
   
    n=n+1;

end
toc;
n  
    % t=n/stepnumber*1/f;
    % disp(['Time to failure is ' num2str(t) ' s.']);
 hold on;
  DamageN=plot ((1:n),D(1:n),'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, ...
   'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'm');
  yield=plot ((1:n),y*(1-D(1:n)).^delta,'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, ...
   'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'm');
% %---------------------Difference between cyclic load calculation and numerical method as function of time-----------------------------
%      hold on
%      Damagediff=plot ((Dcyc(1:n-600)-D(1:n-600)).*Dcyc(1:n-600).^-1,'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, ...
%        'MarkerEdgeColor',  'k', 'MarkerFaceColor' , 'k');
% grid on;
% grid minor;
% set(gca ,'FontSize',25);
% hTitle = title('Relative difference between cyclic load calculation and numerical method' ,'Fontsize' ,30);
% hXLabel = xlabel('Number of steps' ,'Fontsize' ,30);
% hYLabel = ylabel('Relative difference', 'Fontsize' ,30);
% % Adjust font
% set(gca, 'FontName', 'Helvetica')
% set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')
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
% set(gcf, 'PaperPosition', [0 0 1920 1080]); %set(gcf,'PaperPosition',[left,bottom,width,height])
%  saveas(gcf,'F:\Git\PhDreport\4Anew\figures\Damagediff.png');

grid on;
grid minor;
hTitle = title('Damage evolution comparison of three methods' ,'Fontsize' ,30);
hXLabel = xlabel('Number of steps' ,'Fontsize' ,30);
hYLabel = ylabel('Damage', 'Fontsize' ,30);
 hLegend=legend([DamageN,DamageCha,Damagecyc],'Numerical method','Chaboche method',...
   'Cyclic load calculation');
set([hLegend, gca], 'FontSize', 20)
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
set(gcf, 'PaperPosition', [0 0 1920 1080]); %set(gcf,'PaperPosition',[left,bottom,width,height])
% saveas(gcf,'F:\Git\PhDreport\4Anew\figures\damagesin.png');





% %---------------------Plot Trial and Sb evolution-----------------------------
% figure(1);
%  hold on;
%   Trial1=plot ((1:n),normtrial(1:n,1),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize',12, ...
%     'MarkerEdgeColor','r', 'MarkerFaceColor','none');
%   Trial61=plot ((1:n),normtrial(1:n,61),'LineStyle', 'none','LineWidth', 1,'Marker', 's', 'MarkerSize',12, ...
%     'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor','none');
%    Sb1=plot ((1:n),normSb(1:n,1),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 8, ...
%     'MarkerEdgeColor','none', 'MarkerFaceColor',[96 96 96]/255);
%    Sb61=plot ((1:n),normSb(1:n,61),'LineStyle', 'none','LineWidth', 1,'Marker', 's', 'MarkerSize',8, ...
%     'MarkerEdgeColor','none', 'MarkerFaceColor','g');
%   yield1=plot ((1:n),yield(1:n)*s(1).^-1,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 8, ...
%     'MarkerEdgeColor', 'none', 'MarkerFaceColor','b');
%   yield61=plot ((1:n),yield(1:n)*s(61).^-1, 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
%     'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'c');
% %---------------------plot settings-----------------------------
% grid on;
% grid minor;
% set(gca ,'FontSize',20);
% hXLabel = xlabel('t(s)' ,'Fontsize' ,20);
% hYLabel = ylabel('Stress(Pa)', 'Fontsize' ,20);
% hTitle = title('Microscopic stress evolution at 2 scales' ,'Fontsize' ,20);
% set(hTitle, 'FontSize', 20, 'FontWeight' , 'bold')
% hLegend=legend([yield1,Sb1,Trial1,yield8,Sb8,Trial8],'(\sigma_y-\lambda\Sigma_H)/s_1     at scale s_1','||S-b||              at scale s_1',...
%     '||S-b||_{trial}         at scale s_1', '(\sigma_y-\lambda\Sigma_H)/s_8     at scale s_{8}','||S-b||              at scale s_{8}','||S-b||_{trial}         at scale s_{8}');
% set([hLegend, gca], 'FontSize', 20)
% % Adjust font
% set(gca, 'FontName', 'Helvetica')
% set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')
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
% set(gcf, 'PaperPosition', [0 0 1920 1080]); %set(gcf,'PaperPosition',[left,bottom,width,height])
% % saveas(gcf,'trialsin with D.png');


%sp=actxserver('SAPI.SpVoice');
% sp.Speak('Fuck that I finished all this shit finally');
%mail2me('job finished',['Elapsed time is ' num2str(toc) ' seconds. Real test time is ' testtime ' seconds. Number of points to failure is ' NF ' points.']);