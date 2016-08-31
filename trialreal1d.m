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

forceorigin=repmat(forcex,3,1);

ari=10; %  Arithmetic sequence between every recorded points
for i=2:(1*802805) 
force(1+ari*(i-2):1+ari*(i-1))=linspace(forceorigin(i-1),forceorigin(i),ari+1);
end;
%  ari*(i-1)+1 %the number of points
A=1/1e6;
stress=1/A*force;

[x]= [-0.99555697 -0.976663921 -0.942974571 -0.894991998 -0.833442629 -0.759259263 -0.673566368...
-0.57766293 -0.473002731 -0.361172306 -0.243866884 -0.122864693 0 0.122864693 0.243866884 0.361172306...
0.473002731 0.57766293 0.673566368 0.759259263 0.833442629 0.894991998 0.942974571 0.976663921...
0.99555697];
[weight]=[0.011393799	0.026354987	0.040939157	0.054904696	0.068038334	0.0801407	0.091028262...
    0.100535949	0.108519624	0.114858259	0.119455764	0.122242443	0.123176054	0.122242443	0.119455764...
    0.114858259	0.108519624	0.100535949	0.091028262	0.0801407	0.068038334	0.054904696	0.040939157...
    0.026354987	0.011393799];
% [x]=xlsread('Gauss-Legendre Quadrature','Sheet1','b1:z1');
% [weight]=xlsread('Gauss-Legendre Quadrature','Sheet1','b2:z2');


y=6.38e8;           %macroscopic yield stress
lam=0.5;             %hydrostatic pressure sensitivity
E=2e11;              %Young��s modulus
k=6e8;                 %hardening parameter
b=3;                    %weakening scales distribution exponent
nu=0.3;                     %poisson's ratio
tt=2e8;                 %torsion fatigue limit
ff=2.5e8;              %bending fatigue limit
ac=(tt-ff/sqrt(3))/(ff/3); %crossland criterial constant
bc=tt;                     %crossland criterial constant 
sigu=8e8;             %ultimite stress
gam=0.5;                %material parameter from Chaboche law(Wohler curve exponent)
samplerate=256;   %recorded samples per second

%---------------------Main-----------------------------
tic;   
D=0;             %initial damage
n=1;                      %initial recording point
WF=3e6;             %dissipated energy to failure per unit volume
alp=0.8;
step=1/samplerate/ari;
G = (1 - (1 - D).^(gam + 1)).^(1-alp);
t=n*step;
   m(1)=1/3*sum(stress(1)+0+0);
   yield(1)=y-lam*m(1); %macro yield strength considering mean stress effect
   yield(yield<0)=0;
   dev1=[stress(1) 0 0 ;0 0 0 ;0 0 0 ]-m(1)*diag([1,1,1]);
   dev11(1)=dev1(1,1); dev12(1)=dev1(1,2); dev13(1)=dev1(1,3);
   dev21(1)=dev1(2,1); dev22(1)=dev1(2,2); dev23(1)=dev1(2,3);
   dev31(1)=dev1(3,1); dev32(1)=dev1(3,2); dev33(1)=dev1(3,3);
    for i=1:25 %creat a set of Ws[n*25 vector]
    s(i)= (x(i)/2+1/2).^(1/(1-b));

    trial11(1,i)=dev11(1); trial12(1,i)=dev12(1); trial13(1,i)=dev13(1);
    trial21(1,i)=dev21(1); trial22(1,i)=dev22(1); trial23(1,i)=dev23(1);
    trial31(1,i)=dev31(1); trial32(1,i)=dev32(1); trial33(1,i)=dev33(1);
     
    normtrial(1,i)=norm([trial11(1,i), trial12(1,i), trial13(1,i); trial21(1,i), trial22(1,i), trial23(1,i);trial31(1,i), trial32(1,i), trial33(1,i)],'fro');
    eta(1,i)=normtrial(1,i)/yield(1)*s(i)-1;
    eta(eta<0)=0;
    
    Sb11(1,i)=trial11(1,i)*(1+eta(1,i)).^-1;Sb12(1,i)=trial12(1,i)*(1+eta(1,i)).^-1;Sb13(1,i)=trial13(1,i)*(1+eta(1,i)).^-1;
    Sb21(1,i)=trial21(1,i)*(1+eta(1,i)).^-1;Sb22(1,i)=trial22(1,i)*(1+eta(1,i)).^-1;Sb23(1,i)=trial23(1,i)*(1+eta(1,i)).^-1;
    Sb31(1,i)=trial31(1,i)*(1+eta(1,i)).^-1;Sb32(1,i)=trial32(1,i)*(1+eta(1,i)).^-1;Sb33(1,i)=trial33(1,i)*(1+eta(1,i)).^-1;
    normSb(n,i)=norm([Sb11(n,i), Sb12(n,i), Sb13(n,i); Sb21(n,i), Sb22(n,i), Sb23(n,i);Sb31(n,i), Sb32(n,i), Sb33(n,i)],'fro');
           if normtrial(1,i)<= yield(1)*s(i).^-1
              Ws(1,i)=0;
           else
              Ws(1,i)=(E-k)*(1+nu)/(2*E*(E+k*nu))*weight(i)*(normtrial(1,i)-yield(1)*s(i).^-1)*yield(1)*s(i).^-1;  
           end
   end
   W(1)= sum(Ws(1,:)); %integrate all scales in step 1
   G = G+W(1)/WF;
while G<0.02
   m(n+1)=1/3*sum(stress(n+1)+0+0);
   yield(n+1)=y-lam*m(n+1);
   yield(yield<0)=0;
   devn=[stress(n+1) 0 0;0 0 0;0 0 0]-m(n+1)*diag([1,1,1]);
   dev11(n+1)=devn(1,1); dev12(n+1)=devn(1,2); dev13(n+1)=devn(1,3);
   dev21(n+1)=devn(2,1); dev22(n+1)=devn(2,2); dev23(n+1)=devn(2,3);
   dev31(n+1)=devn(3,1); dev32(n+1)=devn(3,2); dev33(n+1)=devn(3,3);
  for i=1:25 %creat a set of Ws[n*25 vector]
    s(i)= (x(i)/2+1/2).^(1/(1-b));
       
    trial11(n+1,i)=Sb11(n,i)+(dev11(n+1)-dev11(n)); trial12(n+1,i)=Sb12(n,i)+(dev12(n+1)-dev12(n)); trial13(n+1,i)=Sb13(n,i)+(dev13(n+1)-dev13(n));
    trial21(n+1,i)=Sb21(n,i)+(dev21(n+1)-dev21(n)); trial22(n+1,i)=Sb22(n,i)+(dev22(n+1)-dev22(n)); trial23(n+1,i)=Sb23(n,i)+(dev23(n+1)-dev23(n));
    trial31(n+1,i)=Sb31(n,i)+(dev31(n+1)-dev31(n)); trial32(n+1,i)=Sb32(n,i)+(dev32(n+1)-dev32(n)); trial33(n+1,i)=Sb33(n,i)+(dev33(n+1)-dev33(n));
    
    normtrial(n+1,i)=norm([trial11(n+1,i), trial12(n+1,i), trial13(n+1,i); trial21(n+1,i), trial22(n+1,i), trial23(n+1,i);trial31(n+1,i), trial32(n+1,i), trial33(n+1,i)],'fro');
    eta(n+1,i)=normtrial(n+1,i)/yield(n+1)*s(i)-1;
    eta(eta<0)=0;
    
    Sb11(n+1,i)=trial11(n+1,i)*(1+eta(n+1,i)).^-1;Sb12(n+1,i)=trial12(n+1,i)*(1+eta(n+1,i)).^-1;Sb13(n+1,i)=trial13(n+1,i)*(1+eta(n+1,i)).^-1;
    Sb21(n+1,i)=trial21(n+1,i)*(1+eta(n+1,i)).^-1;Sb22(n+1,i)=trial22(n+1,i)*(1+eta(n+1,i)).^-1;Sb23(n+1,i)=trial23(n+1,i)*(1+eta(n+1,i)).^-1;
    Sb31(n+1,i)=trial31(n+1,i)*(1+eta(n+1,i)).^-1;Sb32(n+1,i)=trial32(n+1,i)*(1+eta(n+1,i)).^-1;Sb33(n+1,i)=trial33(n+1,i)*(1+eta(n+1,i)).^-1;
    
    normSb(n+1,i)=norm([Sb11(n+1,i), Sb12(n+1,i), Sb13(n+1,i); Sb21(n+1,i), Sb22(n+1,i), Sb23(n+1,i);Sb31(n+1,i), Sb32(n+1,i), Sb33(n+1,i)],'fro');
           if normtrial(n+1,i)<= yield(n+1)*s(i).^-1
              Ws(n+1,i)=0;
           else
              Ws(n+1,i)=(E-k)*(1+nu)/(2*E*(E+k*nu))*weight(i)*(normtrial(n+1,i)-yield(n+1)*s(i).^-1)*yield(n+1)*s(i).^-1;  
           end
 end
  W(n+1)= sum(Ws(n+1,:)); %integrate all scales in step n
  G = G+W(n+1)/WF;
  D=1-(1-G.^(1/(1-alp))).^(1/(gam + 1));
  t=n*step;
%   hold on;
%   yield1=plot (t,yield(n)*s(1).^-1, 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, ...
%     'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'c');
%   Trial1=plot (t,normtrial(n,1),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 6, ...
%     'MarkerEdgeColor','r', 'MarkerFaceColor','r');
%   Sb1=plot (t,normSb(n,1),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 6, ...
%     'MarkerEdgeColor','g', 'MarkerFaceColor','g');
%   yield8=plot (t,yield(n)*s(8).^-1,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 6, ...
%     'MarkerEdgeColor', 'none', 'MarkerFaceColor','b');
%   Trial8=plot (t,normtrial(n,8),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 6, ...
%     'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor',[1 0.5 0]);
%   Sb8=plot (t,normSb(n,8),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 6, ...
%     'MarkerEdgeColor','k', 'MarkerFaceColor','k');

%   DamageN=plot (t,D,'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, ...
%    'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'm');
  n=n+1;
end
toc;
disp(['Number of points to failure is ' num2str(n) ' points.']);

%---------------------plot settings-----------------------------
grid on;
grid minor;
set(gca ,'FontSize',30);
hXLabel = xlabel('t(s)' ,'Fontsize' ,30);

%  hTitle =title('Damage evolution under multidimensional stress' ,'Fontsize' ,25);
%  hYLabel =ylabel('D', 'Fontsize' ,25);

hTitle = title('Microscopic stress evolution at 2 scales' ,'Fontsize' ,30);
hYLabel = ylabel('(S-b)(Pa)', 'Fontsize' ,30);
hLegend=legend([yield1,Sb1,Trial1,yield8,Sb8,Trial8],'(\sigma_y-\lambda\Sigma_H)/s_1     at scale s_1','(S-b)               at scale s_1',...
   '(S-b)_{trial}          at scale s_1', '(\sigma_y-\lambda\Sigma_H)/s_8     at scale s_{8}','(S-b)               at scale s_{8}','(S-b)_{trial}          at scale s_{8}');
set([hLegend, gca], 'FontSize', 30)

% Adjust font
set(gca, 'FontName', 'Helvetica')
set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hXLabel, hYLabel], 'FontSize', 30)
set(hTitle, 'FontSize', 30, 'FontWeight' , 'bold')

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
%  saveas(gcf,'damage1d.png');
saveas(gcf,'trialreal1d.png');


num2str(toc);
sp=actxserver('SAPI.SpVoice');
sp.Speak('Fuck that I finished all this shit finally');

mail2me('job finished',num2str(toc));