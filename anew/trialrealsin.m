% Program to get the Gauss-Legendre Quadrature results
clear;clc;

dbstop if error
format long e
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

E=2e11;               %Young¡¯s modulus
k=6e8;                  %hardening parameter
b=3;                      %weakening scales distribution exponent
nu=0.3;                 %poisson's ratio
tt=2e8;                  %torsion fatigue limit
ff=2.5e8;              %bending fatigue limit
ac=(tt-ff/sqrt(3))/(ff/3); %crossland criterial constant
bc=tt;                    %crossland criterial constant 
sigu=8e8;             %ultimite stress
gam=0.5;              %material parameter from Chaboche law(Wohler curve exponent)
y=6.38e8;            %macroscopic yield stress
WF=3e6;             %dissipated energy to failure per unit volume
load=5e8;            %cyclic load
loadtensor= [load 0 0;0 0 0;0 0 0];
stepnumber=200;        %devide one cycle in 200 parts
f=50;                            %frequency of load

%---------------------numerical method-----------------------------
alp=0.5;
D=0;             %initial damage
n=1;       %initial recording point
G = (1 - (1 - D).^(gam + 1)).^(1-alp);
    stress11(1)=load*sin(2*pi/stepnumber);
    m(1)=1/3*sum(stress11(1)+0+0);
    dev1=[stress11(1) 0 0 ;0 0 0 ;0 0 0 ]-m(1)*diag([1,1,1]);
    dev11(1)=dev1(1,1); dev12(1)=dev1(1,2); dev13(1)=dev1(1,3);
    dev21(1)=dev1(2,1); dev22(1)=dev1(2,2); dev23(1)=dev1(2,3);
    dev31(1)=dev1(3,1); dev32(1)=dev1(3,2); dev33(1)=dev1(3,3);
    for i=1:25 %creat a set of Ws[n*25 vector]
    s(i)= (x(i)/2+1/2).^(1/(1-b));

    trial11(1,i)=dev11(1); trial12(1,i)=dev12(1); trial13(1,i)=dev13(1);
    trial21(1,i)=dev21(1); trial22(1,i)=dev22(1); trial23(1,i)=dev23(1);
    trial31(1,i)=dev31(1); trial32(1,i)=dev32(1); trial33(1,i)=dev33(1);
     
    normtrial(1,i)=norm([trial11(1,i), trial12(1,i), trial13(1,i); trial21(1,i), trial22(1,i), trial23(1,i);trial31(1,i), trial32(1,i), trial33(1,i)],'fro');
    eta(1,i)=normtrial(1,i)/y*s(i)-1;
    eta(eta<0)=0;
    
    Sb11(1,i)=trial11(1,i)*(1+eta(1,i)).^-1;Sb12(1,i)=trial12(1,i)*(1+eta(1,i)).^-1;Sb13(1,i)=trial13(1,i)*(1+eta(1,i)).^-1;
    Sb21(1,i)=trial21(1,i)*(1+eta(1,i)).^-1;Sb22(1,i)=trial22(1,i)*(1+eta(1,i)).^-1;Sb23(1,i)=trial23(1,i)*(1+eta(1,i)).^-1;
    Sb31(1,i)=trial31(1,i)*(1+eta(1,i)).^-1;Sb32(1,i)=trial32(1,i)*(1+eta(1,i)).^-1;Sb33(1,i)=trial33(1,i)*(1+eta(1,i)).^-1;
    normSb(n,i)=norm([Sb11(n,i), Sb12(n,i), Sb13(n,i); Sb21(n,i), Sb22(n,i), Sb23(n,i);Sb31(n,i), Sb32(n,i), Sb33(n,i)],'fro');
    end
tic;    
while G<1
    stress11(n+1)=load*sin((n+1)*2*pi/stepnumber);
    m(n+1)=1/3*sum(stress11(n+1)+0+0);
    
    devn=[stress11(n+1) 0 0;0 0 0;0 0 0]-m(n+1)*diag([1,1,1]);
    dev11(n+1)=devn(1,1); dev12(n+1)=devn(1,2); dev13(n+1)=devn(1,3);
    dev21(n+1)=devn(2,1); dev22(n+1)=devn(2,2); dev23(n+1)=devn(2,3);
    dev31(n+1)=devn(3,1); dev32(n+1)=devn(3,2); dev33(n+1)=devn(3,3);
    for i=1:25 %creat a set of Ws[n*25 vector]
    
    trial11(n+1,i)=Sb11(n,i)+(dev11(n+1)-dev11(n)); trial12(n+1,i)=Sb12(n,i)+(dev12(n+1)-dev12(n)); trial13(n+1,i)=Sb13(n,i)+(dev13(n+1)-dev13(n));
    trial21(n+1,i)=Sb21(n,i)+(dev21(n+1)-dev21(n)); trial22(n+1,i)=Sb22(n,i)+(dev22(n+1)-dev22(n)); trial23(n+1,i)=Sb23(n,i)+(dev23(n+1)-dev23(n));
    trial31(n+1,i)=Sb31(n,i)+(dev31(n+1)-dev31(n)); trial32(n+1,i)=Sb32(n,i)+(dev32(n+1)-dev32(n)); trial33(n+1,i)=Sb33(n,i)+(dev33(n+1)-dev33(n));
    
    normtrial(n+1,i)=norm([trial11(n+1,i), trial12(n+1,i), trial13(n+1,i); trial21(n+1,i), trial22(n+1,i), trial23(n+1,i);trial31(n+1,i), trial32(n+1,i), trial33(n+1,i)],'fro');
    eta(n+1,i)=normtrial(n+1,i)/y*s(i)-1;
    eta(eta<0)=0;
    
    Sb11(n+1,i)=trial11(n+1,i)*(1+eta(n+1,i)).^-1;Sb12(n+1,i)=trial12(n+1,i)*(1+eta(n+1,i)).^-1;Sb13(n+1,i)=trial13(n+1,i)*(1+eta(n+1,i)).^-1;
    Sb21(n+1,i)=trial21(n+1,i)*(1+eta(n+1,i)).^-1;Sb22(n+1,i)=trial22(n+1,i)*(1+eta(n+1,i)).^-1;Sb23(n+1,i)=trial23(n+1,i)*(1+eta(n+1,i)).^-1;
    Sb31(n+1,i)=trial31(n+1,i)*(1+eta(n+1,i)).^-1;Sb32(n+1,i)=trial32(n+1,i)*(1+eta(n+1,i)).^-1;Sb33(n+1,i)=trial33(n+1,i)*(1+eta(n+1,i)).^-1;
    
    normSb(n+1,i)=norm([Sb11(n+1,i), Sb12(n+1,i), Sb13(n+1,i); Sb21(n+1,i), Sb22(n+1,i), Sb23(n+1,i);Sb31(n+1,i), Sb32(n+1,i), Sb33(n+1,i)],'fro');
         if normtrial(n+1,i)<= y*s(i).^-1
             Ws(n+1,i)=0;
         else
             Ws(n+1,i)=(E-k)*(1+nu)/(2*E*(E+k*nu))*weight(i)*(normtrial(n+1,i)-y*s(i).^-1)*y*s(i).^-1;  
         end
    end
  W(1)=0;
  W(n+1)= sum(Ws(n+1,:)); %integrate all scales in step n 
  G = G+W(n+1)/WF;
  D=1-(1-G.^(1/(1-alp))).^(1/(gam + 1));
  t=n/stepnumber*1/f;
%  hold on;
%    yield1=plot (n,y*s(1).^-1, 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, ...
%     'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'c');
%    Trial1=plot (n,sign(trial11(n,1))*normtrial(n,1),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 6, ...
%     'MarkerEdgeColor','r', 'MarkerFaceColor','r');
%     Sb1=plot (n,sign(Sb11(n,1))*normSb(n,1),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 6, ...
%     'MarkerEdgeColor','g', 'MarkerFaceColor','g');
%    yield8=plot (n,y*s(8).^-1,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 6, ...
%      'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b');
%    Trial8=plot (n,sign(trial11(n,8))*normtrial(n,8),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 6, ...
%     'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor',[1 0.5 0]);
%    Sb8=plot (n,sign(Sb11(n,8))*normSb(n,8),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 6, ...
%     'MarkerEdgeColor','k', 'MarkerFaceColor','k');

%   DamageN=plot (t,D,'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, ...
%    'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'm');

%---------------------Difference between cyclic load calculation and numerical method as function of time-----------------------------
%  Gcyc = Gcyc+Wcyc/stepnumber/WF
%  Dcyc=1-(1-Gcyc.^(1/(1-alp))).^(1/(gam + 1));
%  hold on
%  Damagecyc=plot (t,D-Dcyc,'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, ...
%    'MarkerEdgeColor',  'k', 'MarkerFaceColor' , 'k');
  n=n+1;
end;
toc;
disp=(['Number of points to failure is ' num2str(n) ' points.'])
% plot(normtrial)
% %---------------------chaboche method-----------------------------
% alp=0.5;
% D=0;             %initial damage
% n=1;       %initial recording point
%   G = (1 - (1 - D).^(gam + 1)).^(1-alp);
%   m=1/3*sum(diag(loadtensor));
%   S1=loadtensor-m*diag([1,1,1]);
%   sqrj1=1/2*sqrt(1/2)*norm(S1,'fro');
%   M=ff^1.233*(1-3*m/sigu);
%   while G<1
%   NF=1/((gam+1)*(1-alp))*(sqrj1/M)^(-gam);
%   G = G+1/stepnumber/NF
%   D=1-(1-G.^(1/(1-alp))).^(1/(gam + 1));
%   t=n/stepnumber*1/f;
%   hold on;
%   DamageC=plot (t,D, 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, ...
%    'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'g');
%   n=n+1;
%   end
% %---------------------Cyclic load calculation-----------------------------
% Dcyc=0;
% n=1;
% Gcyc = (1 - (1 - Dcyc).^(gam + 1)).^(1-alp);
% Wcyc=4*(E-k)*(1+nu)*(b-1)/(E*(E+k*nu)*b*(b+1))*norm(loadtensor-(1/3*sum(diag(loadtensor)))*diag([1,1,1]),'fro').^(b+1)*y.^(1-b) ;
% while Gcyc< 1
%  Gcyc = Gcyc+Wcyc/stepnumber/WF
%  Dcyc=1-(1-Gcyc.^(1/(1-alp))).^(1/(gam + 1));
%  t=n/stepnumber*1/f;
%  hold on
%  Damagecyc=plot (t,Dcyc,'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, ...
%    'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'b');
%  n=n+1;
% end

%---------------------plot settings-----------------------------
grid on;
grid minor;
% axis([0 0.49 -0.04 0.04]);
set(gca ,'FontSize',30);
hXLabel = xlabel('t(s)' ,'Fontsize' ,30);

 hTitle =title('Microscopic stress evolution at 2 different scales ' ,'Fontsize' ,30);
 hYLabel =ylabel('Stress(Pa)', 'Fontsize' ,30);
hLegend=legend([yield1,Sb1,Trial1,yield8,Sb8,Trial8],'(\sigma_y-\lambda\Sigma_H)/s_1     at scale s_1','(S-b)               at scale s_1',...
   '(S-b)_{trial}          at scale s_1', '(\sigma_y-\lambda\Sigma_H)/s_8     at scale s_{8}','(S-b)               at scale s_{8}','(S-b)_{trial}          at scale s_{8}');


% hTitle = title('Damage evolution comparison of three methods' ,'Fontsize' ,30);
% hYLabel = ylabel('Damage', 'Fontsize' ,30);
%  hLegend=legend([DamageN,DamageC,Damagecyc],'Numerical method','Chaboche method',...
%    'Cyclic load calculation');

% hTitle = title('Difference between cyclic load calculation and numerical method as function of time(time step=1/15000)' ,'Fontsize' ,30);

% Adjust font
set(gca, 'FontName', 'Helvetica')
set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')

set([hXLabel, hYLabel], 'FontSize', 30)
set(hTitle, 'FontSize', 30, 'FontWeight' , 'bold')
set([hLegend, gca], 'FontSize', 30)
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
saveas(gcf,'trialsin.png');
% saveas(gcf,'damagesin.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('Fuck that I finished all this shit finally');

mail2me('job finished',['Elapsed time is ' num2str(toc) ' seconds.']);
