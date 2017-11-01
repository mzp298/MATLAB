% Program to get the Gauss-Legendre Quadrature results (Vectorized)
clear;clc;close all;
dbstop if error
format long e
load('gaussian.mat');


E=72e9;               %Young's modulus
nu=0.3;                 %poisson's ratio
k=6e8;                  %hardening parameter
b=1.1;                      %weakening scales distribution exponent (between 1 and 2)
y=230e6;            %macroscopic yield stress
a=0.4;
W0=3.22e4;             %dissipated energy to failure per unit volume
lam=0.1;               %hydrostatic pressure sensitivity
m=0;
load=2.25e8;            %cyclic load
loadtensor= [load 0 0;0 0 0;0 0 0];
stepnumber=400;        %devide one cycle in X parts
delta_alp=0.05;
alp_ref=zeros(1,stepnumber+1);
n_ref=zeros(1,stepnumber+1);
W_ref=zeros(1,stepnumber+1);
cycles=5;
%---------------------2 Numerical method with alp_optimal_stepnumber-----------------------------
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
normtrial(1)=norm([trial11, trial12, trial13; trial21, trial22, trial23;trial31, trial32, trial33],'fro');
s= (x/2+1/2).^(1/(1-b)); 
yield=y-lam*hydro;
eta=bsxfun(@minus,bsxfun(@times,normtrial(1)/yield,s),1); 
eta(eta<0)=0; 
Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
normSb=sqrt(sum(Sbtensor.^2)); 
Ws=(bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))<=0).*...
    (0)+...
    (bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))>0).*...
    ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s)),yield),s)));
W= sum(Ws);
W_accumulate=W;
sequence=((Smax*yield^-1)*(1-Smax*yield^-1)^-1)^b;
sequence(sequence<0)=0;
alp=1-a*sequence;
alp_ref(1)=alp;
n_ref(1)=n;
W_ref(1)=W_accumulate;

i=1; %first reference alp, W, n index
while n<cycles*stepnumber
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
    normtrial=sqrt(sum(trialtensor.^2));
    yield(n+1)=y-lam*hydro;
    eta=bsxfun(@minus,bsxfun(@times,normtrial/yield(n+1),s),1); 
    eta(eta<0)=0;
    Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
    Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
    Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
    Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
    normSb=sqrt(sum(Sbtensor.^2)); 
    Ws=(bsxfun(@minus,normtrial,bsxfun(@rdivide, yield(n+1),s))<=0).*...
        (0)+...
        (bsxfun(@minus,normtrial,bsxfun(@rdivide, yield(n+1),s))>0).*...
        ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*...
        bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial,bsxfun(@rdivide, yield(n+1),s)),yield(n+1)),s)));
    W= sum(Ws);
    W_accumulate(n+1)=W_accumulate(n)+W;
    sequence=((Smax*yield(n+1)^-1)*(1-Smax*yield(n+1)^-1)^-1)^b;
    sequence(sequence<0)=0;
    alp(n+1)=1-a*sequence;
    
    if abs(alp(n+1)-alp_ref(i))>delta_alp %----giving scalar value to iteration after the addaptation cycle(decrease time step)
        alp_ref(i+1)=alp(n+1);
        n_ref(i+1)=n+1;
        W_ref(i+1)=W_accumulate(n+1)-W_accumulate(n_ref(i));
        i=i+1;
    end
    
    n=n+1;
end
alp_ref(alp_ref==0) = [];
n_ref(n_ref==0) = [];
W_ref(W_ref==0) = [];

e=1; %first D index
j=ceil(length(alp_ref)/cycles); %first reference alp, W, n index when iterate(after adaptation)
D= 1e-16; 
while D(e)<1 %-----------the optimal time steps can be iterated with scalar
    D(e+1)=D(e)+D(e)^alp_ref(j+1)*W_ref(j+1)/W0;
    j=j+1;
    if j+1>=length(alp_ref)
        j=ceil(length(alp_ref)/cycles);
    end
    e=e+1;
end
NF_optimal=e/(length(n_ref)/cycles)
% 
%---------------------1 Numerical method with alp_stepnumber to compare with optimal step number results-----------------------------
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
normtrial(1)=norm([trial11, trial12, trial13; trial21, trial22, trial23;trial31, trial32, trial33],'fro');
s= (x/2+1/2).^(1/(1-b)); 
yield=y-lam*hydro;
eta=bsxfun(@minus,bsxfun(@times,normtrial(1)/yield,s),1); 
eta(eta<0)=0; 
Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
normSb=sqrt(sum(Sbtensor.^2)); 
Ws=(bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))<=0).*...
    (0)+...
    (bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))>0).*...
    ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s)),yield),s)));
W= sum(Ws);
W_accumulate=W;
sequence=((Smax*yield(1)^-1)*(1-Smax*yield^-1)^-1)^b;
sequence(sequence<0)=0;
alp=1-a*sequence;
D_change_alp(2)=D_change_alp(1)+D_change_alp(1)^alp*W/W0;

i=1;
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
    normtrial=sqrt(sum(trialtensor.^2));
    yield(n+1)=y-lam*hydro;
    eta=bsxfun(@minus,bsxfun(@times,normtrial/yield(n+1),s),1); 
    eta(eta<0)=0;
    Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
    Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
    Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
    Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
    normSb=sqrt(sum(Sbtensor.^2)); 
    Ws=(bsxfun(@minus,normtrial,bsxfun(@rdivide, yield(n+1),s))<=0).*...
        (0)+...
        (bsxfun(@minus,normtrial,bsxfun(@rdivide, yield(n+1),s))>0).*...
        ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*...
        bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial,bsxfun(@rdivide, yield(n+1),s)),yield(n+1)),s)));
    W= sum(Ws);
    W_accumulate(n+1)=W_accumulate(n)+W;
    sequence=((Smax*yield(n+1)^-1)*(1-Smax*yield(n+1)^-1)^-1)^b;
    sequence(sequence<0)=0;
    alp(n+1)=1-a*sequence;
    D_change_alp(n+1)=D_change_alp(n)+D_change_alp(n)^alp(n+1)*W/W0;
    n=n+1;
end
NF_stepnumber=n/stepnumber 

plot_alp_num=plot(alp(1:stepnumber),'c','Marker','s','MarkerSize',20,'MarkerFaceColor','none','LineStyle','none', 'LineWidth', 2);
hold on
n_ref(n_ref>stepnumber) = [];
alp_ref=wkeep(alp_ref,length(n_ref),'l');
plot_alp_ref=plot(n_ref,alp_ref,'r','Marker','o','MarkerSize',20,'MarkerFaceColor','r','LineStyle','none', 'LineWidth', 2);
hLegend=legend([plot_alp_num,plot_alp_ref],...
    ['Original $\alpha$ evolution in unit cycle'],...
    ['Optimal $\alpha$ evolution in unit cycle($\Delta\alpha=0.05$)'],...
    'Location','Southeast');
set(hLegend,'Interpreter','latex');
set([hLegend, gca], 'FontSize', 30)
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white
grid on;
grid minor;
hXLabel = xlabel('$N_F$');
hYLabel =ylabel('$\alpha(s_{min})$');
set(hXLabel ,'Interpreter','latex');
set(hYLabel ,'Interpreter','latex');
set(gca, 'FontName', 'Helvetica')
set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hXLabel, hYLabel], 'FontSize',40)
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1600 1000]); %set(gcf,'PaperPosition',[left,bottom,width,height])
% saveas(gcf,'F:\Git\Anew\figures\alpha_opt_vs_alpha_stepnumber.png');
sp=actxserver('SAPI.SpVoice');
sp.Speak('xuelin');

