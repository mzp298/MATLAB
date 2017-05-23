clear;clc;
close all
format long e

% fid = fopen('F:\Git\Cetim\ep_b_11\Acqui_CV.txt');
fid = fopen('/home/ma/MATLAB/Cetim/ep_b_11/Acqui_CV.txt');
[force]=textscan(fid,'%*s%*s%s%*s%*s','headerlines',5);
area=2.98*10.01*1e-6; %meter square
stress11=1000*str2double(strrep(force{1,1},',','.')).*area^-1; %Pa
repetition=xlsread('alpha_varying.xlsx',1,'G12');
nbrep=xlsread('alpha_varying.xlsx',1,'F12');
copy=nbrep+700;
stress11=repmat(stress11(1:repetition),copy,1);
clear force;

%% 
%---------------------Verified parameters in random loading case-----------------------------
b=1.5;                    %weakening scales distribution exponent
a=0.1;                %sensitivity of sequence effect(control alp>0)
WF=5.58e8;            %dissipated energy to failure per unit volume
X=1.065;                 %Major damage threshold (meaning S_{max} must be greater than X/2*yield to activate magnification)
pb=1.4;                %Major damage power
%---------------------Verified parameters in constant loading case-----------------------------
y=230e6;           %macroscopic yield stress
E=72e9;              %Young's modulus
k=6e8;                 %hardening parameter
nu=0.3;                     %poisson's ratio
sigu=320e6;             %ultimite stress
n0=1;                   %number of initial local defects
lam=0.3;               %hydrostatic pressure sensitivity
gam=0.1;             %material parameter from Chaboche law(Wohler curve exponent)


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


%---------------------Vecterization-----------------------------
D=1e-16;                    %initial damage
n=1;                      %initial recording point

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

sequence=((Smax*yield^-1)*(X-Smax*yield^-1)^-1)^pb;
sequence(sequence<0)=0;
alp=1-a*sequence;
D=D+(1-(1-D)^(gam+1))^alp*(1-D)^-gam*W/WF;


tic;
while D<1
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
     
    
    sequence=((Smax*yield^-1)*(X-Smax*yield^-1)^-1)^pb;
    sequence(sequence<0)=0;
    alp=1-a*sequence;
    D=D+(1-(1-D)^(gam+1))^alp*(1-D)^-gam*W/WF;
    
    n=n+1;
end
toc;
disp(['Number of test points is ' num2str(n) ' points.']);
NF=4684211;
disp(['Error is ' sprintf('%2.2f%%', ((NF-n)/NF)*100) ' .']);

%% Initialisation of POI Libs
% Add Java POI Libs to matlab javapath
javaaddpath('poi_library/poi-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
javaaddpath('poi_library/xmlbeans-2.3.0.jar');
javaaddpath('poi_library/dom4j-1.6.1.jar');
javaaddpath('poi_library/stax-api-1.0.1.jar');

xlwrite('alpha_varying.xlsx',max(stress11),1,'C12');
xlwrite('alpha_varying.xlsx',n,1,'D12');
p70=length(find(stress11(1:repetition)>70e6))/repetition;
xlwrite('/home/ma/b_b1.5_X1.065_a0.1_pb1.4/greater_stress_proportion.xlsx',p70,1,'B13');
p90=length(find(stress11(1:repetition)>90e6))/repetition;
xlwrite('/home/ma/b_b1.5_X1.065_a0.1_pb1.4/greater_stress_proportion.xlsx',p90,1,'C13');
p110=length(find(stress11(1:repetition)>110e6))/repetition;
xlwrite('/home/ma/b_b1.5_X1.065_a0.1_pb1.4/greater_stress_proportion.xlsx',p110,1,'D13');
p130=length(find(stress11(1:repetition)>130e6))/repetition;
xlwrite('/home/ma/b_b1.5_X1.065_a0.1_pb1.4/greater_stress_proportion.xlsx',p130,1,'E13');
p150=length(find(stress11(1:repetition)>150e6))/repetition;
xlwrite('/home/ma/b_b1.5_X1.065_a0.1_pb1.4/greater_stress_proportion.xlsx',p150,1,'F13');
p170=length(find(stress11(1:repetition)>170e6))/repetition;
xlwrite('/home/ma/b_b1.5_X1.065_a0.1_pb1.4/greater_stress_proportion.xlsx',p170,1,'G13');
p190=length(find(stress11(1:repetition)>190e6))/repetition;
xlwrite('/home/ma/b_b1.5_X1.065_a0.1_pb1.4/greater_stress_proportion.xlsx',p190,1,'H13');