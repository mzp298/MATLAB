clear;clc;
close all
format long e
run('Parameters_ep_a_b')
load('gaussian.mat');
%vpa(x,289);
%vpa(weight,289);

%%----------epa_01-----
% fid = fopen('F:\Git\Cetim\ep_a_02\Acqui_CV.txt');
fid = fopen('/home/ma/MATLAB/Cetim/ep_a_02/Acqui_CV.txt');
area=2.29*9.99*1e-6; %meter square
NF=414298;
% run('Referece_scalar_iteration.m');

%---------------------1 Numerical method with optimal time steps-----------------------------
[force]=textscan(fid,'%*s%*s%s%*s%*s',repetition,'headerlines',5);
fclose(fid);
stress=1000*str2double(strrep(force{1,1},',','.')).*area^-1; %Pa
max(stress)
stress=stress(991:1000)
for i=2:length(stress)
    stress11(1+ari*(i-2):1+ari*(i-1))=linspace(stress(i-1),stress(i),ari+1);
end

n=1;
tensor = [stress11(1) 0 0 ;...
    0 0 0 ;...
    0 0 0 ];
run('Damiter1')
g=1; %first reference alp, W, n index when generate(filtering small alpha variation)
while n<length(stress11)
    tensor = [stress11(n), 0, 0 ;...
        0, 0, 0 ;...
        0, 0, 0 ; ];
    hydro=1/3*trace(tensor);
    dev1=tensor-hydro*eye(3);
    dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
    dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
    dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
    tensor = [stress11(n+1), 0, 0 ;...
        0, 0,  0 ;...
        0, 0,  0 ; ];
    run('Damiter2')
    n=n+1;
end

%----------scalar iteration
e=1; %first D index
j=1; %first reference alp, W, n index when iterate(after adaptation)
D= 0;
while D<1 %-----------the optimal time steps can be iterated with scalar
    D= (D^(1-alp_ref(j))+(1-alp_ref(j)).*W_ref(j)/W0).^(1/(1-alp_ref(j))); % implicit
    j=j+1;
    if j+1>=length(alp_ref)
        j=1;
    end
    e=e+1;
end
N_ref_per_cycle=length(alp_ref)/(length(stress)/2);%1336
NF_num=(e-1)/N_ref_per_cycle;%5.856e+03
nF=NF_num*2
xlwrite('a_fitting.xls',nF,1,'D03');
xlwrite('a_fitting.xls',(NF-nF)/NF,1,'E03');