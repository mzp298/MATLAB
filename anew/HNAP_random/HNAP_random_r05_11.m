
run('Parameters_HNAP');

m=0;
sigm=0;
scentre=[2*sigm/3            0                0 ;...
                 0             -sigm/3                0 ;...
                 0                       0      -sigm/3];
             
sig_xx_original=f_sig_xx(32)*force{1,1}*1e6;%Pa
tau_xy_original=f_tau_xy(32)*force{1,1}*1e6; %Pa
for i=2:length(sig_xx_original)
    sig_xx(1+ari*(i-2):1+ari*(i-1))=linspace(sig_xx_original(i-1),sig_xx_original(i),ari+1);
    tau_xy(1+ari*(i-2):1+ari*(i-1))=linspace(tau_xy_original(i-1),tau_xy_original(i),ari+1);
end
n=1;
tensor = [sig_xx(1) tau_xy(1) 0 ;...
    tau_xy(1) 0 0 ;...
    0 0 0 ];
run('Damiter1')
g=1; %first reference alp, W, n index when generate(filtering small alpha variation)
while n<length(sig_xx)
    tensor = [sig_xx(n), tau_xy(n), 0 ;...
        tau_xy(n), 0, 0 ;...
        0, 0, 0 ; ];
    hydro(n)=1/3*trace(tensor);
    dev1=tensor-hydro(n)*eye(3)-scentre;
    dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
    dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
    dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
    tensor = [sig_xx(n+1), tau_xy(n+1), 0 ;...
        tau_xy(n+1), 0,  0 ;...
        0, 0,  0 ; ];
    run('Damiter2')
    n=n+1;
end
length(alp)
min(alp)

%% 
%----------scalar iteration
e=1; %first D index
j=1; %first reference alp, W, n index when iterate(after adaptation)
D= 1e-16;
tic;
while D<1 %-----------the optimal time steps can be iterated with scalar
%     D= (D^(1-alp_ref(j))+(1-alp_ref(j)).*W_ref(j)/W0).^(1/(1-alp_ref(j))); % implicit
    D=D+D^alp_ref(j)*W_ref(j)/W0;
    j=j+1;
    if j+1>=length(alp_ref)
        j=1;
    end
    e=e+1;
end
toc;

nF_num=(e-1)/ari;
nF_exp=T_exp(32)/184.32*repetition;
error=(nF_exp-nF_num)/nF_exp
xlwrite('HNAP_random.xls',nF_num,1,'G33');
xlwrite('HNAP_random.xls',error,1,'H33');

