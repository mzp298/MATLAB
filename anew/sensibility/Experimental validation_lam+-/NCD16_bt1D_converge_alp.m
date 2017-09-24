%%----------------------converge beta----------------
%%-----------------------1 analytical log form(maybe singularity)------------------------
% parameters=[2e9, 2] %test if the function is one colume(row)
% fun_analytical=@(parameters)log(parameters(1).*(4.*(E-k).*(1+nu)).^-1.*...% approximately 1e-1(to be tuned)
%     (1-[alp_ben,alp_tor]).^-1.*E.*(E+k.*nu).*...% approximately 1e22(fix)
%     parameters(2).*(parameters(2)+1)*(parameters(2)-1)^-1)...
%     +(parameters(2)-1).*log(y-1/3.*lam.*[stressben+m; zeros(size(stresstor))]')...
%     -(parameters(2)+1).*log([Smax_ben, Smax_tor]) ...% approximately 1e-15(to be tuned)
%     -log([NFben; NFtor]');
%-----------------------2 analytical original(without logform)------------------------
fun_analytical=@(parameters)parameters(1).*E.*(E+k.*nu).*parameters(2).*(parameters(2)+1).*...
    (2.*(1-[alp_ben,alp_tor]).*(E-k).*(1+nu).*(parameters(2)-1)).^-1.*...
    ([Smax_ben, Smax_tor].^(parameters(2)+1).*(y-1/3.*lamplus.*[hydroplus'; zeros(size(stresstor))]').^(1-parameters(2))...
    +[Smax_ben, Smax_tor].^(parameters(2)+1).*(y-1/3.*lamminus.*[hydrominus'; zeros(size(stresstor))]').^(1-parameters(2))).^-1-[NFben; NFtor]';

lb  =  [0,     1];
ub =  [1e12,   100];
% %------initial fitting parameters
parameters0=[1e9,3];
W0(1)=parameters0(1);
b_test(1)=parameters0(2);
b=b_test(1);
run('NCD16_bt1D_new_alp.m'); %-to get new alp(1) ---
[parameters_fit,resnorm,exitflag]=lsqnonlin(fun_analytical,parameters0,lb,ub);
%%--------------or we use specific options for lsqnonlin----------------
% options = optimoptions(@ lsqnonlin,'Algorithm','levenberg-marquardt',...
% options = optimoptions(@ lsqnonlin,'Algorithm','trust-region-reflective',...
%     'StepTolerance',1e-30,...
%     'MaxIterations',10000,...
%     'MaxFunctionEvaluations',10000,...
%     'Display','iter');
% [parameters_fit,resnorm]=lsqnonlin(fun_analytical,parameters0,lb,ub,options);
W0(2)=parameters_fit(1); %original fit
b_test(2)=parameters_fit(2); %original fit
%------iterate to b converge-----------------
p=2;
while abs(b_test(p)-b_test(p-1))/b_test(p-1)>1e-4
b=b_test(p);
save('NCD16.mat','W0','b','a','E','k','nu','y','stepnumber','cycles','delta_alp');
run('NCD16_bt1D_new_alp.m'); %--to get new alp(p) using b_test(p)---
parameters0=[W0(p),b_test(p)];  %from pth fit
[parameters_fit,resnorm,exitflag]=lsqnonlin(fun_analytical,parameters0,lb,ub);
% [parameters_fit,resnorm,exitflag]=lsqnonlin(fun_analytical,parameters0,lb,ub,options);
p=p+1;
W0(p)=parameters_fit(1);  %to get new fit
b_test(p)=parameters_fit(2) %to get new fit
if abs(b_test(p)-b_test(2))/b_test(2)>5 %if new fit is too far from original fit(5 times)
    b_test(p)=b_test(2);
    W0(p)=W0(2);
else
    b_test(p)=b_test(p);
    W0(p)=W0(p);
end
end

W0=W0(p);
b=b_test(p);
W0
b
lamplus
lamminus
parameters=[W0,b];
NF_num=fun_analytical(parameters)+[NFben; NFtor]';
NFben_num=NF_num(1:length(NFben));
NFtor_num=NF_num(1+length(NFben):length(NFtor)+length(NFben));
save('NCD16.mat','W0','b','NFben_num','NFtor_num','alp_ben','alp_tor','a','E','k','nu','y','stepnumber','cycles','delta_alp','cycles90'); %save the best fit W0, b, alp_ben, alp_tor