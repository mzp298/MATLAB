%-----------------------2 analytical original(without logform)------------------------
fun_analytical=@(parameters)parameters(1).*E.*(E+k.*nu).*parameters(2).*(parameters(2)+1).*...
    (2.*(1-[alp_ben,alp_tor]).*(E-k).*(1+nu).*(parameters(2)-1)).^-1.*...
    ([Smax_ben, Smax_tor].^(parameters(2)+1).*(y-lamplus.*[hydroplus'; zeros(size(stresstor))]').^(1-parameters(2))...
    +[Smax_ben, Smax_tor].^(parameters(2)+1).*(y-lamminus.*[hydrominus'; zeros(size(stresstor))]').^(1-parameters(2))).^-1-[NFben; NFtor]';
% 
% lb  =  [0,     1];
% ub =  [1e12,   100];
% % %------initial fitting parameters
% parameters0=[8e8,3];
% W0(1)=parameters0(1);
% b_test(1)=parameters0(2);
% b=b_test(1);
% run('AL6082T6_bt1D_new_alp.m'); %-to get new alp(1) ---
% [parameters_fit,resnorm,exitflag]=lsqnonlin(fun_analytical,parameters0,lb,ub);
% W0(2)=parameters_fit(1); %original fit
% b_test(2)=parameters_fit(2); %original fit
% %------iterate to b converge-----------------
% p=2;
% while abs(b_test(p)-b_test(p-1))/b_test(p-1)>1e-4 %converge threshold
% b=b_test(p);
% save('AL6082T6.mat','W0','b','-append');
% run('AL6082T6_bt1D_new_alp.m');  %--use 'W0','b' to get new alp, NF---
% parameters0=[W0(p),b_test(p)];  %from pth fit
% [parameters_fit,resnorm,exitflag]=lsqnonlin(fun_analytical,parameters0,lb,ub);
% p=p+1;
% W0(p)=parameters_fit(1);  %to (p+1)th fit
% b_test(p)=parameters_fit(2)
% if abs(b_test(p)-b_test(2))/b_test(2)>10 %if p+1 fit is too far from original fit
%     b_test(p)=b_test(2);
%     W0(p)=W0(2);
% else
%     b_test(p)=b_test(p);
%     W0(p)=W0(p);
% end
% end
% 
% W0=W0(p);
% b=b_test(p);
% W0
% b
% lamplus
% lamminus

parameters=[W0,b];
NF_num=fun_analytical(parameters)+[NFben; NFtor]';
NFben_num=NF_num(1:length(NFben));
NFtor_num=NF_num(1+length(NFben):length(NFtor)+length(NFben));
save('AL6082T6.mat','W0','b','NFben_num','NFtor_num','alp_ben','alp_tor','-append'); %save the best fit W0, b, alp_ben, alp_tor