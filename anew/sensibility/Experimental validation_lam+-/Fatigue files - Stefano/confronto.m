%%
load Jabbado8.mat
Jabbado_pp = data_preprocessing8(Jabbado,f_1,t_1);
Jabbado84_pp = data_preprocessing8(Jabbado84,f84_1,t84_1);

figure;
plot(Jabbado_pp(:,3),Jabbado_pp(:,1)-Jabbado_pp(:,2), '*');
hold on
plot(Jabbado_pp(1:15,3),Jabbado_pp(1:15,1)-Jabbado_pp(1:15,2), '*r')

figure;
plot(Jabbado84_pp(:,3),Jabbado84_pp(:,1)-Jabbado84_pp(:,2), '*')
hold on
plot(Jabbado84_pp(1:8,3),Jabbado84_pp(1:8,1)-Jabbado84_pp(1:8,2), 'r*')