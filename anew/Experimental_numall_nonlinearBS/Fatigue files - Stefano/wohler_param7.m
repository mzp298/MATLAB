function [time,charge, timepred,chargepred] = wohler_param7(param, Jabbado, n0)
% Wohler function for the model
% Works with param = wohlerfit7 (data, n_services_user_details)
% 
% Param is the array of parameters of the constitutive law
% Jabbado is the data set used
% N0 is the number of observations used to calibrate the model

% Parameters of the constitutive law
Gamma_0 = param(1);
mu = param(2)/2;
c = param(2)/2;
alpha = param(3);
beta = param(4);

% Parameters for the numeric method
y0 = 0;
tmax = 1e10;
tspan = [0 tmax];
options = odeset('RelTol',1e-7,'AbsTol',1e-6,'Events',@events);

% Calculation of break times at load variation
n = size(Jabbado,1);
charge = zeros(1,n);
time   = zeros(1,n);
T = Jabbado(:,1);
T_0 = Jabbado(:,2);
for i = 1:n
    % Constitutive law
    
%     f = @(t,x) 4./(mu + c -...
%         alpha/beta.*abs(x-Gamma_0).^(alpha-1).*sign(x-Gamma_0)).*...
%         ((T(i)-T_0(i))/2 + abs(x-Gamma_0).^alpha./beta);
%     fx = @(x) f(0,x);
% 
%     % Calculation of evolution in a particular case
%     [t_temp,y_temp] = ode45(f,tspan,y0,options);
%     charge(i) = T(i) - T_0(i);
%     time(i) = t_temp(length(t_temp));
    
    tmp = (2/beta/(T(i) - T_0(i)))^((alpha-1)/alpha);
    gamma = (pi*tmp + sqrt(pi*alpha^2/(alpha-1)/(mu + c)/beta*tmp*sin(pi/alpha)))/...
        (pi*tmp - alpha^2/(alpha-1)/beta/(mu + c)*sin(pi/alpha));
    time(i) = (gamma + 1)*pi/(alpha*sin(pi/alpha))*((mu + c)/4)*...
        beta^(1/alpha)/gamma*((T(i)-T_0(i))./2)^((1-alpha)/alpha) -...
        alpha/4/(alpha-1)/(gamma-1) -...
        (mu + c)*beta/4/(alpha - 1)*Gamma_0^(1-alpha);
    charge(i) = T(i) - T_0(i);
end

% Calculation of break times at load variation
n = 32;
chargepred = zeros(1,n);
timepred   = zeros(1,n);
DT = 5*T_0(1);
for i = 1:n
    % Legge costitutiva
    
%     f = @(t,x) 4./(mu + c -...
%         alpha/beta.*abs(x-Gamma_0).^(alpha-1).*sign(x-Gamma_0)).*...
%         ((Tpred-T_0)/2 + abs(x-Gamma_0).^alpha./beta);
%     fx = @(x) f(0,x);
% 
%     % Calcolo della evoluzione in un caso particolare
%     [t_temp,y_temp] = ode45(f,tspan,y0,options);
%     chargepred(i) = Tpred;
%     timepred(i) = t_temp(length(t_temp));

    tmp = (2/beta/(DT))^((alpha-1)/alpha);
    gamma = (pi*tmp + sqrt(pi*alpha^2/(alpha-1)/(mu + c)/beta*tmp*sin(pi/alpha)))/...
        (pi*tmp - alpha^2/(alpha-1)/beta/(mu + c)*sin(pi/alpha));
    timepred(i) = (gamma + 1)*pi/(alpha*sin(pi/alpha))*((mu + c)/4)*...
        beta^(1/alpha)/gamma*(DT./2)^((1-alpha)/alpha) -...
        alpha/4/(alpha-1)/(gamma-1) -...
        (mu + c)*beta/4/(alpha - 1)*Gamma_0^(1-alpha);
	chargepred(i) = DT;
    DT = .8*DT;
end

% Rappresentazione della curva di wohler
%close all;
figure;
loglog([timepred,time],[chargepred,charge],'*');
hold on
%loglog(time,charge,'x');
%axis([1e4 1e6 1 400]);
% loglog([1e4 1e6], [T_0 T_0],'k','LineWidth',2)
plot(Jabbado(:,3), Jabbado(:,1)-Jabbado(:,2),'*r');
plot(Jabbado(1:n0,3), Jabbado(1:n0,1)-Jabbado(1:n0,2),'o');
axmax = max([time'; Jabbado(:,3)]);
axmin = min([time'; Jabbado(:,3)]);
axmax = 1.1*axmax;
axmin = .9*axmin;
aymax = max(Jabbado(:,1) - Jabbado(:,2));
aymin = min(Jabbado(:,1) - Jabbado(:,2));
aymax = 1.5*aymax;
aymin = aymin/1.5;
axis([axmin axmax aymin aymax])
figure;
loglog(Jabbado(:,3), time, 'r*');
hold on
loglog(Jabbado(1:n0,3), time(1:n0), '*');
fplot(@(x) x,[1e3 1e7]);
fplot(@(x) 2*x,[1e3 1e7]);
fplot(@(x) x/2,[1e3 1e7]);
axis([axmin axmax axmin axmax]);

function [value,isterminal,direction] = events(t,y)
% The exact threshold is given by the value for which f = infty.
value = y - (Gamma_0 + .99*(beta*(mu+c)/alpha)^(1/(alpha-1)));
isterminal = 1;
direction = 1;
end
end