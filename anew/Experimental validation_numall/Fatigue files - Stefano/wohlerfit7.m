function opt_param = wohlerfit7(Jabbado, n)
% Fit del modello sui dati contenuti nella matrice in input usando le prime
% n osservazioni.

% Il punto di partenza dell'algoritmo è completamente arbitrario.

Gamma_0 = 66;
mu = 1400000;
c = 1400000;
mu_c = mu + c;
alpha = 4;
beta = 70000;
T_0 = 400;
param = [Gamma_0, mu_c, alpha, beta];
save tmp.mat Jabbado n

opt_param = fminunc(@min_square,param,optimset('Maxiter',10000000,'TolFun', 1e-9, 'TolX', 1e-9, 'MaxFunEvals', 100000000));

function square_error = min_square(param)

Jabbado = load('tmp.mat', 'Jabbado', 'n');
n = Jabbado.n;
Jabbado = Jabbado.Jabbado(1:n,:);

Gamma_0 = param(1);
mu_c = param(2);
alpha = param(3);
beta = param(4);

square_error = 0;

% if (alpha > 6)
    for i = 1:size(Jabbado,1)
        % Legge costitutiva

         T = Jabbado(i,1);
         T_0 = Jabbado(i,2);
         
        tmp = (2/beta/(T - T_0))^((alpha-1)/alpha);
        gamma = (pi*tmp + sqrt(pi*alpha^2/(alpha-1)/(mu + c)/beta*tmp*sin(pi/alpha)))./...
            (pi*tmp - alpha^2/(alpha-1)/beta/(mu + c)*sin(pi/alpha));
        t_approx = (gamma + 1)*pi/(alpha*sin(pi/alpha))*((mu + c)/4)*...
            beta^(1/alpha)/gamma*((T-T_0)/2)^((1-alpha)/alpha) -...
            alpha/4/(alpha-1)/(gamma-1) -...
            (mu + c)*beta/4/(alpha - 1)*Gamma_0^(1-alpha);

% Qui possiamo provare due strade diverse

        disp('Minimizing the sum of squares');
        square_error = square_error + abs(t_approx - Jabbado(i,3));
        
%        disp('Minimizing the sum of relative errors squares');
%        square_error = square_error + ((t_approx - Jabbado(i,3))/Jabbado(i,3))^2;   

    end
square_error

function [value,isterminal,direction] = events(t,y)
% La soglia esatta é data dal valore per cui f = infty.
value = y - (Gamma_0 + .99*(beta*(mu+c)/alpha)^(1/(alpha-1)));
isterminal = 1;
direction = 1; 
end
end
end
