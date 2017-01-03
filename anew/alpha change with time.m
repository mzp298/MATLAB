% Program to get the Gauss-Legendre Quadrature results


%---------------------alpha-----------------------------
% sig1 = [0 0 0;0 0 0;0 0 load*sin(t);] ;                         %Stress tensor
% p1=1/3*sum(diag(sig1));                                           %Hydraustatic stress
% S1=2*sig1-(1/3*sum(diag(sig1)))*diag([1,1,1]);      % deviatoric stress
% sqrj1=1/2*sqrt(1/2*(S1(1,1)^2+S1(2,2)^2+S1(3,3)^2+2*(S1(1,2)^2)+...
% 2*(S1(1,3)^2)+2*(S1(2,3)^2)));                               % second invariant of the deviatoric stress tensor
% sqrj=max(sqrj1);
% p=max(p1);
% cross=sqrj+ac*p-bc;
% frac=(cross*(sigu-2*sqrj).^-1);
% frac(frac<0)=0;                                      %turn negative values to 0
% alp=1-1000*frac ;                     %characterizes non-linearity of damage accumulation
%---------------------loading paths-----------------------------
% D=1e-15;             %initial damage
% 
% W=0;
% t=step;                  %initial time
% while sign((sin(t+step)-sin(t))*(sin(t)-sin(t-step)))==1  % step up until 1st turning point
%  for i=1:25
%      s(i)= (x(i)/2+1/2).^(1/(1-b));
%          if load*abs(sin(t)-0)-yield*s(i).^-1>=0 
%             Ws(i)=weight(i)*(E-k)*(1+nu)*yield/(2*E*(E+k*nu))*s(i).^-1*abs(load*(sin(t+step)-sin(t)))/step;
%          else
%             Ws(i)=0;
%          end
%  end
%  Wrate = sum(Ws);
%  W=W+Wrate*step;
%  D=D+(1 - (1 - D).^(gam + 1)).^alp*W/WF;
%   hold on;
%  %Energyrate= plot (t,Wrate,'*m');
%   %Energy=plot (t,W,'*b');
%   Damagepaths=plot (t,D,'*r');
%   t=t+step;
% end
% turn(1)= sin(t)
% % for i=1:25 %dissipation rate integrate from 1 to infinity
% % Ws(i)=0.5*weight(i)*((x(i)/2+1/2).^(1/(1-b)) ).^-1;
% % end
% % Wrate = sum(Ws)
% % Wrate=(b-1)/b
% 
% 
% j=0;
% while D<1
% j=j+1;
% t=t+step;
%   while sign((sin(t+step)-sin(t))*(sin(t)-sin(t-step)))==1 %step up until jth turning point
%     for i=1:25
%         s(i)= (x(i)/2+1/2).^(1/(1-b));
%          if load*abs(sin(t)-turn(j))-2*yield*s(i).^-1>=0 %under scale x(i) if there is Ep
%            Ws(i)=weight(i)*(E-k)*(1+nu)*yield/(2*E*(E+k*nu))*s(i).^-1*abs(load*(sin(t+step)-sin(t)))/step;
%          else
%             Ws(i)=0;
%          end
%    end
%  Wrate = sum(Ws);
%  W=W+Wrate*step;
%  D=D+(1 - (1 - D).^(gam + 1)).^alp*W/WF
%  hold on;
%  %Energyrate= plot (t,Wrate,'*m');
%  % Energy=plot (t,W,'*b');
%   Damage=plot (t,D,'*r');
%   t=t+step;
%   end
% turn(j+1)=sin(t);
% end
% %Wcyc=Wp1*(4*t/(2*pi)-1)

