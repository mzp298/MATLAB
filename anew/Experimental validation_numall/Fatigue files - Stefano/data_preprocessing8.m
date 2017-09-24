function new_data = data_preprocessing8(data, f, t)
% Given a matrix of data formatted as in Jabbato page 97 and the fatigue
% limits under fully reversed tension/compression (f) and torsion (t) tests,
% it computes for each row:
% - the amplitude of the macroscopic shear stress on the critical plane T
% - the amplitude of the macroscopic shear stress on the critical plane for
%     the limit loading T_0
% - the numeber of cycles to failure N
%
% WARNING: works only with in phase or 90° phase discrepancy
% tension/compression or torsion tests

[n,m] = size(data);
new_data = zeros(0,3);
for i = 1:n
    % torsion/tension in phase test
    if data(i,7) == 0
        T = sqrt(data(i,3)^2/4 + data(i,4)^2);
        N = data(i,5);
        T_0 = (t*f - (t - f/2)*data(i,2))/(f + (t - f/2)/sqrt(data(i,3)^2/4 + data(i,4)^2)*data(i,3));
        new_data = [new_data; T, T_0, N];
    end
    % torsion/tension 90° out-of-phase test
    if data(i,7) == 90
%         if data(i,4) >= data(i,3)
%             T = data(i,4);
%             N = data(i,5);
%             T_0 = (t*f - (t - f/2)*data(i,2))/(t*f + (t - f/2)*data(i,3))*data(i,4);
%             new_data = [new_data; T, T_0, N];
%         else
            tt = data(i,4);
            ss = data(i,3);
%            Sqrt_K = (ss^2 + tt^2)/(2*ss);
%            L = tt*(ss^2 + tt^2)*sqrt(ss^2 - tt^2)/(2*sqrt(2)*ss^2);
%            T = sqrt(Sqrt_K^2/2 + sqrt(Sqrt_K^4/4 - L^2));
            T = sqrt(ss^2 + tt^2)/(2*sqrt(2)*ss)*sqrt(ss^2 + tt^2 + abs(ss^2-3*tt^2));
            N = data(i,5);
%            H = sqrt(2*pi)*sqrt(ss^2+tt^2)/sqrt(ss^2+tt^2+abs(ss^2-3*tt^2));
%            T_lim = (t*f - (t-f/2)^data(i,2))/(f*(ss^2+tt^2) + 2*(t-f/2)*ss^2)*sqrt(pi)*(ss^2+tt^2);
%            T_0 = T_lim/H;
            T_0 = (t*f - (t - f/2)*data(i,2))/(f*(ss^2 + tt^2) + 2*(t-f/2)*ss^2)*sqrt((ss^2+tt^2)/2)*sqrt(ss^2+tt^2+abs(ss^2-3*tt^2));
            new_data = [new_data; T, T_0, N];
%         end
    end
end

[n,m] = size(new_data);
i = 1;
while i<=n
    for j = 1:m
        if isinf(new_data(i,j)) || isnan(new_data(i,j))% || new_data(i,1)-new_data(i,2) > 35
            new_data = [new_data(1:i-1,:); new_data(i+1:n,:)];
            i = i-1;
            n = n-1;
            break
        end
    end
    i = i+1;
end

%new_data(:,2) = sqrt(new_data(:,3).^2/4 + new_data(:,4).^2);
%new_data = [new_data(:,1:2), new_data(:,5)];