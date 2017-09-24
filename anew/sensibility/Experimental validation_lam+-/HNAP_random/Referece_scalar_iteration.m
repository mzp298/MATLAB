%------------repeated code for random loading with blocks as one cycle reference-----------
[force]=textscan(fid,'%*s%*s%s%*s%*s',repetition,'headerlines',5);
fclose(fid);
stress=1000*str2double(strrep(force{1,1},',','.')).*area^-1; %Pa
stress=[stress(1)/ari; stress];
length(stress)
for i=2:length(stress)
    stress11(1+ari*(i-2):1+ari*(i-1))=linspace(stress(i-1),stress(i),ari+1);
end
 length(stress11)
clear force;

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

%% 

%----------scalar iteration
e=1; %first D index
j=1; %first reference alp, W, n index when iterate(after adaptation)
D= 1e-16;
while D<1 %-----------the optimal time steps can be iterated with scalar
%     D= (D^(1-alp_ref(j))+(1-alp_ref(j)).*W_ref(j)/W0).^(1/(1-alp_ref(j))); % implicit
    D=D+D^alp_ref(j)*W_ref(j)/W0;
    j=j+1;
    if j+1>=length(alp_ref)
        j=1;
    end
    e=e+1;
end
remain_index=rem((e-1),length(alp_ref));
nF=floor((e-1)/length(alp_ref))*repetition+n_ref(remain_index);
disp(['Number of test points is ' num2str(nF) ' points.']);
disp(['Error is ' sprintf('%2.2f%%', ((NF-nF)/NF)*100) ' .']);

