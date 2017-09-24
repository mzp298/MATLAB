for  i=1:length(NFben) %experimental points index
    n=1;       
    tensor = [stressben(i)*sind(n*360/stepnumber)+m 0 0 ;...
        0 0 0 ;...
        0 0 0 ];
    run('Damiter1.m')
    g=1; %first reference alp, W, n index when generate
    while n<cycles*stepnumber
        tensor = [stressben(i)*sind(n*360/stepnumber)+m, 0, 0 ;...
            0, 0, 0 ;...
            0, 0, 0 ; ];
        hydro=1/3*trace(tensor);
        dev1=tensor-hydro*eye(3);
        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
        dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
        dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
        tensor = [stressben(i)*sind((n+1)*360/stepnumber)+m, 0, 0 ;...
            0, 0,  0 ;...
            0, 0,  0 ; ];
        run('Damiter2.m')
        n=n+1;
    end
    
%----------scalar iteration
e=1; %first D index
j=ceil(length(alp_ref)/cycles); %first reference alp, W, n index when iterate(after adaptation)
D= 1e-16; 
while D<1 %-----------the optimal time steps can be iterated with scalar
    D= (D^(1-alp_ref(j))+(1-alp_ref(j)).*W_ref(j)/W0).^(1/(1-alp_ref(j))); % implicit
    j=j+1;
    if j+1>=length(alp_ref)
        j=ceil(length(alp_ref)/cycles);
    end
    e=e+1;
end
NFben_num(i)=e/(length(n_ref)/cycles)
end