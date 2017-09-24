Files
=====
Jabbado8.mat 		-> data from Jabbado(Jabbado84=AL6082T6, Jabbado=30NCD16,Dubar=30NCD16 90 degree)
data.mat		-> other data
data_preprocessing8.m 	-> preprocessing
wohlerfit7.m		-> model fitting
wohler_param7.m		-> model predictions and plotting
confronto.m		.> plotting script example for experimental data
.txt files		-> source data already formatted for the program to work in the files data.mat and Jabbado8.mat


% Commands and operations to run the code
% ===========================================
% Load data from 'Jabbado8.mat'
load Jabbado8.mat

% Preprocessing and saving data in the 'Jabbado8.mat' file
Jabbado_pp = data_preprocessing8(Jabbado, f_1, t_1);
Jabbado84_pp = data_preprocessing8(Jabbado84, f84_1, t84_1);
save Jabbado8.mat Jabbado_pp Jabbado84_pp -append

% Designing experimental points by running scripts 'confronto.m'
run('confronto.m')
% Model fitting
% param = wohlerfit7(Jabbado_pp, n_number_of_observations_in_fitting)
% n_number_of_observations_in_fitting = 11 in Jabbado, check in files for this number.
param = wohlerfit7(Jabbado_pp, 11)
% Representation of results
[time,charge] = wohler_param7(param, Jabbado_pp,1);