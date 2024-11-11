%% Parameter optimization
clear all; close all; clc


c_const = 299792458; % In m/s


% Connect with COMSOL
clc
% Launch COMSOL server
Currentdir = 'D:\AlbertHerreroParareda\DBE-SHG\2D_DDFB_SHG\Eigenmode\';
cd ('C:\Program Files\COMSOL\COMSOL53a\Multiphysics\bin\win64');
system('comsolmphserver.exe &');
cd(Currentdir);
% Establish connection with COMSOL
Currentdir = pwd;

cd('C:\Program Files\COMSOL\COMSOL53a\Multiphysics\mli');
mphstart(2036);
cd(Currentdir);

% Useful commands
import com.comsol.model.*;
import com.comsol.model.util.*;

% Load COMSOL models
filename1 = 'COMSOL_FF_effRefrIndexl_v1_MATLAB'; % Fundamental frequency file
% filename2 = 'COMSOL_SHF_effRefrIndexl_v1_MATLAB'; % SH frequency file

% Load and set the parameters for the first model (Fundamental frequency)
model1 = mphload(filename1);


%%
% Run the simulation and extract the effective refractive index for the fundamental frequency
model1.study('std1').run;

%
clc
% Extract kx values for the parametric sweep
kx_f1 = mphglobal(model1, 'kx', 'dataset', 'dset2', 'outersolnum', 'all');
% kx_f1 = kx_f1(1,:);
% Extract the corresponding frequency values
freq_vec_f1 = mphglobal(model1, 'ewfd.freq', 'dataset', 'dset2', 'outersolnum', 'all'); % Ensure the 'outersolnum' is 'all' to get all frequencies.
% freq_vec_f1 = freq_vec_f1(1,:);
% Calculate effective refractive index
effIndex_f1 = kx_f1 ./ (2 * pi * freq_vec_f1 / c_const);

% Display sizes for debugging
size(kx_f1)
size(freq_vec_f1)
size(effIndex_f1)

%%
% Load and set the parameters for the second model (Second harmonic frequency)
model2 = mphload(filename2);
model2.param.set('period', period_str);
model2.param.set('owgWidth', owgWidth_str);
model2.param.set('gap', gap_str);
model2.param.set('teethWidth', teethWidth_str);

% Run the simulation and extract the effective refractive index for the second harmonic frequency
model2.study('std1').run;
kx_f2 = mphglobal(model2, 'kx', 'dataset', 'dset1');
freq_vec_f2 = mphglobal(model2, 'ewfd.freq', 'dataset', 'dset1');
effIndex_f2 = kx_f2./(2*pi*freq_vec_f2/c_const);

% Find the flatness of the dispersion curve at the band edge (DBE condition)

% Compute the slope of the bands at k=pi/d for f1 (DBE) and f2 (SHG)
slope_f1 = diff(freq_vec_f1) / (kx_f1);
slope_f2 = diff(freq_vec_f2) / (kx_f2);

% Find the flatness around DBE frequency
flatness_f1 = min(abs(slope_f1));
flatness_f2 = min(abs(slope_f2));

% Cost function: Minimize the difference in effective index and maximize flatness (minimize slope)
phase_matching_error = abs(effIndex_f1 - effIndex_f2);
cost = phase_matching_error + 0.1*(flatness_f1 + flatness_f2);

%%
close all
figure(1)
hold on
plot(freq_vec_f1, effIndex_f1)
% plot(freq_vec_f2, effIndex_f2)
hold off
grid on
xlabel('f')
ylabel('n eff')