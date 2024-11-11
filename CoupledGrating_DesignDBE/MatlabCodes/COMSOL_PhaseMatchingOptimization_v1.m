%% Parameter optimization
clear all; close all; clc

%% Set up

% Initial Parameters
period = 267; % In nm
owgWidth = 200;
gap = 200;
teethWidth = 214;
teethShift = period/2;
teethThickness = period/2;
claddingLength = 1500;

freq_DBE = 193.19; % In THz

c_const = 299792458; % In m/s

% Convert parameters to string for COMSOL and add units
period_str = [num2str(period), '[nm]'];
owgWidth_str = [num2str(owgWidth), '[nm]'];
gap_str = [num2str(gap), '[nm]'];
teethWidth_str = [num2str(teethWidth), '[nm]'];
teethShift_str = [num2str(teethShift), '[nm]'];
teethThickness_str = [num2str(teethThickness), '[nm]'];
claddingLength_str = [num2str(claddingLength), '[nm]'];

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
filename2 = 'COMSOL_SHF_effRefrIndexl_v1_MATLAB'; % SH frequency file

%% Optimization Function
clc
% Define initial guess for the parameters
x0 = [period, owgWidth, gap, teethWidth]'; 

% Define optimization bounds (can be modified as needed)
lb = [100, 100, 100, 100]; % Lower bounds
ub = [500, 480, 300, 300]; % Upper bounds

global iterNum cost_DBE_vec cost_phaseMatching_vec params_vec k_step_FF k_step_SHF  kd_range_FF kd_range_SHF
iterNum = 0;
cost_DBE_vec = [];
cost_phaseMatching_vec = [];
params_vec = [];
k_step_FF = 101;
k_step_SHF = 101;
kd_range_FF = 0.05;
kd_range_SHF = 0.2;

% Define the optimization function that will be minimized
% options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp','StepTolerance', 0.5);
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp', ...
    'StepTolerance', 1e-2, ...  % Increase step size tolerance
    'OptimalityTolerance', 1e-6, ...  % Increase optimality tolerance
    'FiniteDifferenceStepSize', 1e-5);  % Adjust finite difference step size

[x_opt, fval] = fmincon(@(x) objectiveDBE(x, filename1, filename2), x0, [], [], [], [], lb, ub, [], options);

disp(['Optimal Parameters: ', num2str(x_opt')]);
disp(['Objective Function Value: ', num2str(fval)]);

%% Objective Function Definition
function cost = objectiveDBE(params, filename1, filename2)
    global iterNum cost_DBE_vec cost_phaseMatching_vec params_vec k_step_FF k_step_SHF  kd_range_FF kd_range_SHF
    display(['Iteration number:', num2str(iterNum)]);
    iterNum = iterNum + 1;
    
    c_const = 299792458; % In m/s
    
    period = params(1);
    owgWidth = params(2);
    gap = params(3);
    teethWidth = params(4);
    params_vec = [params_vec, params];
    
    period_str = [num2str(period), '[nm]'];
    owgWidth_str = [num2str(owgWidth), '[nm]'];
    gap_str = [num2str(gap), '[nm]'];
    teethWidth_str = [num2str(teethWidth), '[nm]'];
    
    display(['period = ', period_str]);
    display(['owgWidth = ', owgWidth_str]);
    display(['teethWidth = ', teethWidth_str]);
    display(['gap = ', gap_str]);
    
    % Load and set the parameters for the first model (Fundamental frequency)
    display('Running FF sim')
    model1 = mphload(filename1);
    model1.param.set('period', period_str);
    model1.param.set('owgWidth', owgWidth_str);
    model1.param.set('gap', gap_str);
    model1.param.set('teethWidth', teethWidth_str);
    model1.param.set('k_step',num2str(k_step_FF));
    model1.param.set('kd_range',num2str(kd_range_FF));
    
    % Run the simulation and extract the effective refractive index for the fundamental frequency
    model1.study('std1').run;
    kx_f1 = mphglobal(model1, 'kx', 'dataset', 'dset2', 'outersolnum', 'all');
    freq_vec_f1 = mphglobal(model1, 'ewfd.freq', 'dataset', 'dset2', 'outersolnum', 'all');
    model1.sol('sol1').clearSolution(); % Erase solutions for next iteration
    
    
    % Implement DBE at the FF
    kx_f1_m2 = kx_f1(2,:);
    freq_f1_m2 = freq_vec_f1(2,:);
    
    vg_f1 = diff(freq_f1_m2)./diff(kx_f1_m2);
    GVD_f1 = diff(vg_f1)./diff(kx_f1_m2(2:end));
    
    kx_f1_center = kx_f1_m2(round(length(kx_f1_m2)/2));
    freq_DBE = freq_f1_m2(round(length(kx_f1_m2)/2));
    
    effIndex_f1_m2 = kx_f1_center./(2*pi*freq_DBE/c_const);
    
    % Cost function: Minimize GVD at the DBE frequency
    cost_DBE = min(abs(GVD_f1));
    
    freq_DBE_str = [num2str(freq_DBE*1e-12),'[THz]'];
    
    % Load and set the parameters for the second model (Second harmonic frequency)
    display('Running SHF sim')
    model2 = mphload(filename2);
    model2.param.set('period', period_str);
    model2.param.set('owgWidth', owgWidth_str);
    model2.param.set('gap', gap_str);
    model2.param.set('teethWidth', teethWidth_str);
    model2.param.set('f_DBE',freq_DBE_str); % Update DBE frequency so that SHF=2*f_DBE.
    model2.param.set('k_step',num2str(k_step_SHF));
    model2.param.set('kd_range',num2str(kd_range_SHF));
    
    % Run the simulation and extract the effective refractive index for the second harmonic frequency
    model2.study('std1').run;
    kx_f2 = mphglobal(model2, 'kx', 'dataset', 'dset2', 'outersolnum', 'all');
    freq_vec_f2 = mphglobal(model2, 'ewfd.freq', 'dataset', 'dset2', 'outersolnum', 'all');
    model2.sol('sol1').clearSolution(); % Erase solutions for next iteration
    
    display('Finished running sims')
    
    kx_f2_m2 = kx_f2(2,:);
    kx_f2_center = kx_f2(round(length(kx_f2_m2)/2));
    freq_SH = 2*freq_DBE;
    freq_SH = freq_vec_f2(1,round(length(kx_f2_m2)/2));
    effIndex_f2 = kx_f2_center./(2*pi*freq_vec_f2(:,round(length(kx_f1_m2)/2)+1)/c_const);
    
    % Calculate values at the SHF
    freq_f2_center = freq_vec_f2(:,round(length(kx_f2_m2)/2)+1); % It's a vector. f_center for all modes
    effIndices_f2 = kx_f2_center./(2*pi*freq_f2_center/c_const);
    
    % Implement equal effective modal indices
    [cost_phaseMatching, index_phaseMatching] = min(abs(effIndices_f2 - effIndex_f1_m2));
    freq_phaseMatching = freq_f2_center(index_phaseMatching); % Frequency at which we have the phase matching
    
    % Double-check
    Delta_freq_SHG = 2*freq_DBE - freq_phaseMatching; % Double check SHG occurs at 2*f_DBE
    Delta_kx_SHG = kx_f2_center - pi/(period*1e-9);  % Double check center kx is pi/d

    % Find the index for the frequency closest to 2*freq_DBE
    [~, idx_SH] = min(abs(freq_vec_f2(2, :) - 2*freq_DBE)); % Use mode 2 (2nd row)

    % FIX THIS: Enforce kx_f2 at 2*freq_DBE to be non-zero
    kx_f2_at_SH = kx_f2(2, idx_SH);
    penalty_kx = 1e3 * max(0, 1e-3 - abs(kx_f2_at_SH)); % Add a penalty if kx is too close to zero

    % Total cost with the new penalty
    weight_DBE = 1;
    weight_phaseMatching = 10;
    weight_penalty = 1;
    cost = weight_DBE * cost_DBE + weight_phaseMatching * cost_phaseMatching + weight_penalty * penalty_kx;

    display(['DBE cost: ', num2str(cost_DBE)]);
    display(['Phase matching cost: ', num2str(cost_phaseMatching)]);
    display(['kx penalty: ', num2str(penalty_kx)]);
    display(['Total cost: ', num2str(cost)]);


    % Store global
    cost_DBE_vec = [cost_DBE_vec; cost_DBE];
    cost_phaseMatching_vec = [cost_phaseMatching_vec; cost_phaseMatching];
    
    % Plots
    % figure(1)
    % plot(kx_f1, freq_vec_f1)
    % grid on
    % xlabel('kx')
    % ylabel('f')
    % title('FF Dispersion')
    %
    % figure(2)
    % plot(kx_f2, freq_vec_f2)
    % grid on
    % xlabel('kx')
    % ylabel('f')
    % title('SHF Dispersion')

end
