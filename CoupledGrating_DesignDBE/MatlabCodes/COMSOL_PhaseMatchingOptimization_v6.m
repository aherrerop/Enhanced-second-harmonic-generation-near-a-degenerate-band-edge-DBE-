%% Parameter optimization
clear all; close all; clc

%% Set up

% Initial Parameters
period = 267.2766; % In nm
owgWidth = 200.1254;
gap = 196.7525;
teethWidth = 214.1481;
teethShift = period/2;
teethThickness = period/2;
claddingLength = 1500;

% freq_DBE = 193; % In THz

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


%% Useful commands
clc
cd(Currentdir);
import com.comsol.model.*;
import com.comsol.model.util.*;

% Load COMSOL models
filename1 = 'COMSOL_FF_effRefrIndexl_v1_MATLAB'; % Fundamental frequency file
filename2 = 'COMSOL_SHF_effRefrIndexl_v1_MATLAB'; % SH frequency file

%% Optimization Function
clc
% Define initial guess for the parameters
x0 = [period, owgWidth, gap, teethWidth, teethShift, teethThickness]'; 

% Define optimization bounds (can be modified as needed)
lb = [100, 100, 100, 100, 0, 100]; % Lower bounds
ub = [500, 480, 300, 300, 250, 250]; % Upper bounds

global iterNum cost_DBE_vec cost_phaseMatching_vec params_vec k_step_FF k_step_SHF  kd_range_FF kd_range_SHF
iterNum = 0;
cost_DBE_vec = [];
cost_phaseMatching_vec = [];
params_vec = [];
k_step_FF = 101;
k_step_SHF = 101;
kd_range_FF = 0.05;
kd_range_SHF = 0.9;

% Define the optimization function that will be minimized
% options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp','StepTolerance', 0.5);
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp', ...
    'StepTolerance', 1e-2, ...  % Increase step size tolerance
    'OptimalityTolerance', 1e-6, ...  % Increase optimality tolerance
    'MaxIterations', 1000, ...  % Set a higher max iteration count
    'FiniteDifferenceStepSize', 1e-3);  % Adjust finite difference step size

[x_opt, fval] = fmincon(@(x) objectiveDBE(x, filename1, filename2), x0, [], [], [], [], lb, ub, [], options);

disp(['Parameters : period, owgWidth, gap, teethWidth, teethShift, teethTHickness']);
disp(['Optimal Parameters: ', num2str(x_opt')]);
disp(['Objective Function Value: ', num2str(fval)]);

%% Objective Function Definition
function cost = objectiveDBE(params, filename1, filename2)
    global iterNum cost_DBE_vec cost_phaseMatching_vec params_vec k_step_FF k_step_SHF kd_range_FF kd_range_SHF
    display(['Iteration number:', num2str(iterNum)]);
    iterNum = iterNum + 1;

    c_const = 299792458; % In m/s

    % Extract parameters from input
    period = params(1);
    owgWidth = params(2);
    gap = params(3);
    teethWidth = params(4);
    teethShift = params(5);
    teethThickness = params(6);
    params_vec = [params_vec, params];

    % Convert parameters to strings
    period_str = [num2str(period), '[nm]'];
    owgWidth_str = [num2str(owgWidth), '[nm]'];
    gap_str = [num2str(gap), '[nm]'];
    teethWidth_str = [num2str(teethWidth), '[nm]'];
    teethShift_str = [num2str(teethShift), '[nm]'];
    teethThickness_str = [num2str(teethThickness), '[nm]'];

    % Define the parameter names
    param_names = {'period', 'owgWidth', 'gap', 'teethWidth', 'teethShift', 'teethThickness'};

    % Display the parameter names
    fprintf('Parameters: ');
    fprintf('%s\t', param_names{:});
    fprintf('\n');

    % Display the optimal parameter values with specified formatting
    fprintf('Optimal Parameters: ');
    fprintf('%.2f\t', params);
    fprintf('\n');


    try
        % Load and set parameters for the FF model
        display('Running FF sim')
        model1 = mphload(filename1);
        model1.param.set('period', period_str);
        model1.param.set('owgWidth', owgWidth_str);
        model1.param.set('gap', gap_str);
        model1.param.set('teethWidth', teethWidth_str);
        model1.param.set('teethShift', teethShift_str);
        model1.param.set('teethThickness', teethThickness_str);
        model1.param.set('k_step', num2str(k_step_FF));
        model1.param.set('kd_range', num2str(kd_range_FF));

        % Run the FF simulation and extract the dispersion data
        model1.study('std1').run;
        kx_f1 = mphglobal(model1, 'kx', 'dataset', 'dset2', 'outersolnum', 'all');
        freq_vec_f1 = mphglobal(model1, 'ewfd.freq', 'dataset', 'dset2', 'outersolnum', 'all');
        model1.sol('sol1').clearSolution(); % Clear solutions after running

        % DBE at FF: Extract GVD and minimize it
        kx_f1_m2 = kx_f1(2, :);
        freq_f1_m2 = freq_vec_f1(2, :);
        vg_f1 = diff(freq_f1_m2) ./ diff(kx_f1_m2);
        GVD_f1 = diff(vg_f1) ./ diff(kx_f1_m2(2:end));
        cost_DBE = min(abs(GVD_f1)); % Minimize GVD

        % Calculate effective index at FF
        kx_f1_center = kx_f1_m2(round(length(kx_f1_m2) / 2));
        freq_DBE = freq_f1_m2(round(length(kx_f1_m2) / 2));
        effIndex_f1_m2 = kx_f1_center / (2 * pi * freq_DBE / c_const);

        % SHF Simulation at 2*f_DBE
        freq_DBE_str = [num2str(freq_DBE * 1e-12), '[THz]'];
        display('Running SHF sim')
        model2 = mphload(filename2);
        model2.param.set('period', period_str);
        model2.param.set('owgWidth', owgWidth_str);
        model2.param.set('gap', gap_str);
        model2.param.set('teethWidth', teethWidth_str);
        model2.param.set('teethShift', teethShift_str);
        model2.param.set('teethThickness', teethThickness_str);
        model2.param.set('f_DBE', freq_DBE_str); % Update DBE frequency for SHF
        model2.param.set('k_step', num2str(k_step_SHF));
        model2.param.set('kd_range', num2str(kd_range_SHF));

        % Run the SHF simulation and extract the dispersion data
        model2.study('std1').run;
        kx_f2 = mphglobal(model2, 'kx', 'dataset', 'dset2', 'outersolnum', 'all');
        freq_vec_f2 = mphglobal(model2, 'ewfd.freq', 'dataset', 'dset2', 'outersolnum', 'all');
        model2.sol('sol1').clearSolution(); % Clear solutions after running

        display('Finished running sims')
        
        % RBE at SHF: Extract group velocity and minimize 

        SH_freq = 2*freq_DBE;
        A = find(freq_vec_f2 < SH_freq + 2e12);
        B = find(freq_vec_f2 > SH_freq);
        freq_box = freq_vec_f2(:,A:B);
        kx_box = kx_f2(:,A:B);

        % Calculate v_g near SHF & minimize it
        [numModes, freqLength] = size(freq_box);
        for ii = 1:numModes
            vg_box(ii,:) = diff(freq_box(ii,:)) ./ diff(kx_box(ii,:));
        end
        [cost_RBE, freq_index_RBE] = min(abs(vg_box),[],'all'); % Minimize GVD
        % Convert the linear index to row and column indices
        [row_RBE, col_RBE] = ind2sub(size(vg_box), freq_index_RBE);

        % Retrieve the corresponding effective modal index at the RBE
        kx_f2_RBE = kx_box(row_RBE, col_RBE);
        freq_RBE = freq_box(row_RBE, col_RBE + 1); % +1 to align with diff result
        effIndex_RBE = kx_f2_RBE / (2 * pi * freq_RBE / c_const);
        cost_phaseMatching = abs(effIndex_f1_m2 - effIndex_RBE);

    catch ME
        disp(['Error encountered: ', ME.message]);
        cost_DBE = 4e6;
        cost_RBE = 2e6;
        cost_phaseMatching = 1e6;
    end

    % Total cost
    weight_DBE = 1;
    weight_RBE = 1e-6;
    weight_phaseMatching = 1e2;
    cost = weight_DBE * cost_DBE + weight_RBE * cost_RBE + weight_phaseMatching * cost_phaseMatching;

    display(['DBE cost: ', num2str(cost_DBE)]);
    display(['Phase matching cost: ', num2str(cost_phaseMatching)]);
    display(['RBE cost: ', num2str(cost_RBE)]);
    display(['Total cost: ', num2str(cost)]);

    % Store global values
    cost_DBE_vec = [cost_DBE_vec; cost_DBE];
    cost_phaseMatching_vec = [cost_phaseMatching_vec; cost_phaseMatching];
end
