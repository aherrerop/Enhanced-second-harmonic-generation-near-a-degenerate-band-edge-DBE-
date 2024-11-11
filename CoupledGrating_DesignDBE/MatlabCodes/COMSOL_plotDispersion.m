%% 
clear all; close all; clc

%% Set up

% Initial Parameters
period = 267.2807; % In nm
owgWidth = 200.1254;
gap = 196.7525;
teethWidth = 214.1482;
teethShift = period/2;
teethThickness = period/2;
claddingLength = 1500;

% Initial Parameters
period = 268.93; % In nm
owgWidth = 199.14;
gap = 200.02;
teethWidth = 217.73;
teethShift = 0;
teethThickness = 100;
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

%% Connect with COMSOL
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

%% Useful commands
import com.comsol.model.*;
import com.comsol.model.util.*;

% Load COMSOL models
filename1 = 'COMSOL_FF_effRefrIndexl_v1_MATLAB'; % Fundamental frequency file
filename2 = 'COMSOL_SHF_effRefrIndexl_v1_MATLAB'; % SH frequency file


%%
clc

k_step_FF = 51;
k_step_SHF = 51;
kd_range_FF = 0.05;
kd_range_SHF = 0.9;

period_str = [num2str(period), '[nm]'];
owgWidth_str = [num2str(owgWidth), '[nm]'];
gap_str = [num2str(gap), '[nm]'];
teethWidth_str = [num2str(teethWidth), '[nm]'];

display(['period = ', period_str]);
display(['owgWidth = ', owgWidth_str]);
display(['gap = ', gap_str]);
display(['teethWidth = ', teethWidth_str]);
display(['teethShift = ', teethShift_str]);
display(['teethThickness = ', teethThickness_str]);

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

%% Load and set the parameters for the second model (Second harmonic frequency)
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

%%
close all
figure(1)
hold on
plot(kx_f1(1,:), freq_vec_f1(1,:))
plot(kx_f1(2,:), freq_vec_f1(2,:))
% plot(kx_f1(3,:), freq_vec_f1(3,:))
yline(193e12)
hold off
grid on
xlabel('$kd$','interpreter','latex')
ylabel('$f$','interpreter','latex')

ylim([192.3e12 193.25e12])
%%
close all
figure(1)
hold on
for ii=1:8
plot(kx_f2(ii,:)*(period*1e-9)/pi, freq_vec_f2(ii,:),'.b','linewidth',0.5)
end
yline(2*193e12)
hold off
grid on
xlabel('$kd/\pi$','interpreter','latex')
ylabel('$f$','interpreter','latex')

%%
% Find the index for the frequency closest to 2*freq_DBE
% Define the target frequency
SH_freq = 2 * freq_DBE;

% Tolerance for finding crossings near SH_freq
tolerance = 1e12;  % Adjust this as needed based on your data scale

% Find indices of elements in freq_vec_f2 close to SH_freq
[rows, cols] = find(abs(freq_vec_f2 - SH_freq) < tolerance);

% Retrieve the corresponding kx and frequency values for each crossing
crossing_kx = kx_f2(sub2ind(size(kx_f2), rows, cols));
crossing_freq = freq_vec_f2(sub2ind(size(freq_vec_f2), rows, cols));

% % Display the results
% disp('kx and frequency values closest to SH_freq:');
% for i = 1:length(crossing_kx)
%     disp(['kx: ', num2str(crossing_kx(i)), ', Frequency: ', num2str(crossing_freq(i) * 1e-12), ' THz']);
% end

% Check if there are any crossings near SH_freq
if isempty(crossing_kx)
    % Apply high penalties if no crossings are found
    cost_confinedMode = 1e6;
    cost_phaseMatching = 1e6;
    disp('No crossing found close to SH_freq. Applying penalty.');
else
    % Enforce kx_f2 at 2*freq_DBE to be different from zero
    if any(abs(crossing_kx - pi / (period * 1e-9)) < 1e-3)
        cost_confinedMode = 1e6; % High cost if kx is too close to zero
    else
        cost_confinedMode = 0;
        disp(['Closest kx values: ', num2str(crossing_kx(:)')]);
    end

    % Calculate effective index at the closest point to 2*freq_DBE
    effIndex_f2_at_SH = crossing_kx ./ (2 * pi * crossing_freq / c_const);

    % Calculate Phase Matching Cost
    cost_phaseMatching = min(abs(effIndex_f1_m2 - effIndex_f2_at_SH), [], 'all');
end


% Total cost
weight_DBE = 1;
weight_phaseMatching = 10;
weight_confinedMode = 1;
cost = weight_DBE * cost_DBE + weight_phaseMatching * cost_phaseMatching + weight_confinedMode * cost_confinedMode;

display(['DBE cost: ', num2str(cost_DBE)]);
display(['Phase matching cost: ', num2str(cost_phaseMatching)]);
display(['Confined mode cost: ', num2str(cost_confinedMode)]);
display(['Total cost: ', num2str(cost)]);

%%

close all
figure(4)
for ii=1:101
plot(freq_vec_f2(:,ii), effIndex_f2_at_SH)
end
yline(effIndex_f1_m2)
xline(2*193e12)
grid on
xlabel('$f$','interpreter','latex')
ylabel('$n_{eff}$','interpreter','latex')
