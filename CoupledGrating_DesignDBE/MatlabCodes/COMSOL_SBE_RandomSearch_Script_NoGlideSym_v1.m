%% COMSOL SBE random search script
clear all; close all; clc

%% Set up

% Parameters
period = 350; % In nm
owgWidth = 250;
gap = 300;
teethWidth = 200;
teethShift = period/2;
teethThickness = period/2;
claddingLength = 1000;
kx = 0;
band = 1;
sigma_Ave = 1e20;
sigma_StDev = 1e20;
DBEfreq = 193.54; % In THz

% Parameter ranges for random generation
period_range = [250, 500]; % Example range
owgWidth_range = [100, 500];
gap_range = [100, 500];
teethWidth_range = [100, 500];

%% Connect with COMSOL
clc
% Launch COMSOL server
Currentdir = 'D:\AlbertHerreroParareda\DBE-SHG\2D_DDFB_SHG';
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

% Copy and open file
filename = 'COMSOL_DBE_Search_Model_v1';
fullname = [Currentdir '\' filename];

%% Optimize for the DBE
clc
% Define global variables
global model Function_Counter Results_filename

Results_filename = 'Results_v1.csv';

cd(Currentdir);
model = mphload(filename);

%%

Istrue = 0;
kk = 0;

while Istrue == 0
    % Generate random parameter values within the predefined range
    period = randi(period_range);
    owgWidth = randi(owgWidth_range);
    gap = randi(gap_range);
    teethWidth = randi(teethWidth_range);

    % Ensure teethShift and teethThickness sum to match the period constraint
    teethThickness = randi([50, round(period-1)]); % Assuming teethThickness must be at least 50 nm and up to the period
    teethShift = randi([0, round(period - teethThickness - 1)]); % teethShift <= period - teethThickness

    % Convert parameters to string for COMSOL
    period_str = num2str(period);
    owgWidth_str = num2str(owgWidth);
    gap_str = num2str(gap);
    teethWidth_str = num2str(teethWidth);
    teethShift_str = num2str(teethShift);
    teethThickness_str = num2str(teethThickness);
    claddingLength_str = num2str(claddingLength);

    % Update variables in COMSOL
    model.param.set('owgWidth', owgWidth_str);
    model.param.set('teethWidth', teethWidth_str);
    model.param.set('teethShift', teethShift_str);
    model.param.set('period', period_str);
    model.param.set('gap', gap_str);
    model.param.set('teethThickness', teethThickness_str);
    model.param.set('claddingLength', claddingLength_str);

    % Run the parametric sweep
    model.study('std1').run;

    disp(strcat('Running model: ',num2str(kk)));
    kk = kk + 1;

    % Collect the eigenfrequencies
    Freq = zeros(4,10);
    for ii = 1:10
        a = mpheval(model,'ewfd.freq','dataset','dset2','outersolnum',ii);
        Freq(:,ii) = a.d1(1:4,1);
    end

    % Post-processing
    for jj = 1:4
        Average_F(jj) = mean(Freq(jj,:));
        [F_max(jj,:), index_max(jj,:)] = max(Freq(jj,:));
        [F_min(jj,:), index_min(jj,:)] = min(Freq(jj,:));
        F_ave_max(jj) = F_max(jj,1) - F_min(jj,1);
        F_std_dev(jj) = std(Freq(jj,:));
    end

    % Find peaks in the second eigenfrequency
    [peaks, locs] = findpeaks(Freq(2,:));

    % Check the condition for peaks
    if length(peaks) > 1
        % Store parameter values if multiple peaks are found
        parameters = [period, owgWidth, gap, teethWidth, teethShift, teethThickness];
        % Save or display the parameters as needed
        disp('Parameters with multiple peaks:');
        disp(parameters);
        
        % Extract the frequency at kd = pi (center)
        kd_index = find(Freq(1,:) == pi);
        if ~isempty(kd_index)
            center_freq = Freq(2, kd_index);
            disp(['Frequency at kd = pi: ', num2str(center_freq)]);
        end

        % Calculate delta f (difference between frequencies at the peaks)
        delta_f = max(peaks) - min(peaks);
        disp(['Delta f (difference between peak frequencies): ', num2str(delta_f)]);

        Istrue = 1; % Exit loop if multiple peaks are found
        elseif length(peaks) == 1
    disp('No SBE! Only one peak found.');
else
    disp('No SBE! No peaks found.');
    end
end

disp('Optimization complete.');
