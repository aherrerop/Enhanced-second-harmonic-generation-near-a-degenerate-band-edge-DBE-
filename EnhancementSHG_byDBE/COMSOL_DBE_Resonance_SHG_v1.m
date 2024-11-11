%% DBE resonance search and SHG
clear all; close all; clc

%% Set up

% Parameters
period = 270; % In nm
owgWidth = 200;
gap = 190;
teethWidth = 216.5;
teethShift = period/2;
teethThickness = period/2;
claddingLength = 1500;

freqDBE = 193.043; % In THz

%% Connect with COMSOL
clc
% Launch COMSOL server
Currentdir = 'D:\AlbertHerreroParareda\DBE-SHG\2D_DDFB_SHG\freqDomain\';
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

% Open the file to find the DBE resonance
filename1 = 'freqDom_SHG_DBE_resonance_matlab_v2';
filename2 = 'freqDom_SHG_DBE_calcSHG_matlab_v2';
fullname = [Currentdir '\' filename1];

%% Open the file to apply the SHG at the DBE resonance

cd(Currentdir);
% model1 = mphload(filename1);
% model2 = mphload(filename2);

%%
clc
numberCells_vec = [20:20:200];
freqDBE_resonance = zeros(1,length(numberCells_vec));

for ii = 1:length(numberCells_vec)
    
    numberCells = numberCells_vec(ii);
    numberCells_str = num2str(numberCells);

    display('Number of unit cells is')
    display(num2str(numberCells));

    model1 = mphload(filename1);
    model1.param.set('numberCells', numberCells_str);

    display('Run model 1')
    model1.study('std1').run;
    % Code to extract the DBE resonance: Frequency at which there are
    % maximum fields
    A = mpheval(model1, 'ewfd1.Ez', 'dataset', 'dset1', 'outersolnum', 'all');
    freq_vec = mphglobal(model1, 'ewfd1.freq', 'dataset', 'dset1');
    S21_vec = mphglobal(model1, 'ewfd1.S21', 'dataset', 'dset1');
    numFreqs = length(freq_vec);

    Ez = A.d1;

    % Find peaks below the DBE frequency
    valid_indices = freq_vec <= freqDBE * 1e12;
    Ez_valid = Ez(valid_indices, :);
    S21_valid = S21_vec(valid_indices, :);
    freq_valid = freq_vec(valid_indices);

    [~, peakIdx] = findpeaks(max(abs(Ez_valid), [], 2));
    [~, closestIdx] = min(abs(freq_valid(peakIdx) - freqDBE * 1e12));

    maxFrequency = freq_valid(peakIdx(closestIdx)) * 1e-12; % In THz

    [~, peakIdx_S21] = findpeaks(max(abs(S21_valid), [], 2));
    [~, closestIdx_S21] = min(abs(freq_valid(peakIdx_S21) - freqDBE * 1e12));

    maxFrequency_S21 = freq_valid(peakIdx_S21(closestIdx_S21)) * 1e-12; % In THz

    freqDBE_resonance(ii) = maxFrequency;
    freqDBE_resonance_S21(ii) = maxFrequency_S21;
    freqDBE_resonance_str = num2str(freqDBE_resonance(ii));

    model2 = mphload(filename2);
    model2.param.set('numberCells', numberCells_str);
    model2.param.set('f1',freqDBE_resonance_str);

    display('Run model 2')
    model2.study('std1').run;

    % Store results in cell arrays
    numCellsArray(ii) = numberCells;
    freqDBEResonanceArray(ii) = freqDBE_resonance(ii);
    freqDBEResonanceArray_S21(ii) = freqDBE_resonance_S21(ii);
    fundamentalEzArray{ii} = mpheval(model2,'ewfd1.Ez').d1;
    secondHarmonicEzArray{ii} = mpheval(model2,'ewfd2.Ez').d1;
    positionMatrixArray{ii} = mpheval(model2,'ewfd1.Ez').p;

end

% Save all data to a .mat file
save('comsol_results_v4.mat', 'numCellsArray', 'freqDBEResonanceArray', ...
    'fundamentalEzArray', 'secondHarmonicEzArray', 'positionMatrixArray','freqDBEResonanceArray_S21');
