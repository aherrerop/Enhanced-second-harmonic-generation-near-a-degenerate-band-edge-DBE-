%% COMSOL_DBE_SHG_Postprocessing
clear all; close all; clc

% Load the data from the .mat file
load('comsol_results_v4.mat');

% The variables now available in the workspace are:
% - numCellsArray
% - freqDBEResonanceArray
% - fundamentalEzArray
% - secondHarmonicEzArray
% - positionMatrixArray


% Post-processing
E_in = 1;

for ii = 1:length(numCellsArray)
    convEff_xy{ii} = abs(secondHarmonicEzArray{ii})./abs(E_in);
    [convEff_max{ii}, index_xy_max{ii}] = max(convEff_xy{ii});
    convEff_x{ii} = sum(convEff_xy{ii}, 1);
    x{ii} = positionMatrixArray{ii}(1,:);

end

figure(1)
plot(numCellsArray,cell2mat(convEff_max))
grid on
ylabel('Conversion efficiency vs N')


figure(2)
plot(freqDBEResonanceArray, cell2mat(convEff_max))
grid on
ylabel('Conversion efficiency vs frequency')

figure(3)
plot(cell2mat(x), cell2mat(convEff_x))
grid on
ylabel('Conversion efficiency vs position')


figure(4)
plot(numCellsArray, freqDBEResonanceArray)
grid on
xlabel('Number of unit cells')
ylabel('DBE resonance frequency')

%%
%% Now we start

% N = numCellsArray;
% Q = freqDBEResonanceArray;
% Begin = 1; 
% x_cut = numCellsArray(Begin:end);
% y_cut = freqDBEResonanceArray(Begin:end);
% 
% % ft = fittype('a*x^7+b*x^4+c');
% ft = fittype('a4*x^4 + a0');
% [cf, g] = fit(x_cut',y_cut',ft)

figure(4)
hold on
plot(numCellsArray, freqDBEResonanceArray,'o-')
plot(numCellsArray, freqDBEResonanceArray_S21,'--')
hold off
grid on
xlabel('Number of unit cells')
ylabel('DBE resonance frequency')
legend('frequency with highest fields','frequency with highest transmission (S21)')
set(gca,'fontsize',20)
