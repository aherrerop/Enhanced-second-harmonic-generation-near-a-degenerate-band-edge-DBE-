%% Connect with COMSOL
clear all; close all; clc
% % Launch COMSOL server
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
% path=pwd; % Get the folder path of this m.file
% filename = strcat('Eigenmode_MATLAB_CG_DBE_v1.mph'); % Name of new COMSOL file
% status = copyfile('CG_Eigenmode_v1.mph', filename); % Name of the template. DO NOT CHANGE
filename = 'COMSOL_DBE_Search_Model_v1';
fullname=[Currentdir '\' filename];


%% Optimize for the DBE
clc
    % Define global variables
global model Function_Counter Results_filename

Results_filename = 'Results_v1.csv';

cd(Currentdir);
model = mphload(filename);

%%
clc
mpheval(model,'ewfd.freq','dataset','dset2','outersolnum',1);

A = mpheval(model,'ewfd.Ex','dataset','dset2','outersolnum',1);

B = mphsolinfo(model);

C = mpheval(model,'ewfd.Ez')

% The value of the field is in C.d1 (it is evaluated at each mesh node, for
% each eigenfrequency)





