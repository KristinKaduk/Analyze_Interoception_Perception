clear all; 
pathname = 'Y:\Projects\Interoception\Perceptual_Certainty\Data\behavior\rawData';

files =dir([pathname  filesep '*_*.mat']); % checking only for mat files in the specified subfolder
filenames={files.name};

filename = 'KEHA031_2018-12-20_01_02.mat';
load([pathname filesep filename]);
current_trial = trial; 
filename = 'KEHA031_2018-12-20_05.mat';
load([pathname filesep filename]);

i = 1;
for file_index=numel(current_trial)+1: (numel(current_trial)+ numel(trial))
current_trial(file_index) = trial(i);
i = i+1;
end
trial = current_trial;
filename = 'KEHA031_2018-12-20_01_02_05.mat';
save([pathname filesep filename],'SETTINGS','trial');



%