clear all; 
pathname = 'C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Data\behavior\rawData';

files =dir([pathname  filesep '*_*.mat']); % checking only for mat files in the specified subfolder
filenames={files.name};

filename = 'UTNI101_2018-09-07_01.mat';
load([pathname filesep filename]);
current_trial = trial; 
filename = 'UTNI102_2018-09-08_01.mat';
load([pathname filesep filename]);

i = 1;
for file_index=numel(current_trial)+1: (numel(current_trial)-1 + numel(trial))
current_trial(file_index) = trial(i);
i = i+1;
end
trial = current_trial;
filename = 'UTNI101_2018-09-07_01_02.mat';
save([pathname filesep filename],'SETTINGS','trial');



%