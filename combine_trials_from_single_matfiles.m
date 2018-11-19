function combine_trials_from_single_matfiles(folder,subfolder)
global SETTINGS

subfolder_dir=dir([folder filesep subfolder filesep '*_*.mat']); % checking only for mat files in the specified subfolder
subfolder_content={subfolder_dir.name};
subfolder_content=sort(subfolder_content);

for file_index=1:numel(subfolder_content)
    load([folder filesep subfolder filesep subfolder_content{file_index}]);
    if file_index==numel(subfolder_content) %&& (~isfield(trial,'state') || isempty(trial.state))
        break;
    end
    trial(file_index)=current_trial;
end

if exist('trial','var')
    save([folder filesep subfolder],'SETTINGS','trial');
end
rmdir([folder filesep subfolder],'s');

%%%
folder = 'C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\PerceptualM2S_Wagering\Data\M2S\behavior';
filename = 'ANKO28_02_2018-06-23_01.mat';
i = 1;
for file_index=99: 98+ numel(trial)
current_trial(file_index) = trial(i);
i = i+1;
end
trial = current_trial;
save([folder filesep filename],'SETTINGS','trial');
