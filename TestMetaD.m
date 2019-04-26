pathname = 'Y:\Projects\Interoception\Perceptual_Certainty\Data\behavior'; 
subfolder_dir=dir([pathname filesep  '*.mat']); % checking only for mat files in the specified subfolder
subfolder_content={subfolder_dir.name};
subfolder_content=sort(subfolder_content);

file_index= 6 % 19 %6
for file_index=   1 : numel(subfolder_content)


load([pathname filesep subfolder_content{file_index}]);
selected_Left                        =[trial.selected_Left];
selected_Right                       =[trial.selected_Right];

Stimuli_LeftOrRight    =[trial.Stimuli_LeftOrRight];

unsuccess              =[trial.unsuccess];
success = zeros(1,length(selected_Left)); 
success((Stimuli_LeftOrRight == 2 & selected_Right == 1 | Stimuli_LeftOrRight == 1 & selected_Left == 1)) = 1; 

wagering_or_controll_wagering_post   =[trial.wagering_or_controll_wagering_post]; %postwagering =2
n_wagers                             =SETTINGS.NrOfWagers;
wager_choosen_post                   =[trial.wager_choosen_post];

%% General Performance
%% General Performace
completedTrials_success = success; 
completedTrials_success(unsuccess ~= 0) = NaN;
completedTrials_unsuccess = unsuccess; 
completedTrials_unsuccess(unsuccess ~= 0) = NaN;
CompletedTrials_postwagering = wager_choosen_post; 
CompletedTrials_postwagering(unsuccess ~= 0) = NaN;
CompletedTrials_postwagering(wagering_or_controll_wagering_post == 1) = NaN; %delete all Pre-Wagering


General(file_index).Subject      = subfolder_content{file_index}(1:6);  %HARD CODED!!!
General(file_index).NrTrial      = nansum(completedTrials_unsuccess == 0); 
General(file_index).Success      = nansum(completedTrials_success== 1) ; 
General(file_index).Unsuccess    = nansum(completedTrials_success == 0) ; 
General(file_index).Performance  = General(file_index).Success / (General(file_index).Success +General(file_index).Unsuccess) ;
General(file_index).AverageCertainty  = round(nansum(CompletedTrials_postwagering)/sum(~isnan(CompletedTrials_postwagering)),2) ;


%% Post-Wagering
Cer = 1:6 ;% unique(wager_choosen_post(unsuccess == 0));
CompletedTrials = wager_choosen_post; 
CompletedTrials(unsuccess ~= 0) = NaN;
CompletedTrials(wagering_or_controll_wagering_post == 1) = NaN; %delete all Pre-Wagering

for i = 1: length(Cer)
    PostCertain(file_index).Subject         = subfolder_content{file_index}(1:6); 
    PostCertain(file_index).CertaintyLevels(i)   = Cer(i); 
    PostCertain(file_index).NrTrial(i)      = sum(CompletedTrials == Cer(i)); 
    PostCertain(file_index).NrTrials_Suc(i)      = sum(success(find(CompletedTrials == Cer(i))) == 1) ; 
    PostCertain(file_index).NrTrials_Unsuc(i)    = sum(success(find(CompletedTrials == Cer(i))) == 0) ; 
    PostCertain(file_index).Performance(i)       = PostCertain(file_index).NrTrials_Suc(i)/ (PostCertain(file_index).NrTrials_Suc(i) +PostCertain(file_index).NrTrials_Unsuc(i)) ; 
    PostCertain(file_index).GeneralPerformance   = sum(PostCertain(file_index).NrTrials_Suc) /  sum(PostCertain(file_index).NrTrial);

end

PostCertain(file_index).PercSuc = PostCertain(file_index).NrTrials_Suc / sum(PostCertain(file_index).NrTrials_Suc)*100;
PostCertain(file_index).PercUnsuc = PostCertain(file_index).NrTrials_Unsuc / sum(PostCertain(file_index).NrTrials_Unsuc)*100;
PostCertain(file_index).PercAll = PostCertain(file_index).NrTrial / sum(PostCertain(file_index).NrTrial)*100;

%% META-D

out = metaD_PerSubject(n_wagers, wagering_or_controll_wagering_post, completedTrials_success,selected_Left, wager_choosen_post,selected_Right);
General(file_index).d       = out.da;
General(file_index).metaD   = out.meta_da;
General(file_index).metaEfficiency = out.M_ratio;
General(file_index).metaEfficiency_Difference = out.M_diff;
General(file_index).meta_c1     =  out.type2_fit.meta_c1; 
General(file_index).c_1         =  out.c_1; 
General(file_index).pHit         =  out.HR1; 
General(file_index).pFA         =  out.FAR1; 


% We get the same results for dprime with Igor's functioin using the hitrate and falsealarm 
testsim_dprime(General(file_index).pHit,General(file_index).pFA)

end