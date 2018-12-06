clear all; close all; 
c = 1; Subject = {}; 
pathname = 'C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Data\behavior'; 
cd('C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Analyze')

subfolder_dir=dir([pathname filesep  '*.mat']); % checking only for mat files in the specified subfolder
subfolder_content={subfolder_dir.name};
subfolder_content=sort(subfolder_content);
Plotcolor_D = [{'g.' };{'m.'};{'b.' }; {'c.' }; {'r.' }; {'k.' }];
Plotcolor_L = [{'g.-' };{'m.-'};{'b.-' }; {'c.-' };{'r.-' }; {'k.-' }];
Plotcolor_Line = [{'g-' };{'m-'};{'b-' }; {'c-' };{'r-' }; {'k' }];
Plotcolor_LinedDots = [{'g.--' };{'m.--'};{'b.--' }; {'c.--' };{'r.--' }; {'k.--' }];

Plotcolor_Star = [{'g*' };{'m*'};{'b*' }; {'c*' }; {'r*' }; {'k*' }];

for file_index= 4 1 : numel(subfolder_content)
    %intiate Variables
load([pathname filesep subfolder_content{file_index}]);
unsuccess             =[trial.unsuccess];
refAngle_Sample       =[trial.refAngle_Sample];
refAngle_M2S          =[trial.refAngle_M2S];
Stimuli_LeftOrRight =[trial.Stimuli_LeftOrRight];
selected_Right         =[trial.selected_Right];
selected_Left         =[trial.selected_Left];
Time_matchtosample_appeared_run                     = [trial.Time_matchtosample_appeared_run];
Time_match_to_sample_button_released_run            = [trial.Time_match_to_sample_button_released_run];
Time_match_to_sample_selected_1stPress_run          = [trial.Time_match_to_sample_selected_1stPress_run];
Time_match_to_sample_selected_run                   = [trial.Time_match_to_sample_selected_run];
Time_match_to_sample_ReactionTime = Time_match_to_sample_button_released_run - Time_matchtosample_appeared_run;
Time_match_to_sample_MovementTime = Time_match_to_sample_selected_run - Time_matchtosample_appeared_run;
position_wager_chosen_pre   =[trial.position_wager_chosen_pre];
position_wager_chosen_post  =[trial.position_wager_chosen_post];
wager_choosen_pre           =[trial.wager_choosen_pre];
wager_choosen_post          =[trial.wager_choosen_post];
wagering_or_controll_wagering_post   =[trial.wagering_or_controll_wagering_post]; %postwagering =2
wagering_or_controll_wagering_pre    =[trial.wagering_or_controll_wagering_pre];%

n_wagers        =SETTINGS.NrOfWagers;
n_trials        =SETTINGS.n_trials;
control_post    =[trial.wagering_or_controll_wagering_post];
control_pre     =[trial.wagering_or_controll_wagering_pre];
selected_Right  =[trial.selected_Right];
Stimuli_LeftOrRight     =[trial.Stimuli_LeftOrRight];
selected_Left           =[trial.selected_Left];
wager_choosen_post      =[trial.wager_choosen_post];
wager_choosen_pre       =[trial.wager_choosen_pre];
correctSelection_M2S    =[trial.correctSelection_M2S];

success = zeros(1,length(Time_match_to_sample_MovementTime)); 
success((Stimuli_LeftOrRight == 2 & selected_Right == 1 | Stimuli_LeftOrRight == 1 & selected_Left == 1)) = 1; 
nansum(success)/length(success);
%ANHA30_01_02_2018-06-18.mat , in total: 311 trials, 102 trials with
%fixation break etc mistakes; 83 successful; in total: 40% success rate
diff = abs(refAngle_M2S - refAngle_Sample ); 

%% General Performace
completedTrials_success = success; 
completedTrials_success(unsuccess ~= 0) = NaN;
completedTrials_unsuccess = unsuccess; 
completedTrials_unsuccess(unsuccess ~= 0) = NaN;
CompletedTrials_postwagering = wager_choosen_post; 
CompletedTrials_postwagering(unsuccess ~= 0) = NaN;
CompletedTrials_postwagering(wagering_or_controll_wagering_post == 1) = NaN; %delete all Pre-Wagering

%unsuccess(unsuccess ~= 0) = [];
General(file_index).Subject      = subfolder_content{file_index}(1:6);  %HARD CODED!!!
General(file_index).NrTrial      = nansum(completedTrials_unsuccess == 0); 
General(file_index).Success      = nansum(completedTrials_success== 1) ; 
General(file_index).Unsuccess    = nansum(completedTrials_success == 0) ; 
General(file_index).Performance  = General(file_index).Success / (General(file_index).Success +General(file_index).Unsuccess) ;
General(file_index).AverageCertainty  = round(nansum(CompletedTrials_postwagering)/sum(~isnan(CompletedTrials_postwagering)),2) ;



%% Difficulty Level
% How many trials per difficulty? Performance per Difficulty? 
Diff_completedTrials = diff;
Diff_completedTrials(unsuccess ~= 0) = NaN;
y = unique(Diff_completedTrials(unsuccess == 0));
for i = 1: length(y)
    Diff.Subject         = subfolder_content{file_index}(1:6);  %HARD CODED!!!
    Diff.DiffLevels(i)   = y(i); 
    Diff.NrTrial(i)      = nansum(Diff_completedTrials == y(i)); 
    Diff.Success(i)      = nansum(success(find(Diff_completedTrials == y(i))) == 1) ; 
    Diff.Unsuccess(i)    = nansum(success(find(Diff_completedTrials == y(i))) == 0) ; 
    Diff.Performance(i)  = Diff.Success(i)/ (Diff.Success(i) +Diff.Unsuccess(i)) ; 
    Diff.RT_mean(i)      = nanmean(Time_match_to_sample_ReactionTime(find(Diff_completedTrials == y(i)))); 
    Diff.RT_sd(i)        = nanstd(Time_match_to_sample_ReactionTime(find(Diff_completedTrials == y(i))));     
    Diff.MVT_mean(i)     = nanmean(Time_match_to_sample_MovementTime(find(Diff_completedTrials == y(i)))); 
    Diff.MVT_sd(i)       = nanstd(Time_match_to_sample_MovementTime(find(Diff_completedTrials == y(i))));
    unique(refAngle_Sample(find(Diff_completedTrials == y(i)))) ;
    unique(refAngle_M2S(find(Diff_completedTrials == y(i)))) ;

    [refAngle_Sample(find(Diff_completedTrials == y(i))); refAngle_M2S(find(Diff_completedTrials == y(i)))];
    
end

%% Sample
Samples = unique(refAngle_Sample(unsuccess == 0));
CompletedTrials = refAngle_Sample; 
CompletedTrials(unsuccess ~= 0) = NaN;

for i = 1: length(Samples)
    Sample.Subject         = subfolder_content{file_index}(1:6); 
    Sample.Samples(i)   = Samples(i); 
    Sample.NrTrial(i)      = nansum(CompletedTrials == Samples(i)); 
    Sample.Success(i)      = nansum(success(find(CompletedTrials == Samples(i))) == 1) ; 
    Sample.Unsuccess(i)    = nansum(success(find(CompletedTrials == Samples(i))) == 0) ; 
    Sample.Performance(i)  = Sample.Success(i)/ (Sample.Success(i) +Sample.Unsuccess(i)) ; 
    Sample.RT_mean(i)      = nanmean(Time_match_to_sample_ReactionTime(find(CompletedTrials == Samples(i)))); 
    Sample.RT_sd(i)        = nanstd(Time_match_to_sample_ReactionTime(find(CompletedTrials == Samples(i))));     
    Sample.MVT_mean(i)     = nanmean(Time_match_to_sample_MovementTime(find(CompletedTrials == Samples(i)))); 
    Sample.MVT_sd(i)       = nanstd(Time_match_to_sample_MovementTime(find(CompletedTrials == Samples(i))));
end

%% Samples per Difficulty Level
Diff_completedTrials = diff;
Diff_completedTrials(unsuccess ~= 0) = NaN;
Dif= unique(Diff_completedTrials(unsuccess == 0));
for i_diff = 1: length(Dif)
    Diff_Ind = find(Diff_completedTrials == Dif(i_diff));
    Samples = unique(refAngle_Sample(Diff_Ind));
    CompletedTrials = refAngle_Sample(Diff_Ind);
for i = 1: length(Samples)
    Diff_Sample(i_diff).Subject         = subfolder_content{file_index}(1:6); 
    Diff_Sample(i_diff).Samples(i)   = Samples(i); 
    Diff_Sample(i_diff).NrTrial(i)      = nansum(CompletedTrials == Samples(i)); 
    Diff_Sample(i_diff).Success(i)      = nansum(success(find(CompletedTrials == Samples(i))) == 1) ; 
    Diff_Sample(i_diff).Unsuccess(i)    = nansum(success(find(CompletedTrials == Samples(i))) == 0) ; 
    Diff_Sample(i_diff).Performance(i)  = Diff_Sample(i_diff).Success(i)/ (Diff_Sample(i_diff).Success(i) +Diff_Sample(i_diff).Unsuccess(i)) ; 
    Diff_Sample(i_diff).RT_mean(i)      = nanmean(Time_match_to_sample_ReactionTime(find(CompletedTrials == Samples(i)))); 
    Diff_Sample(i_diff).RT_sd(i)        = nanstd(Time_match_to_sample_ReactionTime(find(CompletedTrials == Samples(i))));     
    Diff_Sample(i_diff).MVT_mean(i)     = nanmean(Time_match_to_sample_MovementTime(find(CompletedTrials == Samples(i)))); 
    Diff_Sample(i_diff).MVT_sd(i)       = nanstd(Time_match_to_sample_MovementTime(find(CompletedTrials == Samples(i))));
end
end

%% Post-Wagering
Cer = 1:6 ;% unique(wager_choosen_post(unsuccess == 0));
CompletedTrials = wager_choosen_post; 
CompletedTrials(unsuccess ~= 0) = NaN;
CompletedTrials(wagering_or_controll_wagering_post == 1) = NaN; %delete all Pre-Wagering

for i = 1: length(Cer)
    PostCertain.Subject         = subfolder_content{file_index}(1:6); 
    PostCertain.CertaintyLevels(i)   = Cer(i); 
    PostCertain.NrTrial(i)      = sum(CompletedTrials == Cer(i)); 
    PostCertain.NrTrials_Suc(i)      = sum(success(find(CompletedTrials == Cer(i))) == 1) ; 
    PostCertain.NrTrials_Unsuc(i)    = sum(success(find(CompletedTrials == Cer(i))) == 0) ; 
    PostCertain.Performance(i)       = PostCertain.NrTrials_Suc(i)/ (PostCertain.NrTrials_Suc(i) +PostCertain.NrTrials_Unsuc(i)) ; 
    PostCertain.GeneralPerformance   = sum(PostCertain.NrTrials_Suc) /  sum(PostCertain.NrTrial);

end

PostCertain.PercSuc = PostCertain.NrTrials_Suc / sum(PostCertain.NrTrials_Suc)*100;
PostCertain.PercUnsuc = PostCertain.NrTrials_Unsuc / sum(PostCertain.NrTrials_Unsuc)*100;
PostCertain.PercAll = PostCertain.NrTrial / sum(PostCertain.NrTrial)*100;


Cer = 1:6 ;%unique(wager_choosen_pre(unsuccess == 0));
CompletedTrials = wager_choosen_pre; 
CompletedTrials(unsuccess ~= 0) = NaN;
CompletedTrials(wagering_or_controll_wagering_post == 2) = NaN;

for i = 1: length(Cer)
    PreCertain.Subject         = subfolder_content{file_index}(1:6); 
    PreCertain.CertaintyLevels(i)   = Cer(i); 
    PreCertain.NrTrial(i)      = sum(CompletedTrials == Cer(i)); 
    PreCertain.NrTrials_Suc(i)      = sum(success(find(CompletedTrials == Cer(i))) == 1) ; 
    PreCertain.NrTrials_Unsuc(i)    = sum(success(find(CompletedTrials == Cer(i))) == 0) ; 
    PreCertain.Performance(i)  = PreCertain.NrTrials_Suc(i)/ (PreCertain.NrTrials_Suc(i) +PreCertain.NrTrials_Unsuc(i)) ; 
    PreCertain.GeneralPerformance = sum(PreCertain.NrTrials_Suc) /  sum(PreCertain.NrTrial);
end

PreCertain.PercSuc = PreCertain.NrTrials_Suc / sum(PreCertain.NrTrials_Suc)*100;
PreCertain.PercUnsuc = PreCertain.NrTrials_Unsuc / sum(PreCertain.NrTrials_Unsuc)*100;
PreCertain.PercAll = PreCertain.NrTrial / sum(PreCertain.NrTrial)*100;


out = metaD_PerSubject(n_wagers, wagering_or_controll_wagering_post, completedTrials_success,selected_Left, wager_choosen_post,selected_Right);
General(file_index).d       = out.da;
General(file_index).metaD   = out.meta_da;
General(file_index).metaEfficiency = out.M_ratio;
General(file_index).Typ1_criterion = out.Typ1_criterion;
%% SLOPE-BASED METACOGNITION
% linear fit of all correct trials as function of Wagers of one subject
        p_conf_post      = polyfit(1:n_wagers,PostCertain.PercSuc,1);
		r_conf_post      = p_conf_post(1) .* [1:n_wagers] + p_conf_post(2);
		p_err_post       = polyfit(1:n_wagers,PostCertain.PercUnsuc,1);
		r_err_post       = p_err_post(1) .* [1:n_wagers] + p_err_post(2);
		
        p_conf_pre       = polyfit(1:n_wagers,PreCertain.PercSuc,1);
		r_conf_pre       = p_conf_pre(1) .* [1:n_wagers] + p_conf_pre(2);
		p_err_pre        = polyfit(1:n_wagers,PreCertain.PercUnsuc,1);
		r_err_pre        = p_err_pre(1) .* [1:n_wagers] + p_err_pre(2);
        r_conf_corrected = r_conf_post(:) - r_conf_pre(:);
		r_err_corrected  = r_err_post(:)  - r_err_pre(:);
		
		meta_conf          = p_conf_post(1) - p_conf_pre(1);
		meta_err           = -p_err_post(1) + p_err_pre(1);
		meta_conf_pre      = p_conf_post(1) - p_conf_pre(1);
		meta_err_pre       = -p_err_post(1) + p_err_pre(1);
		General(file_index).meta_slope               = meta_conf(:) + meta_err(:);
		meta_pre           = meta_conf_pre + meta_err_pre;
        
%%
WritingLabelAxis_Size = 18;
% Plot performance per difficulty for the three subjects
figure(1); set(gcf,'Name','performance per difficulty');
subplot(1,2,1);
plot(  Diff.DiffLevels, Diff.Performance, Plotcolor_L{file_index}, 'MarkerSize',25) 
set(gca,'YTick', 0.4:0.1:1);set(gca,'ylim',[0.4 1.05]); 
set(gca,'XTick', Diff.DiffLevels); set(gca,'xlim',[0 14]); 
ylabel('perceptual performance (%)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('degree of rotation (difficult to easy)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('performance per difficulty level','fontsize',14,'fontweight','b' );
Subject{c} = Diff.Subject; 
legend(Subject);hold on

subplot(1,2,2);
plot(  Diff.DiffLevels, Diff.NrTrial, Plotcolor_L{file_index}, 'MarkerSize',25)       
set(gca,'XTick', Diff.DiffLevels);set(gca,'xlim',[0 14]);
ylabel('Nr. trials','fontsize',12,'fontweight','b' );
xlabel('degree of rotation (difficult to easy)','fontsize',14,'fontweight','b' );
title('Nr. of trials per difficulty level','fontsize',14,'fontweight','b' );
hold on


figure(2); set(gcf,'Name','response Times per difficulty');
subplot(1,2,1);
%plot(  Diff.DiffLevels, Diff.RT_mean, Plotcolor{file_index}, 'MarkerSize',20) 
errorbar(Diff.DiffLevels,Diff.RT_mean,Diff.RT_sd,Plotcolor_D{file_index}, 'MarkerSize',25)
set(gca,'YTick', 0.3:0.1:1.2);set(gca,'ylim',[0.25 1.2]); 
set(gca,'XTick', Diff.DiffLevels);set(gca,'xlim',[0 14]);
ylabel('reaction time  (s)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('degree of rotation (similar,difficult to easy)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('reaction time per difficulty level','fontsize',14,'fontweight','b' );
legend(Subject);hold on


subplot(1,2,2);
%plot(Diff.DiffLevels, Diff.MVT_mean, Plotcolor{file_index}, 'MarkerSize',20)   
errorbar(Diff.DiffLevels,Diff.MVT_mean,Diff.MVT_sd,Plotcolor_D{file_index}, 'MarkerSize',25)
set(gca,'YTick', 0.3:0.1:1.2);set(gca,'ylim',[0.25 1.2]); 
set(gca,'XTick', Diff.DiffLevels);set(gca,'xlim',[0 14]);
ylabel('movement time (s)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('degree of rotation (similar,difficult to easy)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('movement time per difficulty level','fontsize',14,'fontweight','b' );
legend(Subject); hold on

figure(3); set(gcf,'Name','performance per Sample');
subplot(2,1,1);
plot(  Sample.Samples, Sample.Performance, Plotcolor_L{file_index}, 'MarkerSize',25)       
set(gca,'YTick', 0:0.1:1);set(gca,'ylim',[-0.05 1.05]); 
set(gca,'XTick', Sample.Samples(1):4:Sample.Samples(end));%set(gca,'xlim',[-1 21]); 
ylabel('perceptual performance (%)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('rotation of the Sample','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('performance per Sample presentd during M2S','fontsize',24,'fontweight','b' );
Subject{c} = Diff.Subject; 
legend(Subject);hold on

subplot(2,1,2);
plot(  Sample.Samples, Sample.NrTrial, Plotcolor_L{file_index}, 'MarkerSize',25)  
set(gca,'XTick', Sample.Samples(1):4:Sample.Samples(end));
ylabel('Nr. tials','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('rotation of the Sample','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('Nr. of tials per Sample presentd during M2S','fontsize',14,'fontweight','b' );
hold on


figure(4); set(gcf,'Name','averageCertaintyRating (Post) vs perceptual Performance');
plot(  out.da, General(file_index).AverageCertainty, Plotcolor_D{file_index}, 'MarkerSize',25)   
set(gca,'YTick', 0:1:6);set(gca,'ylim',[0 6]); 
set(gca,'XTick', 0:1:5);set(gca,'xlim',[0 5]); 
ylabel('averageCertaintyRating (Post)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('perceptual Performance (dprime)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
Subject{c} = Diff.Subject; 
legend(Subject);hold on



figure(5); set(gcf,'Name','averagePerformance Post vs PreCertainty');
plot(  PostCertain.GeneralPerformance, PreCertain.GeneralPerformance, Plotcolor_D{file_index}, 'MarkerSize',25)   
plot(0:5,0:5, 'k-')
set(gca,'YTick', 0:0.1:1);set(gca,'ylim',[0 1]); 
set(gca,'XTick', 0:0.1:1);set(gca,'xlim',[0 1]); 
ylabel('averagePerformance Pre-Certainty','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('averagePerformance Post-Certainty ','fontsize',WritingLabelAxis_Size,'fontweight','b' );
Subject{c} = Diff.Subject; 
legend(Subject);hold on

figure(6); set(gcf,'Name','performance per Certainty Level (Post)');
subplot(1,2,1);
plot(  PostCertain.CertaintyLevels, PostCertain.Performance, Plotcolor_L{file_index}, 'MarkerSize',25)       
set(gca,'YTick', 0:0.1:1);set(gca,'ylim',[-0.05 1.05]); 
set(gca,'XTick', PostCertain.CertaintyLevels(1):1:PostCertain.CertaintyLevels(6));set(gca,'xlim',[0 7]); 
ylabel('performance (%)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('Certainty Levels ','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('performance per Post-Certainty Level','fontsize',14,'fontweight','b' );
Subject{c} = Diff.Subject; 
legend(Subject);hold on


subplot(1,2,2);
plot(  PostCertain.CertaintyLevels, PostCertain.NrTrial, Plotcolor_L{file_index}, 'MarkerSize',25)       
set(gca,'XTick', PostCertain.CertaintyLevels(1):1:PostCertain.CertaintyLevels(6));set(gca,'xlim',[0 7]); 
ylabel('Nr. trials','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('Certainty Levels ','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('Nr. of trials per Post-Certainty Level','fontsize',14,'fontweight','b' );
hold on


figure(7); set(gcf,'Name','performance per Certainty Level (Pre)');
subplot(1,2,1);
plot(  PreCertain.CertaintyLevels, PreCertain.Performance, Plotcolor_L{file_index}, 'MarkerSize',25)       
set(gca,'YTick', 0:0.1:1);set(gca,'ylim',[-0.05 1.05]); 
set(gca,'XTick', PreCertain.CertaintyLevels(1):1:PreCertain.CertaintyLevels(6));set(gca,'xlim',[0 7]); 
ylabel('performance (%)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('Certainty Levels ','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('performance per Pre-Certainty Level','fontsize',14,'fontweight','b' );
Subject{c} = Diff.Subject; 
legend(Subject);hold on

subplot(1,2,2);
plot(  PreCertain.CertaintyLevels, PreCertain.NrTrial, Plotcolor_L{file_index}, 'MarkerSize',25)       
set(gca,'XTick', PreCertain.CertaintyLevels(1):1:PreCertain.CertaintyLevels(6));set(gca,'xlim',[0 7]); 
ylabel('Nr. trials','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('Certainty Levels ','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('Nr. of trials per Pre-Certainty Level','fontsize',14,'fontweight','b' );
hold on

figure(8)% perceptual performance  & d'
plot( out.da, General(file_index).Performance, Plotcolor_D{file_index}, 'MarkerSize',25)
set(gca,'YTick', 0:0.1:1);set(gca,'ylim',[0 1]); 
set(gca,'XTick', 0:1:5);set(gca,'xlim',[0 5]); 
ylabel('Performance (%)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('d (Type 1)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('','fontsize',14,'fontweight','b' );
Subject{c} = Diff.Subject; 
legend(Subject);hold on


figure(9)% Meta-D as function of D-prime
plot( out.da, out.meta_da, Plotcolor_D{file_index}, 'MarkerSize',25);hold on
plot(0:5,0:5, 'k-')
set(gca,'YTick', 0:1:5);set(gca,'ylim',[0 5]); 
set(gca,'XTick', 0:1:5);set(gca,'xlim',[0 5]); 
ylabel('perceptual metacognition (meta d)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('d (Type 1)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('','fontsize',14,'fontweight','b' );
Subject{c} = Diff.Subject; 
legend(Subject);hold on



figure(10)%
subplot(2,1,1)
plot(1:n_wagers, r_conf_pre,  Plotcolor_Line{file_index}, 'MarkerSize',25);hold on
plot(1:n_wagers, PreCertain.PercSuc,  Plotcolor_LinedDots{file_index}, 'MarkerSize',15)       
set(gca,'ylim',[-50 100]); 
set(gca,'XTick', 0:1:n_wagers);set(gca,'xlim',[0 (n_wagers +1)]); 
ylabel('Percent of trials(correct)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('Wagers','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('PreCertainty Linear fits for correct trials','fontsize',14,'fontweight','b' );
Subject{c} = Diff.Subject; 
legend(Subject);hold on

subplot(2,1,2)
plot( 1:n_wagers, r_err_pre, Plotcolor_Line{file_index}, 'MarkerSize',25)  
plot(1:n_wagers, PreCertain.PercUnsuc,  Plotcolor_LinedDots{file_index}, 'MarkerSize',15)  
set(gca,'ylim',[-50 100]); 
set(gca,'XTick', 0:1:n_wagers);set(gca,'xlim',[0 (n_wagers +1)]); 
ylabel('Percent of trials(incorrect)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('Wagers','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('Linear fits for incorrect trials','fontsize',14,'fontweight','b' );
Subject{c} = Diff.Subject; 
legend(Subject);hold on

 
figure(11)% Meta-D as function of D-prime
subplot(2,1,1)
plot(1:n_wagers, r_conf_post,  Plotcolor_Line{file_index}, 'MarkerSize',25)
plot(1:n_wagers, PostCertain.PercSuc,  Plotcolor_LinedDots{file_index}, 'MarkerSize',15)  
set(gca,'ylim',[-50 100]); 
set(gca,'XTick', 0:1:n_wagers);set(gca,'xlim',[0 (n_wagers +1)]); 
ylabel('Percent of trials (correct)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('Wagers','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('PostCertainty Linear fits for correct trials','fontsize',14,'fontweight','b' );
Subject{c} = Diff.Subject; 
legend(Subject);hold on

subplot(2,1,2)
plot( 1:n_wagers, r_err_post, Plotcolor_Line{file_index}, 'MarkerSize',25)  
plot(1:n_wagers, PostCertain.PercUnsuc,  Plotcolor_LinedDots{file_index}, 'MarkerSize',15)
set(gca,'ylim',[-50 100]); 
set(gca,'XTick', 0:1:n_wagers);set(gca,'xlim',[0 (n_wagers +1)]); 
ylabel('Percent of trials (incorrect)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('Wagers','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('Linear fits for incorrect trials','fontsize',14,'fontweight','b' );
Subject{c} = Diff.Subject; 
legend(Subject);hold on

figure(12)% Meta-D as function of D-prime
subplot(2,1,1)
plot(1:n_wagers, r_conf_corrected,  Plotcolor_L{file_index}, 'MarkerSize',25)
set(gca,'XTick', 0:1:n_wagers);set(gca,'xlim',[0 (n_wagers +1)]); 
ylabel('Percent of trials (correct)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('Bi-directional Certainty (PDW slopes minus PreDW slopes)','fontsize',14,'fontweight','b' );
legend(Subject);hold on
subplot(2,1,2)
plot(1:n_wagers, r_err_corrected,  Plotcolor_L{file_index}, 'MarkerSize',25)
set(gca,'XTick', 0:1:n_wagers);set(gca,'xlim',[0 (n_wagers +1)]); 
ylabel('Percent of trials (incorrect)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
legend(Subject);hold on

figure(13)% Slope-based metacognitive ability & Meta-D
plot( out.meta_da, General(file_index).meta_slope, Plotcolor_L{file_index}, 'MarkerSize',25)       
%set(gca,'YTick', 0:1:5);set(gca,'ylim',[0 5]); 
set(gca,'XTick', 0:1:5);set(gca,'xlim',[0 5]); 
xlabel('perceptual metacognition (meta d)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
ylabel('slope-based metacognitive ability','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('','fontsize',14,'fontweight','b' );
Subject{c} = Diff.Subject; 
legend(Subject);hold on

% 
% 
% figure(6)
% for i_diff = 1: length(Dif)
% ax1 = subplot(2,3,i_diff);
% plot( Diff_Sample(i_diff).Samples, Diff_Sample(i_diff).Performance, Plotcolor_L{file_index}, 'MarkerSize',25)       
% set(gca,'YTick', 0:0.1:1);set(gca,'ylim',[-0.05 1.05]); 
% set(gca,'XTick', Diff_Sample(i_diff).Samples(1):4:Diff_Sample(i_diff).Samples(end));%set(gca,'xlim',[-1 21]); 
% ylabel('performance (%)','fontsize',12,'fontweight','b' );
% xlabel('rotation of the Sample','fontsize',12,'fontweight','b' );
% title(['DifficultyLevel  ' num2str(Dif(i_diff)) ' °'],'fontsize',14,'fontweight','b' );
% Subject{c} = Diff.Subject; 
% legend(Subject);hold on
% end


% 
% figure(7)
% for i_diff = 1: length(Dif)
% ax1 = subplot(2,3,i_diff);
% plot( Diff_Sample(i_diff).Samples, Diff_Sample(i_diff).NrTrial, Plotcolor_L{file_index}, 'MarkerSize',25)       
% set(gca,'YTick', 0:1:10);set(gca,'ylim',[-0.05 10]); 
% set(gca,'XTick', Diff_Sample(i_diff).Samples(1):4:Diff_Sample(i_diff).Samples(end));%set(gca,'xlim',[-1 21]); 
% ylabel('Nr. of Trials','fontsize',12,'fontweight','b' );
% xlabel('rotation of the Sample','fontsize',12,'fontweight','b' );
% title(['DifficultyLevel  ' num2str(Dif(i_diff)) ' °'],'fontsize',14,'fontweight','b' );
% Subject{c} = Diff.Subject; 
% legend(Subject);hold on
% end
c = c +1; 

%% create Tables to combine Data from all participants
General


PreCertain  = [];
PostCertain = []; 
Diff        = []; 
Sample      = []; 
Diff_Sample = []; 
end
% save mat-figures to file
path_SaveFig = ['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results\']; 
path_SaveFig = ['Y:\Personal\Kristin\Talks\ScienceMeeting_Okt2018\Figures\']; 

for i = 1:12
    h(i) = figure(i); 
end  
savefig(h, [path_SaveFig  'Graphs_PerceptualCertainty_102108.fig'])
save([path_SaveFig  'Graphs_PerceptualCertainty_102108'], 'png')

print(h,[path_SaveFig  'Graphs_PerceptualCertainty_102108'], '-dpng')


%%% Create Table from GeneralStructure
TableGeneral = struct2table(General);
path_save = 'C:\Users\kkaduk\Dropbox\promotion\Projects\Metacogition_Interoception_Human\M2S_Certainty\Data\TableEachRowParticipant\';
save([path_save, 'Table_PercDiscTask_CertaintyScale_EachParticipantOneRow' ],'TableGeneral');
writetable(TableGeneral, [path_save, 'Table_PercDiscTask_CertaintyScale_EachParticipantOneRow'], 'Delimiter', ' ')






% Plot performance per Sample for the three subjects

% Plot RT per Sample for the three subjects


% ideas how to use other options as a for-loop
groupedDifficulty = arrayfun(@(y)find(x == 0), unique(x), 'UniformOutput',false);
Rotation_Sample = cell2mat(arrayfun(@(data) vertcat(data.hnd.cue(1).shape.rotation),data,'uni',0));
diff(groupedDifficulty{1})

diff(success == 1)/ diff(success == 0);
for i = 1: length(y)
    groupedDifficulty(i) = sum(x == y(6)); 
end