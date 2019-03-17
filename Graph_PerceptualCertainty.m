clear all; close all; close all
%%What is the function doing?
%1) load the dataset created in the function
% CreateSaveVariables_PerceptualCertainty.m
%2) Graphs
pathname = 'Y:\Projects\Interoception\Perceptual_Certainty\Data\behavior'; 

subfolder_dir=dir([pathname filesep  '*.mat']); % checking only for mat files in the specified subfolder
subfolder_content={subfolder_dir.name};
subfolder_content=sort(subfolder_content);

pathname_Data = 'Y:\Projects\Interoception\Perceptual_Certainty\Results\DiscriminationPerformance_Certainty';
Files =dir([pathname_Data filesep  '*.mat']); % checking only for mat files in the specified subfolder

load([pathname_Data filesep Files.name])
%% Properties
cmap = jet(numel(Diff));
WritingLabelAxis_Size = 18;
%%
c = 1;
for file_index=   1 : numel(Diff)
% Plot performance per difficulty for the three subjects
figure(1); set(gcf,'Name','performance per difficulty');
subplot(1,2,1);
plot(  Diff(file_index).DiffLevels, Diff(file_index).Performance, '.-','color',cmap(file_index,:) ,'MarkerSize',25) 
set(gca,'YTick', 0.4:0.1:1);set(gca,'ylim',[0.4 1.05]); 
set(gca,'XTick', Diff(file_index).DiffLevels); set(gca,'xlim',[0 14]); 
ylabel('perceptual performance (%)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('degree of rotation (difficult to easy)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('performance per difficulty level','fontsize',14,'fontweight','b' );
Subject{c} = Diff(file_index).Subject; 
%legend(Subject);
hold on
% add the mean of all participants
%plot( Diff(file_index).DiffLevels), [Diff(:).Performance](1), 'k.-', ,'MarkerSize',35) 


subplot(1,2,2);
plot(  Diff(file_index).DiffLevels, Diff(file_index).NrTrial, '.-','color',cmap(file_index,:), 'MarkerSize',25)       
set(gca,'XTick', Diff(file_index).DiffLevels);set(gca,'xlim',[0 14]);
ylabel('Nr. trials','fontsize',12,'fontweight','b' );
xlabel('degree of rotation (difficult to easy)','fontsize',14,'fontweight','b' );
title('Nr. of trials per difficulty level','fontsize',14,'fontweight','b' );
hold on


figure(2); set(gcf,'Name','response Times per difficulty');
subplot(1,2,1);
%plot(  Diff.DiffLevels, Diff.RT_mean, Plotcolor{file_index}, 'MarkerSize',20) 
errorbar(Diff(file_index).DiffLevels,Diff(file_index).RT_mean,Diff(file_index).RT_sd,'.','color',cmap(file_index,:), 'MarkerSize',25)
set(gca,'YTick', 0.3:0.1:1.2);set(gca,'ylim',[0.25 1.2]); 
set(gca,'XTick', Diff(file_index).DiffLevels);set(gca,'xlim',[0 14]);
ylabel('reaction time  (s)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('degree of rotation (similar,difficult to easy)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('reaction time per difficulty level','fontsize',14,'fontweight','b' );
%legend(Subject);
hold on


subplot(1,2,2);
%plot(Diff.DiffLevels, Diff.MVT_mean, Plotcolor{file_index}, 'MarkerSize',20)   
errorbar(Diff(file_index).DiffLevels,Diff(file_index).MVT_mean,Diff(file_index).MVT_sd,'.','color',cmap(file_index,:), 'MarkerSize',25)
set(gca,'YTick', 0.3:0.1:1.2);set(gca,'ylim',[0.25 1.2]); 
set(gca,'XTick', Diff(file_index).DiffLevels);set(gca,'xlim',[0 14]);
ylabel('movement time (s)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('degree of rotation (similar,difficult to easy)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('movement time per difficulty level','fontsize',14,'fontweight','b' );
%legend(Subject); 
hold on

figure(3); set(gcf,'Name','performance per Sample');
subplot(2,1,1);
plot(  Sample(file_index).Samples, Sample(file_index).Performance, '.','color',cmap(file_index,:), 'MarkerSize',25)       
set(gca,'YTick', 0:0.1:1);set(gca,'ylim',[-0.05 1.05]); 
set(gca,'XTick', Sample(file_index).Samples(1):4:Sample(file_index).Samples(end));%set(gca,'xlim',[-1 21]); 
ylabel('perceptual performance (%)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('rotation of the Sample','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('performance per Sample presentd during M2S','fontsize',24,'fontweight','b' );
Subject{c} = Diff(file_index).Subject; 
%legend(Subject);
hold on

subplot(2,1,2);
plot(  Sample(file_index).Samples, Sample(file_index).NrTrial, '.','color',cmap(file_index,:), 'MarkerSize',25)  
set(gca,'XTick', Sample(file_index).Samples(1):4:Sample(file_index).Samples(end));
ylabel('Nr. tials','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('rotation of the Sample','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('Nr. of tials per Sample presentd during M2S','fontsize',14,'fontweight','b' );
hold on


figure(4); set(gcf,'Name','averageCertaintyRating (Post) vs perceptual Performance');
plot(  General(file_index).d, General(file_index).AverageCertainty, '.','color',cmap(file_index,:), 'MarkerSize',25)   
set(gca,'YTick', 0:1:6);set(gca,'ylim',[0 6]); 
set(gca,'XTick', 0:1:5);set(gca,'xlim',[0 5]); 
ylabel('averageCertaintyRating (Post)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('perceptual Performance (dprime)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
Subject{c} = Diff(file_index).Subject; 
text(General(file_index).d+0.05, General(file_index).AverageCertainty,Diff(file_index).Subject)

%legend(Subject);
hold on



figure(5); set(gcf,'Name','averagePerformance Post vs PreCertainty');
plot(  PostCertain(file_index).GeneralPerformance, PreCertain(file_index).GeneralPerformance, '.','color',cmap(file_index,:), 'MarkerSize',25)   
plot(0:5,0:5, 'k-')
set(gca,'YTick', 0:0.1:1);set(gca,'ylim',[0 1]); 
set(gca,'XTick', 0:0.1:1);set(gca,'xlim',[0 1]); 
ylabel('averagePerformance Pre-Certainty','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('averagePerformance Post-Certainty ','fontsize',WritingLabelAxis_Size,'fontweight','b' );
Subject{c} = Diff(file_index).Subject; 
text(PostCertain(file_index).GeneralPerformance+0.05, PreCertain(file_index).GeneralPerformance,Diff(file_index).Subject)

%legend(Subject);
hold on

figure(6); set(gcf,'Name','performance per Certainty Level (Post)');
subplot(1,2,1);
plot(  PostCertain(file_index).CertaintyLevels, PostCertain(file_index).Performance, '.-','color',cmap(file_index,:), 'MarkerSize',25)       
set(gca,'YTick', 0:0.1:1);set(gca,'ylim',[-0.05 1.05]); 
set(gca,'XTick', PostCertain(file_index).CertaintyLevels(1):1:PostCertain(file_index).CertaintyLevels(6));set(gca,'xlim',[0 7]); 
ylabel('performance (%)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('Certainty Levels ','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('performance per Post-Certainty Level','fontsize',14,'fontweight','b' );
Subject{c} = Diff(file_index).Subject; 
legend(Subject);hold on


subplot(1,2,2);
plot(  PostCertain(file_index).CertaintyLevels, PostCertain(file_index).NrTrial, '.','color',cmap(file_index,:), 'MarkerSize',25)       
set(gca,'XTick', PostCertain(file_index).CertaintyLevels(1):1:PostCertain(file_index).CertaintyLevels(6));set(gca,'xlim',[0 7]); 
ylabel('Nr. trials','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('Certainty Levels ','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('Nr. of trials per Post-Certainty Level','fontsize',14,'fontweight','b' );
hold on


figure(7); set(gcf,'Name','performance per Certainty Level (Pre)');
subplot(1,2,1);
plot(  PreCertain(file_index).CertaintyLevels, PreCertain(file_index).Performance, '.-','color',cmap(file_index,:), 'MarkerSize',25)       
set(gca,'YTick', 0:0.1:1);set(gca,'ylim',[-0.05 1.05]); 
set(gca,'XTick', PreCertain(file_index).CertaintyLevels(1):1:PreCertain(file_index).CertaintyLevels(6));set(gca,'xlim',[0 7]); 
ylabel('performance (%)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('Certainty Levels ','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('performance per Pre-Certainty Level','fontsize',14,'fontweight','b' );
Subject{c} = Diff(file_index).Subject; 
legend(Subject);hold on

subplot(1,2,2);
plot(  PreCertain(file_index).CertaintyLevels, PreCertain(file_index).NrTrial, '.-','color',cmap(file_index,:), 'MarkerSize',25)       
set(gca,'XTick', PreCertain(file_index).CertaintyLevels(1):1:PreCertain(file_index).CertaintyLevels(6));set(gca,'xlim',[0 7]); 
ylabel('Nr. trials','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('Certainty Levels ','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('Nr. of trials per Pre-Certainty Level','fontsize',14,'fontweight','b' );
hold on

figure(8); set(gcf,'Name','performance  vs dprime');% perceptual performance  & d'
plot( General(file_index).d, General(file_index).Performance, '.','color',cmap(file_index,:), 'MarkerSize',25)
set(gca,'YTick', 0:0.1:1);set(gca,'ylim',[0 1]); 
set(gca,'XTick', 0:1:5);set(gca,'xlim',[0 5]); 
ylabel('Performance (%)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('d (Type 1)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('','fontsize',14,'fontweight','b' );
Subject{c} = Diff(file_index).Subject; 
%legend(Subject);
hold on
text(General(file_index).d+0.05, General(file_index).Performance,Diff(file_index).Subject)


figure(9);  set(gcf,'Name','dprime  vs metaD');% Meta-D as function of D-prime
plot( General(file_index).d, General(file_index).metaD, '.','color',cmap(file_index,:), 'MarkerSize',25); hold on; 
set(gca,'YTick', -2:1:5);set(gca,'ylim',[-2 5]); 
set(gca,'XTick',-2:1:5);set(gca,'xlim',[-2 5]); 
ylabel('perceptual metacognition (meta d)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('d (Type 1)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('','fontsize',14,'fontweight','b' );
Subject{c} = Diff(file_index).Subject; 
text(General(file_index).d+0.05, General(file_index).metaD,Diff(file_index).Subject)
legend(Subject);hold on; 
if file_index == numel(subfolder_content)
    plot(-2:5,-2:5, 'k-'); 
plot([-2; 5].', [0; 0].', 'k--');
plot([0;0].', [-2; 5].', 'k--');
end




figure(10);  set(gcf,'Name','PreCertainty -linear fits');
subplot(2,1,1)
plot(PreCertain(file_index).CertaintyLevels, SlopeBasedMetacognition(file_index).r_conf_pre, '-','color',cmap(file_index,:), 'MarkerSize',25);hold on
plot(PreCertain(file_index).CertaintyLevels, PreCertain(file_index).PercSuc, '.--','color',cmap(file_index,:), 'MarkerSize',15)       
set(gca,'ylim',[-50 100]); 
set(gca,'XTick', 0:1:max(PreCertain(file_index).CertaintyLevels));set(gca,'xlim',[0 (max(PreCertain(file_index).CertaintyLevels) +1)]); 
ylabel('Percent of trials(correct)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('Wagers','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('PreCertainty Linear fits for correct trials','fontsize',14,'fontweight','b' );
Subject{c} = Diff(file_index).Subject; 
legend(Subject);hold on

subplot(2,1,2)
plot( PreCertain(file_index).CertaintyLevels, SlopeBasedMetacognition(file_index).r_err_pre, '-','color',cmap(file_index,:), 'MarkerSize',25)  
plot(PreCertain(file_index).CertaintyLevels, PreCertain(file_index).PercUnsuc, '.--','color',cmap(file_index,:), 'MarkerSize',15)  
set(gca,'ylim',[-50 100]); 
set(gca,'XTick', 0:1:max(PreCertain(file_index).CertaintyLevels));set(gca,'xlim',[0 (max(PreCertain(file_index).CertaintyLevels) +1)]); 
ylabel('Percent of trials(incorrect)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('Wagers','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('Linear fits for incorrect trials','fontsize',14,'fontweight','b' );
Subject{c} = Diff(file_index).Subject; 
legend(Subject);hold on

 
figure(11);  set(gcf,'Name','PostCertainty -linear fits');
subplot(2,1,1)
plot(PostCertain(file_index).CertaintyLevels,  SlopeBasedMetacognition(file_index).r_conf_post,  '-','color',cmap(file_index,:), 'MarkerSize',25)
plot(PostCertain(file_index).CertaintyLevels, PostCertain(file_index).PercSuc,  '.--','color',cmap(file_index,:), 'MarkerSize',15)  
set(gca,'ylim',[-50 100]); 
set(gca,'XTick', 0:1:max(PostCertain(file_index).CertaintyLevels));set(gca,'xlim',[0 (max(PostCertain(file_index).CertaintyLevels) +1)]); 
ylabel('Percent of trials (correct)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('Wagers','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('PostCertainty Linear fits for correct trials','fontsize',14,'fontweight','b' );
Subject{c} = Diff(file_index).Subject; 
legend(Subject);hold on

subplot(2,1,2)
plot( PostCertain(file_index).CertaintyLevels,  SlopeBasedMetacognition(file_index).r_err_post, '-','color',cmap(file_index,:), 'MarkerSize',25)  
plot(PostCertain(file_index).CertaintyLevels, PostCertain(file_index).PercUnsuc,  '.--','color',cmap(file_index,:), 'MarkerSize',15)
set(gca,'ylim',[-50 100]); 
set(gca,'XTick', 0:1:max(PostCertain(file_index).CertaintyLevels));set(gca,'xlim',[0 (max(PostCertain(file_index).CertaintyLevels) +1)]); 
ylabel('Percent of trials (incorrect)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
xlabel('Wagers','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('Linear fits for incorrect trials','fontsize',14,'fontweight','b' );
Subject{c} = Diff(file_index).Subject; 
legend(Subject);hold on

figure(12);  set(gcf,'Name','Bi-directional Certainty');
subplot(2,1,1)
plot(1:max(PostCertain(file_index).CertaintyLevels),  SlopeBasedMetacognition(file_index).r_conf_corrected,  '.-','color',cmap(file_index,:), 'MarkerSize',25)
set(gca,'XTick', 0:1:max(PostCertain(file_index).CertaintyLevels));set(gca,'xlim',[0 (max(PostCertain(file_index).CertaintyLevels) +1)]); 
ylabel('Percent of trials (correct)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('Bi-directional Certainty (PDW slopes minus PreDW slopes)','fontsize',14,'fontweight','b' );
legend(Subject);hold on
subplot(2,1,2)
plot(1:max(PostCertain(file_index).CertaintyLevels), SlopeBasedMetacognition(file_index).r_err_corrected,  '.-','color',cmap(file_index,:), 'MarkerSize',25)
set(gca,'XTick', 0:1:max(PostCertain(file_index).CertaintyLevels));set(gca,'xlim',[0 (max(PostCertain(file_index).CertaintyLevels) +1)]); 
ylabel('Percent of trials (incorrect)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
legend(Subject);hold on

figure(13);  set(gcf,'Name','Slope-Meta vs Meta-D');
plot(  General(file_index).metaD, General(file_index).meta_slope, '.-','color',cmap(file_index,:), 'MarkerSize',25)       
set(gca,'YTick', 0:1:8);set(gca,'ylim',[0 8]); 
set(gca,'XTick', 0:1:5);set(gca,'xlim',[0 5]); 
xlabel('perceptual metacognition (meta d)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
ylabel('slope-based metacognitive ability','fontsize',WritingLabelAxis_Size,'fontweight','b' );
title('','fontsize',14,'fontweight','b' );
text(General(file_index).metaD+0.03, General(file_index).meta_slope,Diff(file_index).Subject)
Subject{c} = Diff(file_index).Subject; 
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




%Diff        = []; 
Diff_Sample = []; 
end
% save mat-figures to file
path_SaveFig = ['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results']; 
path_SaveFig = ['Y:\Personal\Kristin\Talks\ScienceMeeting_Okt2018\Figures']; 

for i = 1:12
    h = [];
    h(1) = figure(i); 
    print(h,[path_SaveFig filesep 'Graphs_PerceptualCertainty_',num2str(i)], '-dpng')

end  

for i = 1:12
        h(i) = figure(i); 
end
savefig(h, [path_SaveFig filesep 'Graphs_PerceptualCertainty_102108.fig'])

%% create Tables to combine Data from all participants
%%% Create Table from GeneralStructure
TableGeneral = struct2table(General);
path_save = 'C:\Users\kkaduk\Dropbox\promotion\Projects\Metacogition_Interoception_Human\M2S_Certainty\Data\TableEachRowParticipant\';
save([path_save, 'Table_PercDiscTask_CertaintyScale_EachParticipantOneRow' ],'TableGeneral');
writetable(TableGeneral, [path_save, 'Table_PercDiscTask_CertaintyScale_EachParticipantOneRow'], 'Delimiter', ' ')


TableDiff = struct2table(Diff);


% 
% 
% % Plot performance per Sample for the three subjects
% 
% % Plot RT per Sample for the three subjects
% 
% 
% % ideas how to use other options as a for-loop
% groupedDifficulty = arrayfun(@(y)find(x == 0), unique(x), 'UniformOutput',false);
% Rotation_Sample = cell2mat(arrayfun(@(data) vertcat(data.hnd.cue(1).shape.rotation),data,'uni',0));
% diff(groupedDifficulty{1})
% 
% diff(success == 1)/ diff(success == 0);
% for i = 1: length(y)
%     groupedDifficulty(i) = sum(x == y(6)); 
% end