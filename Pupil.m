% Main Script using the pupil-size-master toolbox
%1. create a mat-filefor each trial with the data-structure needed for the
%toolbox
%2.

%% TODO
% 1. check which eye is better
% 2. very long gaps cut out
%% Housekeeping:
clc; clear;
%%
WritingLabelAxis_Size = 20;
preprocessing.AddPath               = 0;
preprocessing.CreateMatFiles        = 1;
preprocessing.compare2Procedures    = 0;
preprocessing.analysis = 1;
preprocessing.GraphsPreprocessing = 0;

path_beh = 'C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Data\behavior\'; 
%%
if preprocessing.AddPath
    % Add paths:
    addpath(genpath('..\..\helperFunctions\'));
    addpath('..\..\dataModels\');
    
    % Check if the code needs to run in legacy mode (if MATLAB is older than
    % v2013b and the 'table' datatype has not yet been introduced, a folder
    % containing a table spoofer is added to the path, see
    % ..\..\helperFunctions\LEGACY\):
    legacyModeCheck();
end
%% prepare the data recorded in the variable as mat-file
if preprocessing.CreateMatFiles
    path_beh = 'C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Data\behavior\';
    files = dir([path_beh ,'*.mat']);
    for fileIndx = 2 :numel(files)
        Data = load([path, files(fileIndx).name]);
        succesfull_Trials = find([Data.trial.unsuccess]== 0);
        mkdir(['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results\Pupil\matfiles\',files(fileIndx).name(1:6)] )
        
        for TrialIndex =   1: length(succesfull_Trials) ; %TrialIndex = 65
            TrialIndx = succesfull_Trials(TrialIndex);
            matFilename = [files(fileIndx).name(1:end-4),'_Trial_',num2str(TrialIndx) ,'.mat'];
            %raw data
            L_raw           = Data.trial(TrialIndx).save_Pupil_w(:,1);
            R_raw           = Data.trial(TrialIndx).save_Pupil_w(:,2);
            t_ms            = Data.trial(TrialIndx).save_Eye_Time(:,1)*1000; %convert from sec to ms
            
            % gab in the beginning of the raw data because of breaks ->
            % usually 8-20s difference
            %Difference = (vel(t_ms*1000));
            % t_ms(find(Difference > 60000))'
            %median(Difference)
            % remove NaNs - specifically, there are a lot on the end of the vector from initiating the vector
            L_raw = L_raw(~isnan(L_raw)); %delete NaNs
            R_raw = R_raw(~isnan(R_raw)); %delete NaNs
            t_ms =  t_ms(~isnan(t_ms));
            
            % Right eye was not recorded & left eye has less variance
            if all(R_raw == intmin('int16')) || nanstd(Data.trial(TrialIndx).save_Pupil_w(:,1))< nanstd(Data.trial(TrialIndx).save_Pupil_w(:,2))
                L = L_raw;
                curEye = 'left';
                Eye = 1;
            else
                L = R_raw;
                curEye = 'right';
                Eye = 2;
            end
            
            % Remove the 0 samples:
            L(L==0) = NaN;
            
            %each trial starts from zero defined by its own first sample
            %(using zeroTime_ms) or by start of the trial
            zeroTime_ms = t_ms(1);
            t_ms_backup = t_ms;
            t_ms        = (t_ms - (Data.trial(TrialIndx).Time_pressed_rest_Fix1Base_run*1000)); %Data.trial(TrialIndx).Trial_start_time_run
            
            diameterUnit = 'px';
            diameter     = struct('t_ms',t_ms,'L',L,'R',[]);
            
            %plot raw data
            if preprocessing.GraphsPreprocessing ;
                
                p0 = figure(2);
                plot(t_ms, L_raw,'.');
                ylabel('Pupil Dilation ','fontsize',WritingLabelAxis_Size,'fontweight','b' );
                xlabel('Time (ms)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
                onset_Ind_TrialStart = Data.trial(TrialIndx).Time_pressed_rest_Fix1Base_run*1000  - Data.trial(TrialIndx).Time_pressed_rest_Fix1Base_run*1000 ;
                onset_IndSample         = Data.trial(TrialIndx).Time_sample_appeared_run*1000   -Data.trial(TrialIndx).Time_pressed_rest_Fix1Base_run*1000;
                onset_Ind_M2S           = Data.trial(TrialIndx).Time_matchtosample_appeared_run*1000   -Data.trial(TrialIndx).Time_pressed_rest_Fix1Base_run*1000;
                onset_Ind_M2S_selected  = Data.trial(TrialIndx).Time_match_to_sample_selected_run*1000   -Data.trial(TrialIndx).Time_pressed_rest_Fix1Base_run*1000;
                onset_Ind_PostCertainty = Data.trial(TrialIndx).Wager_start_time_wagering_post_run*1000  -Data.trial(TrialIndx).Time_pressed_rest_Fix1Base_run*1000;
                onset_Ind_PreCertainty  = Data.trial(TrialIndx).Wager_start_time_wagering_pre_run*1000   -Data.trial(TrialIndx).Time_pressed_rest_Fix1Base_run*1000;
                yPosition = -2.5;
                line([onset_IndSample onset_IndSample], get(gca,'YLim'),'Color','black','LineStyle','--');
                txt1 = 'Sample'; text(typecast(double(onset_IndSample), 'double'),min(L),txt1)
                line([onset_Ind_M2S onset_Ind_M2S], get(gca,'YLim'),'Color','black','LineStyle','--');
                txt1 = 'M2S';text(typecast(double(onset_Ind_M2S), 'double'),min(L),txt1)
                line([onset_Ind_M2S_selected onset_Ind_M2S_selected], get(gca,'YLim'),'Color','black','LineStyle','--');
                line([onset_Ind_PostCertainty onset_Ind_PostCertainty], get(gca,'YLim'),'Color','black','LineStyle','--');
                txt1 = 'Post Cer';text(typecast(double(onset_Ind_PostCertainty), 'double'),min(L),txt1)
                line([onset_Ind_PreCertainty onset_Ind_PreCertainty], get(gca,'YLim'),'Color','black','LineStyle','--');
                txt1 = 'Pre Cer';text(typecast(double(onset_Ind_PreCertainty), 'double'),min(L),txt1)
                line([onset_Ind_TrialStart onset_Ind_TrialStart], get(gca,'YLim'),'Color','black','LineStyle','--');
                txt1 = 'Trial Start';text(typecast(double(onset_Ind_TrialStart), 'double'),min(L),txt1)
                title(['Pupil Diameter   ', files(fileIndx).name(1:6), '   NrTrial:',num2str(TrialIndx) ,'   (',num2str(TrialIndex),')'])
                %
                print(p0,['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results\Pupil\Graph\Graph_rawData\',...
                    files(fileIndx).name(1:end-4),'_', num2str(TrialIndx),'_', num2str(TrialIndex)],'-dpng','-r0') %dpng
                close all;
                
                p1 = figure(2);
                plot((t_ms), vel(t_ms),'.');
                ylabel('difference to the previous sample (ms)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
                xlabel('Time (s)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
                print(p1,['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results\Pupil\Graph\Graph_SamplingRate\',...
                    files(fileIndx).name(1:end-4),'_', num2str(TrialIndx),'_VelocityOfSamplingTime'],'-dpng','-r0') %dpng
            end
            
            %% extract trial information
            % field_names = fieldnames(Data.trial); %All VariableNames
            %
            % % Process events (only extract relevant events):
            % evtRows = arrayfun(@(f) ~isempty(f.message),rawEDF.FEVENT);
            % eventData.t    = double((vertcat(rawEDF.FEVENT(evtRows).sttime)...
            %     -zeroTime_ms))/1000;
            % eventData.name = {rawEDF.FEVENT(evtRows).message}';
            %
            % % Find trials start and end times:
            % trialStartRows = find(strcmp(Data.trial,'save_Eye_y'));
            %
            %
            % Data.trial(j).find(strcmp(field_names, 'Trial_start_time_run')))
            % assert(length(trialStartRows)==36);
            % segmentStart   = eventData.t(trialStartRows);
            % segmentEnd     = eventData.t(trialStartRows+9);
            %
            % % Extract trial metadata:
            % trialNum = regexp(eventData.name...
            %     ,'!V TRIAL_VAR trial (\d{2,3})','tokens','once');
            % assert(isequal(find(~cellfun(@isempty,trialNum)),trailStartRows+5))
            % trial = str2double([trialNum{:}]');
            %
            % % Extract CA metadata:
            % trialCA = regexp(eventData.name...
            %     ,'!V TRIAL_VAR CA (\d|?)','tokens','once');
            % assert(isequal(find(~cellfun(@isempty,trialCA)),trailStartRows+7))
            % CA = str2double([trialCA{:}]');
            %
            % % Extract Condition metadata:
            % trialCondition = regexp(eventData.name...
            %     ,'!V TRIAL_VAR Condition (\d)','tokens','once');
            % assert(isequal(find(~cellfun(@isempty,trialCondition)),trailStartRows+8))
            % Condition = str2double([trialCondition{:}]');
            %
            % % Extract picture metadata:
            % trialPicture = regexp(eventData.name...
            %     ,'!V TRIAL_VAR picture (.)+$','tokens','once');
            % assert(isequal(find(~cellfun(@isempty,trialPicture)),trailStartRows+6))
            % picture = [trialPicture{:}]';
            %
            % % Make table:
            % segmentName  = strcat('pictureSeq_'...
            %     ,strrep(cellstr(num2str((1:36)')),' ','0'));
            % SegmentSource      = segmentName;
            % [~,justFileName,~] = fileparts(edfFilename);
            % SegmentSource(:)   = {justFileName};
            % fileType           = SegmentSource;
            % fileType(:)        = regexp(edfFilename...
            %     ,'.+(F|B|iCB1|iCB2).edf','tokens','once');
            % eyeUsed           = segmentName;
            % eyeUsed(:)        = {curEye};
            % segmentData = table(...
            %     segmentStart,segmentEnd,segmentName,SegmentSource...
            %     ,fileType,trial,CA,Condition,picture,eyeUsed);
            %
            %
            
            
            segmentData = table();
            %...
            %     segmentStart,segmentEnd,segmentName,SegmentSource...
            %     ,fileType,trial,CA,Condition,picture,eyeUsed);
            
            cd(['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results\Pupil\matfiles\',files(fileIndx).name(1:6)])
            
            % Build RawFileModel instance, which saves the data to a mat file that is
            % compatible with the other data models:
            RawFileModel(diameterUnit,diameter,segmentData...
                ,zeroTime_ms,matFilename);
            
            
            
        end
    end %participants
end

%% Preprocess the pupil-data for each trial

path       = 'C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results\Pupil\matfiles\';
Folder     = dir(path );
Folder(1:2)=[];
for i_Folder =3 1:length(Folder)
    folderName = [path ,Folder(i_Folder).name, '\'];
    rawFiles         = dir([folderName '*.mat']);
    Names = sort_nat({rawFiles.name});
    rawFiles = cell2struct(Names, 'name', 1);
    %for fileIndx = 1:length(rawFiles)
    
    close all;
    % Get settings the standard setttings:
    customSettings  = PupilDataModel.getDefaultSettings();
    
    % Customize settings (the current dataset features arbitrary units; as
    % such, no absolute maximum can be applied. Instead, the max is set to
    % infinite, and the min to a really low values):
    customSettings.raw.PupilDiameter_Max = inf;
    customSettings.raw.PupilDiameter_Min = 0.1;
    
    % The following flag determines if the raw data filter saves information
    % about each intermediate filter step. Enabling it uses considerably more
    % memory, and makes plotting multiple files slow. Use it only when
    % designing and tweaking the filter parameters.
    customSettings.raw.keepFilterData    = true;
    % Isolated sample filter criteria:
    
    % 'Sample-islands' are clusters of samples that are temporally seperated
    % from other samples. The minimum distance used to consider
    % samples 'separated' is specified below:
    customSettings.islandFilter_islandSeperation_ms    = 400;   %40[ms]
    
    % When clusters are seperated, they need to have the minimum size
    % specified below. Sample islands smaller than this temporal width and
    % separated from other samples by the distance mentioned above are
    % marked as invalid.
    customSettings.islandFilter_minIslandWidth_ms      = 500;   %50[ms]
    
    %remove gaps - gaps_ms > obj.settings.interp_maxGap
    customSettings.valid.interp_maxGap = 20000; %250ms
    % Construct one PupilDataModel instance per file using the batchConstructor
    % method:
    hPupilData = PupilDataModel.batchConstructor(...
        folderName,rawFiles,customSettings);
    
    %% plot the Sampling rate
    
    if preprocessing.GraphsPreprocessing
        
        p0 = figure(1);
        plot((hPupilData.timestamps_RawData_ms/1000), vel(hPupilData.timestamps_RawData_ms),'.');
        ylabel('difference to the previous sample (ms)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
        xlabel('Time (s)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
        %print(p0,['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Analyze\Pupil\Graph_SamplingRate\',...
        %rawFiles(fileIndx).name(1:end-4),'_VelocityOfSamplingTime'],'-dpng','-r0') %dpng
    end
    %% Filter
    %[valOut,speedFiltData,devFiltData]  =rawDataFilter(hPupilData.timestamps_RawData_ms, hPupilData.leftPupil_RawData.rawSample,customSettings.raw);
    hPupilData.filterRawData();
    if preprocessing.GraphsPreprocessing
        path_beh = 'C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Data\behavior\';
        files = dir([path_beh ,'*.mat']);
        Data = load([path_beh, files(fileIndx).name]);
        succesfull_Trials = find([Data.trial.unsuccess]== 0);
        
        for NrTrial =  1: length(hPupilData)
            p1 =figure(2);
            diameter.L_filtered = hPupilData(NrTrial).leftPupil_RawData.rawSample;
            diameter.L_filtered(~hPupilData(NrTrial).leftPupil_RawData.isValid  ) = nan;
            plot((hPupilData(NrTrial).timestamps_RawData_ms/1000), hPupilData(NrTrial).leftPupil_RawData.rawSample, 'b.'); hold on;
            plot((hPupilData(NrTrial).timestamps_RawData_ms/1000), diameter.L_filtered, 'r.'); hold on;
            ylabel('Pupil Dilation ','fontsize',WritingLabelAxis_Size,'fontweight','b' );
            xlabel('Time (s)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
            title('Blue - raw Data; Red -> removed Samples to Filter procedure','fontsize',14,'fontweight','b' );
            onset_Ind_TrialStart = Data.trial(succesfull_Trials(NrTrial)).Time_pressed_rest_Fix1Base_run - Data.trial(succesfull_Trials(NrTrial)).Time_pressed_rest_Fix1Base_run;
            onset_IndSample         = Data.trial(succesfull_Trials(NrTrial)).Time_sample_appeared_run -Data.trial(succesfull_Trials(NrTrial)).Time_pressed_rest_Fix1Base_run;
            yPosition = -2.5;
            line([onset_IndSample onset_IndSample], get(gca,'YLim'),'Color','black','LineStyle','--');
            txt1 = 'Sample'; text(typecast(double(onset_IndSample), 'double'),min(hPupilData(NrTrial).leftPupil_RawData.rawSample),txt1)
            
            line([onset_Ind_TrialStart onset_Ind_TrialStart], get(gca,'YLim'),'Color','black','LineStyle','--');
            txt1 = 'Trial Start';text(typecast(double(onset_Ind_TrialStart), 'double'),min(hPupilData(NrTrial).leftPupil_RawData.rawSample),txt1)
            
            print(p1,['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results\Pupil\Graph\Graph_Filtered\',...
                rawFiles(fileIndx).name(1:end-6),'_', num2str(NrTrial), '_Filtered'],'-dpng','-r0') %dpng
            close all;
        end
    end
    %% Interpolating and Filtering the Valid Samples:
    % Once the valid samples are identified, they can be used to create a
    % smooth high-resolution pupil size signal through interpolation and
    % low-pass filtering, which is what the processValidSamples method does.
    % Formula..    diaInterp = interp1(valid_t_ms./1000 ,validDiams,t_upsampled,'linear');
    
    % Process the valid samples:
    hPupilData.processValidSamples();
    
    if preprocessing.GraphsPreprocessing
        
        %Graph to see the smoothing and interpolation
        for NrTrial =  1: length(hPupilData)
            p2 =figure(4);
            plot((hPupilData(NrTrial).timestamps_RawData_ms/1000), hPupilData(NrTrial).leftPupil_RawData.rawSample, 'b.'); hold on;
            %plot((hPupilData(NrTrial).leftPupil_ValidSamples.signal.t),hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap ,'g.')
            plot((hPupilData(NrTrial).leftPupil_ValidSamples.signal.t),hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter ,'g.')
            
            onset_Ind_TrialStart = Data.trial(succesfull_Trials(NrTrial)).Time_pressed_rest_Fix1Base_run - Data.trial(succesfull_Trials(NrTrial)).Time_pressed_rest_Fix1Base_run;
            onset_IndSample         = Data.trial(succesfull_Trials(NrTrial)).Time_sample_appeared_run -Data.trial(succesfull_Trials(NrTrial)).Time_pressed_rest_Fix1Base_run;
            yPosition = -2.5;
            line([onset_IndSample onset_IndSample], get(gca,'YLim'),'Color','black','LineStyle','--');
            txt1 = 'Sample'; text(typecast(double(onset_IndSample), 'double'),min(hPupilData(NrTrial).leftPupil_RawData.rawSample),txt1)
            
            line([onset_Ind_TrialStart onset_Ind_TrialStart], get(gca,'YLim'),'Color','black','LineStyle','--');
            txt1 = 'Trial Start';text(typecast(double(onset_Ind_TrialStart), 'double'),min(hPupilData(NrTrial).leftPupil_RawData.rawSample),txt1)
            
            
            ylabel('Pupil Dilation ','fontsize',WritingLabelAxis_Size,'fontweight','b' );
            xlabel('Time (ms)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
            title('Blue: raw Data; Green: Interpolated','fontsize',14,'fontweight','b' );
            print(p2,['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results\Pupil\Graph\Graph_SmoothedInterpolated\',...
                rawFiles(fileIndx).name(1:end-6),'_', num2str(NrTrial), '_smoothed_Interpolated'],'-dpng','-r0') %dpng
            close all;
        end
    end
     
    %%
    files = dir([path_beh ,'*.mat']);
    Data = load([path_beh, files(i_Folder).name]);
    succesfull_Trials = find([Data.trial.unsuccess]== 0);
    
    %% plot average timeCourse for the 5 difficulty levels
    % Which trials belong to which difficulty level?
    NrTrial = 0;
    DaTab = [];
    for indtrial = succesfull_Trials;
        NrTrial =    NrTrial +1;
        DaTab(NrTrial).NrTrial              = NrTrial; % how many completed trials..
        DaTab(NrTrial).Ind_CompletedTrials  = indtrial; %index of the trial in the group of all started trials
        DaTab(NrTrial).difficultyLevel      = abs([Data.trial(indtrial).refAngle_M2S] - [Data.trial(indtrial).refAngle_Sample] );
        DaTab(NrTrial).Correct              = [Data.trial(indtrial).unsuccess] == 0 ;
        DaTab(NrTrial).Incorrect            = [Data.trial(indtrial).unsuccess] == 1 ;
        DaTab(NrTrial).PostCertaintyLevel   = [Data.trial(indtrial).wager_choosen_post]; 
        DaTab(NrTrial).PreCertaintyLevel    = [Data.trial(indtrial).wager_choosen_pre];
        DaTab(NrTrial).RatingsOrControl  = [Data.trial(indtrial).wagering_or_controll_wagering_post]; %1- rated preCertainty, 2%rated postCertainty, 3- both ratings
    end
    %Pupil Diameter.. Average all trials belonging to one difficulty level for each time point
    % Which time points are identical for each trial to be averaged?
    % did we change the number of samples per trial? yes
    
    for NrTrial =  1: length(DaTab);
        % time for events is continous increasing from trial to trial
        DaTab(NrTrial).Time_pressed_rest_Fix1Base_run   = Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run;
        DaTab(NrTrial).Trial_finish_time_run            = Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Trial_finish_time_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run;
        DaTab(NrTrial).Trial_finish_time_run_pupil      = hPupilData(NrTrial).leftPupil_ValidSamples.signal.t(end);
        DaTab(NrTrial).Time_sample_appeared_run         = Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_sample_appeared_run;
        DaTab(NrTrial).Time_sample_appeared_run_trial   = Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_sample_appeared_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run;
        DaTab(NrTrial).Wager_start_time_wagering_post_run   = Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Wager_start_time_wagering_post_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run;
        DaTab(NrTrial).Wager_end_time_wagering_post_run     = Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Wager_end_time_wagering_post_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run;
        DaTab(NrTrial).Time_pressed_rest_Fix3_MTS_run       = Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix3_MTS_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run;

        index_SampleAppears = [];  index_SampleDisappears = []; index_PercChoiceAppears = []; %- Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Trial_start_time_run
        [c ,DaTab(NrTrial).index_TrialStart]          = min(abs(hPupilData(NrTrial).leftPupil_ValidSamples.signal.t - (Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run)));
        [c ,DaTab(NrTrial).index_SampleAppears]       = min(abs(hPupilData(NrTrial).leftPupil_ValidSamples.signal.t - (Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_sample_appeared_run       -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run)));
        [c ,DaTab(NrTrial).index_SampleDisappears]    = min(abs(hPupilData(NrTrial).leftPupil_ValidSamples.signal.t - (Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_sample_disappeared_run    -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run )));
        [c ,DaTab(NrTrial).index_PercChoiceAppears]   = min(abs(hPupilData(NrTrial).leftPupil_ValidSamples.signal.t - (Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_matchtosample_appeared_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run)));
        [c ,DaTab(NrTrial).index_PercChoice_Selected] = min(abs(hPupilData(NrTrial).leftPupil_ValidSamples.signal.t - (Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_match_to_sample_selected_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run)));
        [c ,DaTab(NrTrial).index_PercChoice_Selected_TimeLater] = min(abs(hPupilData(NrTrial).leftPupil_ValidSamples.signal.t - ((Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_match_to_sample_selected_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run) +4)));
        
        [c ,DaTab(NrTrial).index_StartRest_afterM2S] = min(abs(hPupilData(NrTrial).leftPupil_ValidSamples.signal.t - ((Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix3_MTS_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run) )));
        [c ,DaTab(NrTrial).index_pressed_rest_Fix2_Pre] = min(abs(hPupilData(NrTrial).leftPupil_ValidSamples.signal.t - ((Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix2_Pre_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run) )));

        [c ,DaTab(NrTrial).index_PostCertaintyAppeared]    = min(abs(hPupilData(NrTrial).leftPupil_ValidSamples.signal.t - ((Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Wager_start_time_wagering_post_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run))));
        [c ,DaTab(NrTrial).index_PostCertaintyDisappeared] = min(abs(hPupilData(NrTrial).leftPupil_ValidSamples.signal.t - ((Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Wager_end_time_wagering_post_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run))));
        [c ,DaTab(NrTrial).index_PostCertainty_1streleasedButton] = min(abs(hPupilData(NrTrial).leftPupil_ValidSamples.signal.t - ((Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Wager_time_1streleasedButton_post_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run))));

        [c ,DaTab(NrTrial).index_Rest_PostCert]       = min(abs(hPupilData(NrTrial).leftPupil_ValidSamples.signal.t - ((Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix4_Post_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run) )));
        [c ,DaTab(NrTrial).index_TrialEnd]            = min(abs(hPupilData(NrTrial).leftPupil_ValidSamples.signal.t - (Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Trial_finish_time_run - Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run)));
        
        DaTab(NrTrial).SamplesBetweenSample_PercChoice_Selected_TimeLater  = DaTab(NrTrial).index_PercChoice_Selected_TimeLater - DaTab(NrTrial).index_SampleAppears ;
        DaTab(NrTrial).SamplesBetweenRestM2S_EndTrial = DaTab(NrTrial).index_TrialEnd - DaTab(NrTrial).index_StartRest_afterM2S ;
        DaTab(NrTrial).Start_EndTrial = DaTab(NrTrial).index_TrialEnd - DaTab(NrTrial).index_TrialStart ;
        DaTab(NrTrial).SamplesBetweenSample_PercChoiceSelected_TimeLater = DaTab(NrTrial).index_PercChoice_Selected_TimeLater - DaTab(NrTrial).index_SampleAppears ;

    end
    
    %% plot  the pupil data with lines for important events
  count = 1;  Idx_NoSignal = NaN; 
    for NrTrial =  1: length(hPupilData) %Nr. of completed trials
        if length(hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter) < 2000 || nanstd(hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter)< 20
        DaTab(NrTrial).difficultyLevel = NaN; 
        DaTab(NrTrial).PostCertaintyLevel = NaN; 
        DaTab(NrTrial).SamplesBetweenSample_PercChoice_Selected_TimeLater = [];
        DaTab(NrTrial).SamplesBetweenRestM2S_EndTrial = [];
        Idx_NoSignal(count) = NrTrial; 
        printToConsole(2, ['  no pupil signal: delete the following trials: ', num2str(NrTrial)]);
        count = count+1; 
        end
    end
    %%
    if preprocessing.GraphsPreprocessing
        
        %Graph to see the smoothing and interpolation
        for NrTrial =  1: length(hPupilData) %Nr. of completed trials
            % plot the
            p4 =figure(6);
            plot((hPupilData(NrTrial).leftPupil_ValidSamples.signal.t),hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter ,'r.')
            
            title(['Pupil Diameter   ', files(fileIndx).name(1:6), '   NrTrial:',num2str(NrTrial)])
            ylabel('Pupil Dilation ','fontsize',WritingLabelAxis_Size,'fontweight','b' );
            xlabel('Time (s)','fontsize',WritingLabelAxis_Size,'fontweight','b' );
            onset_Ind_Time_TrialStart = Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Trial_start_time_run - Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Trial_start_time_run;
            onset_Ind_Time_Fix1Base = Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run - Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run;
            onset_IndSample         = Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_sample_appeared_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run;
            onset_Ind_M2S           = Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_matchtosample_appeared_run-Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run;
            onset_Ind_M2S_selected  = Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_match_to_sample_selected_run-Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run;
            onset_Ind_PostCertainty = Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Wager_start_time_wagering_post_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run;
            onset_Ind_PreCertainty  = Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Wager_start_time_wagering_pre_run -Data.trial(DaTab(NrTrial).Ind_CompletedTrials).Time_pressed_rest_Fix1Base_run;
            yPosition = -2.5;
            line([onset_IndSample onset_IndSample], get(gca,'YLim'),'Color','black','LineStyle','--');
            txt1 = 'Sample'; text(typecast(double(onset_IndSample), 'double'),min(hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap),txt1)
            line([onset_Ind_M2S onset_Ind_M2S], get(gca,'YLim'),'Color','black','LineStyle','--');
            line([onset_Ind_M2S_selected onset_Ind_M2S_selected], get(gca,'YLim'),'Color','black','LineStyle','--');
            line([onset_Ind_PostCertainty onset_Ind_PostCertainty], get(gca,'YLim'),'Color','black','LineStyle','--');
            txt1 = 'Post Cer';text(typecast(double(onset_Ind_PostCertainty), 'double'),min(hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap),txt1)
            line([onset_Ind_PreCertainty onset_Ind_PreCertainty], get(gca,'YLim'),'Color','black','LineStyle','--');
            txt1 = 'Pre Cer';text(typecast(double(onset_Ind_PreCertainty), 'double'),min(hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap),txt1)
            line([onset_Ind_Time_TrialStart onset_Ind_Time_TrialStart], get(gca,'YLim'),'Color','black','LineStyle','--');
            txt1 = 'Trial Start';text(typecast(double(onset_Ind_Time_TrialStart), 'double'),min(hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap),txt1)
            line([onset_Ind_Time_Fix1Base onset_Ind_Time_Fix1Base], get(gca,'YLim'),'Color','black','LineStyle','--');
            txt1 = 'Fix1';text(typecast(double(onset_Ind_Time_Fix1Base), 'double'),min(hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap)+20,txt1)
            
            
            print(p4,['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results\Pupil\Graph\Graph_SmoothedInterpolated_MarkedEvents\',...
                rawFiles(fileIndx).name(1:end-6),'_', num2str(NrTrial), '_smoothed_Interpolated'],'-dpng','-r0') %dpng
            close all; 
        end
    end
    
    
    %% BaselineCorrection - Difficulty Levels
   BaselineCorrection = 'SingleTimePoint';
    BaselineCorrection = 'Average';
    for NrTrial = 1:length(hPupilData)
       p2 = [];  p1 = [];  t2 = []; 

        if ~ismember(Idx_NoSignal, NrTrial)
        switch BaselineCorrection
            case 'SingleTimePoint'
                hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap_Baselinecorected_Diff   = hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap   -  hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap(DaTab(NrTrial).index_PercChoiceAppears);
                t2 = hPupilData(NrTrial).leftPupil_ValidSamples.signal.t(DaTab(NrTrial).index_PercChoiceAppears);
                hPupilData(NrTrial).leftPupil_ValidSamples.signal.t_Baselinecorected_Diff   = hPupilData(NrTrial).leftPupil_ValidSamples.signal.t  ...
                    -  t2;
            case 'Average'
                p2 = hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap(DaTab(NrTrial).index_PercChoiceAppears -200);
                p1 = hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap(DaTab(NrTrial).index_PercChoiceAppears);
                t2 = hPupilData(NrTrial).leftPupil_ValidSamples.signal.t(DaTab(NrTrial).index_PercChoiceAppears);

                if p2 > p1
                   Av =  mean(p1:p2); 
                elseif p2 < p1 
                    Av =  mean(p2:p1); 
                else
                    Av =  mean(p1:p2); 
                    printToConsole(2, ['  p1 and p2 equal size ', num2str(NrTrial)]);
                end
                
                hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap_Baselinecorected_Diff   = (hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap  ...
                    -  Av )/Av;
                hPupilData(NrTrial).leftPupil_ValidSamples.signal.t_Baselinecorected_Diff   = hPupilData(NrTrial).leftPupil_ValidSamples.signal.t  ...
                    -  t2;      
        end
        end
    end
    
    % Prepare the plots for each Difficulty Level
    difficultyLevels                            = unique([DaTab.difficultyLevel]);
    difficultyLevels(isnan(difficultyLevels))   =[]; 
    %initiave the vector
    averagePupildiamter                 = NaN( min([DaTab.SamplesBetweenSample_PercChoice_Selected_TimeLater]),length(difficultyLevels));
    averagePupildiamter_Corrected       = NaN( min([DaTab.SamplesBetweenSample_PercChoice_Selected_TimeLater]),length(difficultyLevels));
    Time_Diff                           = NaN( min([DaTab.SamplesBetweenSample_PercChoice_Selected_TimeLater]),length(difficultyLevels));
    %average the pupildiameter for each sample per difficulty level
    DiffLevel = []; 

    for i = 1:length(difficultyLevels)
        counter = 0;
%         if min([DaTab.SamplesBetweenSample_PercChoice_Selected_TimeLater]) == 0
%            idx =  find([DaTab(:).SamplesBetweenSample_PercChoice_Selected_TimeLater] == 0)
%         end
        for i_Sample = 1: min([DaTab.SamplesBetweenSample_PercChoice_Selected_TimeLater])
            DiffLevel(i,:).ind                                  = find([DaTab.difficultyLevel] == difficultyLevels(i));
            DiffLevel(i,:).index_SampleAppears                  = arrayfun(@(x)x.index_SampleAppears, DaTab(DiffLevel(i,:).ind));
            DiffLevel(i,:).index_PercChoice_Selected_TimeLater  = arrayfun(@(x)x.index_PercChoice_Selected_TimeLater, DaTab(DiffLevel(i,:).ind ));
            DiffLevel(i,:).indCounter                           = DiffLevel(i,:).index_SampleAppears +counter; %
            counter = counter +1;
            % Sample Appears until Selection +4s
            averagePupildiamter(i_Sample,i)             = nanmean(arrayfun(@(x,y) x.leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap(y), hPupilData(DiffLevel(i,:).ind),DiffLevel(i,:).indCounter'));
            averagePupildiamter_Corrected(i_Sample,i)   = nanmean(arrayfun(@(x,y) x.leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap_Baselinecorected_Diff(y), hPupilData(DiffLevel(i,:).ind),DiffLevel(i,:).indCounter'));
            Time_Diff(i_Sample,i)                       = nanmean(arrayfun(@(x,y) x.leftPupil_ValidSamples.signal.t_Baselinecorected_Diff(y), hPupilData(DiffLevel(i,:).ind),DiffLevel(i,:).indCounter'));
        end
    end
    
    %% subplots for each difficulty - time period specified
    % baselineCorrection
    p5 =figure(5);
    for i_diff = 1:length(difficultyLevels)
        subplot(length(difficultyLevels),1,i_diff);
        for NrTrial =  DiffLevel(i_diff,:).ind
            xTime = hPupilData(NrTrial).leftPupil_ValidSamples.signal.t_Baselinecorected_Diff(DaTab(NrTrial).index_SampleAppears:DaTab(NrTrial).index_PercChoice_Selected_TimeLater)-hPupilData(NrTrial).leftPupil_ValidSamples.signal.t(DaTab(NrTrial).index_SampleAppears);
            yPuDiameter = hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap_Baselinecorected_Diff(DaTab(NrTrial).index_SampleAppears:DaTab(NrTrial).index_PercChoice_Selected_TimeLater)-hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap(DaTab(NrTrial).index_PercChoiceAppears);
            plot(xTime, yPuDiameter);hold on;
        end        
    end
   ylabel('Relative pupil diameter','fontsize',WritingLabelAxis_Size,'fontweight','b' );

    print(p5,['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results\Pupil\Graph\Graph_DifficultyLevels\',...
        files(i_Folder).name(1:end-4), 'DifficultyLevels_EachTrial2'],'-dpng','-r0') %dpng
    %% Plot average difficulty
    Plotcolor_D = [{'.' };{'.'};{'.' }; {'.' }; {'.' }]; %difficult to easy
    Plotcolor_D2 = [[1,0,0] ;[0.7,0.4,0 ] ;[0.9,0.7,0.1];[0.5,0.7,0.2] ;[0, 0.5, 0 ]]; %difficult to easy

    p6 =figure(6);
    for i_diff = 1:length(difficultyLevels)
        plot(  Time_Diff(:,i_diff) , averagePupildiamter(:,i_diff) , Plotcolor_D{i_diff},'MarkerEdge',Plotcolor_D2(i_diff,:) , 'MarkerFace',Plotcolor_D2(i_diff,:) , 'MarkerSize',6,'DisplayName',num2str(difficultyLevels(i_diff)))   ;
        Diff{i_diff} = num2str(difficultyLevels(i_diff));
        legend(Diff, 'FontSize',12);hold on  
    end
       ylabel('Relative pupil diameter','fontsize',WritingLabelAxis_Size,'fontweight','b' );
       xlabel('Time from stimulus onset','fontsize',WritingLabelAxis_Size,'fontweight','b' );

    line([0 0], get(gca,'YLim'),'Color','black','LineStyle','--');
    txt1 = '2nd Stimuli Appears '; text(typecast(double(0), 'double'),min(min(averagePupildiamter_Corrected)),txt1)
    line([-1 -1], get(gca,'YLim'),'Color','black','LineStyle','--');
    txt1 = 'mask';text(typecast(double(-1), 'double'),min(min(averagePupildiamter_Corrected)),txt1)
    line([-2 -2], get(gca,'YLim'),'Color','black','LineStyle','--');
    txt1 = 'Sample';text(typecast(double(-2), 'double'),min(min(averagePupildiamter_Corrected)),txt1)
    % perceptual reaction time
    print(p6,['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results\Pupil\Graph\Graph_DifficultyLevels\',...
        files(i_Folder).name(1:end-4), 'DifficultyLevels_AverageBaselineCorrectBefore2ndStimuli2'],'-dpng','-r0') %dpng
    
    
    
    
    
    p7 =figure(7);
    for i_diff = 1:length(difficultyLevels)
        plot(  Time_Diff(:,i_diff) , averagePupildiamter_Corrected(:,i_diff) , Plotcolor_D{i_diff},'MarkerEdge',Plotcolor_D2(i_diff,:) , 'MarkerFace',Plotcolor_D2(i_diff,:) , 'MarkerSize',6,'DisplayName',num2str(difficultyLevels(i_diff)))   ;
        Diff{i_diff} = num2str(difficultyLevels(i_diff));
        legend(Diff, 'FontSize',12);hold on  
    end
 ylabel('Relative pupil diameter','fontsize',WritingLabelAxis_Size,'fontweight','b' );
       xlabel('Time from stimulus onset','fontsize',WritingLabelAxis_Size,'fontweight','b' );

    line([0 0], get(gca,'YLim'),'Color','black','LineStyle','--');
    txt1 = '2nd Stimuli Appears '; text(typecast(double(0), 'double'),min(min(averagePupildiamter_Corrected)),txt1)
    line([-1 -1], get(gca,'YLim'),'Color','black','LineStyle','--');
    txt1 = 'mask';text(typecast(double(1000), 'double'),min(min(averagePupildiamter_Corrected)),txt1)
    line([-2 -2], get(gca,'YLim'),'Color','black','LineStyle','--');
    txt1 = 'Sample';text(typecast(double(2000), 'double'),min(min(averagePupildiamter_Corrected)),txt1)
    % perceptual reaction time
    print(p7,['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results\Pupil\Graph\Graph_DifficultyLevels\',...
        files(i_Folder).name(1:end-4), 'DifficultyLevels_NoBaseline2'],'-dpng','-r0') %dpng
  
    
    close all;
   % end %all trials per subject
%% Difficulty Level - show the complete TRIAL
  % Prepare the plots for each Difficulty Level
    difficultyLevels                            = unique([DaTab.difficultyLevel]);
    difficultyLevels(isnan(difficultyLevels))   =[]; 
    %initiave the vector
    averagePupildiamter                 = NaN( min([DaTab.Start_EndTrial]),length(difficultyLevels));
    averagePupildiamter_Corrected       = NaN( min([DaTab.Start_EndTrial]),length(difficultyLevels));
    Time_Diff                           = NaN( min([DaTab.Start_EndTrial]),length(difficultyLevels));

    %average the pupildiameter for each sample per difficulty level
    DiffLevel_complTrial = []; 

    for i = 1:length(difficultyLevels)
        counter = 0;
%         if min([DaTab.SamplesBetweenSample_PercChoice_Selected_TimeLater]) == 0
%            idx =  find([DaTab(:).SamplesBetweenSample_PercChoice_Selected_TimeLater] == 0)
%         end
        for i_Sample = 1: min([DaTab.Start_EndTrial])
            DiffLevel_complTrial(i,:).ind                                  = find([DaTab.difficultyLevel] == difficultyLevels(i));
            DiffLevel_complTrial(i,:).index_TrialStart                 = arrayfun(@(x)x.index_TrialStart, DaTab(DiffLevel_complTrial(i,:).ind));
            DiffLevel_complTrial(i,:).index_TrialEnd  = arrayfun(@(x)x.index_TrialEnd, DaTab(DiffLevel_complTrial(i,:).ind ));
            DiffLevel_complTrial(i,:).indCounter                           = DiffLevel_complTrial(i,:).index_TrialStart +counter; %
            counter = counter +1;
            % Sample Appears until Selection +4s
            averagePupildiamter(i_Sample,i)             = nanmean(arrayfun(@(x,y) x.leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap(y), hPupilData(DiffLevel_complTrial(i,:).ind),DiffLevel_complTrial(i,:).indCounter'));
            averagePupildiamter_Corrected(i_Sample,i)   = nanmean(arrayfun(@(x,y) x.leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap_Baselinecorected_Diff(y), hPupilData(DiffLevel_complTrial(i,:).ind),DiffLevel_complTrial(i,:).indCounter'));
            Time_Diff(i_Sample,i)                       = nanmean(arrayfun(@(x,y) x.leftPupil_ValidSamples.signal.t_Baselinecorected_Diff(y), hPupilData(DiffLevel_complTrial(i,:).ind),DiffLevel_complTrial(i,:).indCounter'));
        end
    end
    
   %% Plot average difficulty for the complete trial
    Plotcolor_D = [{'.' };{'.'};{'.' }; {'.' }; {'.' }]; %difficult to easy
    Plotcolor_D2 = [[1,0,0] ;[0.7,0.4,0 ] ;[0.9,0.7,0.1];[0.5,0.7,0.2] ;[0, 0.5, 0 ]]; %difficult to easy

    p6 =figure(6);
    for i_diff = 1:length(difficultyLevels)
        plot(  Time_Diff(:,i_diff) , averagePupildiamter(:,i_diff) , Plotcolor_D{i_diff},'MarkerEdge',Plotcolor_D2(i_diff,:) , 'MarkerFace',Plotcolor_D2(i_diff,:) , 'MarkerSize',6,'DisplayName',num2str(difficultyLevels(i_diff)))   ;
        Diff{i_diff} = num2str(difficultyLevels(i_diff));
        legend(Diff, 'FontSize',12);hold on  
    end
       ylabel('Relative pupil diameter','fontsize',WritingLabelAxis_Size,'fontweight','b' );
       xlabel('Time from stimulus onset','fontsize',WritingLabelAxis_Size,'fontweight','b' );

    line([0 0], get(gca,'YLim'),'Color','black','LineStyle','--');
    txt1 = '2nd Stimuli Appears '; text(typecast(double(0), 'double'),min(min(averagePupildiamter_Corrected)),txt1)
    line([-1 -1], get(gca,'YLim'),'Color','black','LineStyle','--');
    txt1 = 'mask';text(typecast(double(1000), 'double'),min(min(averagePupildiamter_Corrected)),txt1)
    line([-2 -2], get(gca,'YLim'),'Color','black','LineStyle','--');
    txt1 = 'Sample';text(typecast(double(2000), 'double'),min(min(averagePupildiamter_Corrected)),txt1)
    % perceptual reaction time
    print(p6,['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results\Pupil\Graph\Graph_DifficultyLevels\',...
        files(i_Folder).name(1:end-4), 'DifficultyLevels_NotbaselineCorrected2ndStimuli_CompleteTrial'],'-dpng','-r0') %dpng
    
    
    
    
    
    p7 =figure(7);
    for i_diff = 1:length(difficultyLevels)
        plot(  Time_Diff(:,i_diff) , averagePupildiamter_Corrected(:,i_diff) , Plotcolor_D{i_diff},'MarkerEdge',Plotcolor_D2(i_diff,:) , 'MarkerFace',Plotcolor_D2(i_diff,:) , 'MarkerSize',6,'DisplayName',num2str(difficultyLevels(i_diff)))   ;
        Diff{i_diff} = num2str(difficultyLevels(i_diff));
        legend(Diff, 'FontSize',12);hold on  
    end
 ylabel('Relative pupil diameter','fontsize',WritingLabelAxis_Size,'fontweight','b' );
       xlabel('Time from stimulus onset','fontsize',WritingLabelAxis_Size,'fontweight','b' );

    line([0 0], get(gca,'YLim'),'Color','black','LineStyle','--');
    txt1 = '2nd Stimuli Appears '; text(typecast(double(0), 'double'),min(min(averagePupildiamter_Corrected)),txt1)
    line([-1 -1], get(gca,'YLim'),'Color','black','LineStyle','--');
    txt1 = 'mask';text(typecast(double(1000), 'double'),min(min(averagePupildiamter_Corrected)),txt1)
    line([-2 -2], get(gca,'YLim'),'Color','black','LineStyle','--');
    txt1 = 'Sample';text(typecast(double(2000), 'double'),min(min(averagePupildiamter_Corrected)),txt1)
    % perceptual reaction time
    print(p7,['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results\Pupil\Graph\Graph_DifficultyLevels\',...
        files(i_Folder).name(1:end-4), 'DifficultyLevels_baselineCorrected2ndStimuli_CompleteTrial'],'-dpng','-r0') %dpng
  
    
    close all; 
   
   
   
%% BaselineCorrection - Certainty Levels
% Subtract the average value for a specific period of one trial from the
% PupilData of this trial
   BaselineCorrection = 'SingleTimePoint';
    BaselineCorrection = 'Average';
    for NrTrial = 1:length(hPupilData)
       if ~ismember(Idx_NoSignal, NrTrial)

        switch BaselineCorrection
            case 'SingleTimePoint'
                hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap_Baselinecorected_PostCert   = ...
                    hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap   -  ...
                    hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap(DaTab(NrTrial).index_PostCertaintyAppeared);
                
                t2 = hPupilData(NrTrial).leftPupil_ValidSamples.signal.t(DaTab(NrTrial).index_PostCertaintyAppeared);
                hPupilData(NrTrial).leftPupil_ValidSamples.signal.t_Baselinecorected_PostCert   = hPupilData(NrTrial).leftPupil_ValidSamples.signal.t  ...
                    -  t2; 
            case 'Average'
                p1 = hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap((DaTab(NrTrial).index_PostCertaintyAppeared -200));
                p2 = hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap(DaTab(NrTrial).index_PostCertaintyAppeared);
                t2 = hPupilData(NrTrial).leftPupil_ValidSamples.signal.t(DaTab(NrTrial).index_PostCertaintyAppeared);

              if p2 > p1
                   Av =  mean(p1:p2); 
                elseif p2 < p1 
                    Av =  mean(p2:p1); 
                else
                    Av =  mean(p1:p2); 
                    printToConsole(2, ['  p1 and p2 equal size ', num2str(NrTrial)]);
              end
                hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap_Baselinecorected_PostCert   = hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap  ...
                    -  Av;  
                hPupilData(NrTrial).leftPupil_ValidSamples.signal.t_Baselinecorected_PostCert   = hPupilData(NrTrial).leftPupil_ValidSamples.signal.t  ...
                    -  t2; 
        end
       end
    end
    
   % Prepare the plots for each CertaintyLevel
    PostCertaintyLevels = 1:6; %unique([DaTab.PostCertaintyLevel]);
    %initiave the vector
    averagePupildiamter                     = NaN( min([DaTab.SamplesBetweenRestM2S_EndTrial]),length(PostCertaintyLevels));
    averagePupildiamter_Corrected_PostCer   = NaN( min([DaTab.SamplesBetweenRestM2S_EndTrial]),length(PostCertaintyLevels));
    Time_PostCert                           = NaN( min([DaTab.SamplesBetweenRestM2S_EndTrial]),length(PostCertaintyLevels));
    %average the pupildiameter for each sample per difficulty level
    Certainty = []; 
    for i = 1:length(PostCertaintyLevels)
        counter = 0;
        % starting from Rest_M2S as first Sample, average all sample
        for i_Sample = 1: min([DaTab.SamplesBetweenRestM2S_EndTrial])
            Certainty(i,:).indPost                                  = find([DaTab.PostCertaintyLevel] == PostCertaintyLevels(i) & [DaTab.RatingsOrControl] == 2);
            Certainty(i,:).index_StartRest_afterM2S                           = arrayfun(@(x)x.index_StartRest_afterM2S, DaTab(Certainty(i,:).indPost));
            Certainty(i,:).index_TrialEnd                           = arrayfun(@(x)x.index_TrialEnd, DaTab(Certainty(i,:).indPost ));
            Certainty(i,:).indCounter                               = Certainty(i,:).index_StartRest_afterM2S + counter; %changes with each Sample
            counter = counter +1;
            % Sample Appears until Selection +4s
            averagePupildiamter(i_Sample,i)                     = nanmean(arrayfun(@(x,y) x.leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap(y),                            hPupilData(Certainty(i,:).indPost),Certainty(i,:).indCounter'));
            averagePupildiamter_Corrected_PostCer(i_Sample,i)   = nanmean(arrayfun(@(x,y) x.leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap_Baselinecorected_PostCert(y),  hPupilData(Certainty(i,:).indPost),Certainty(i,:).indCounter'));
            Time_PostCert(i_Sample,i)                           = nanmean(arrayfun(@(x,y) x.leftPupil_ValidSamples.signal.t_Baselinecorected_PostCert(y),                           hPupilData(Certainty(i,:).indPost),Certainty(i,:).indCounter'));
            Time(i_Sample,1)                   = nanmean(arrayfun(@(x,y) x.leftPupil_ValidSamples.signal.t(y)                                                            , hPupilData(Certainty(i,:).indPost),Certainty(i,:).indCounter'));
        end
    end
   %% subplots for each certainty level - showing each trial
    p5 =figure(5);
    for i_cert = 1:length(PostCertaintyLevels)
        subplot(length(PostCertaintyLevels),1,i_cert);
        for NrTrial =  Certainty(i_cert,:).indPost
            xTime = hPupilData(NrTrial).leftPupil_ValidSamples.signal.t_Baselinecorected_PostCert(DaTab(NrTrial).index_StartRest_afterM2S:DaTab(NrTrial).index_TrialEnd)        -       hPupilData(NrTrial).leftPupil_ValidSamples.signal.t(DaTab(NrTrial).index_StartRest_afterM2S);
            yPuDiameter = hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap_Baselinecorected_PostCert(DaTab(NrTrial).index_StartRest_afterM2S:DaTab(NrTrial).index_TrialEnd)       -       hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap(DaTab(NrTrial).index_PostCertaintyAppeared);
            plot(xTime, yPuDiameter);hold on;
            
        end        
    end
    print(p5,['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results\Pupil\Graph\Graph_PostCertainty\',...
        files(i_Folder).name(1:end-4), 'PostCertainty_EachTrial'],'-dpng','-r0') %dpng
    %% Plot average Certainty Levels
    Plotcolor_D = [{'.' };{'.'};{'.' }; {'.' }; {'.' }; {'.'}]; %difficult to easy
    Plotcolor_D2 = [[1,0,0.4] ;[0.7,0,1];[0.5,0.2,0.5] ;[0,0,1];[0,0.4,0.7] ;[0,1,1]]; %difficult to easy

    p6 =figure(6);
    for i_Postcert = 1:length(PostCertaintyLevels)
        plot(  Time_PostCert(:,i_Postcert) , averagePupildiamter(:,i_Postcert) , Plotcolor_D{i_Postcert},'MarkerEdge',Plotcolor_D2(i_Postcert,:) , 'MarkerFace',Plotcolor_D2(i_Postcert,:) , 'MarkerSize',6,'DisplayName',num2str(PostCertaintyLevels(i_Postcert)))   ;
       % Indx = find(round(Time_PostCert(:,i_Postcert),3) == 0); 
       % plot( Time_PostCert(Indx,i_Postcert) , averagePupildiamter_Corrected_PostCer(Indx,i_Postcert), 'ro')        
        %plot( 1:length(averagePupildiamter(:,i_Postcert)), averagePupildiamter_Corrected_PostCer(:,i_Postcert) , Plotcolor_D{i_Postcert}, 'MarkerSize',2,'DisplayName',num2str(PostCertaintyLevels(i_Postcert)))   ;
        PostCert{i_Postcert} = num2str(PostCertaintyLevels(i_Postcert));
        legend(PostCert, 'FontSize',12);hold on   
    end
    line([0 0], get(gca,'YLim'),'Color','black','LineStyle','--');
    txt1 = 'PostCert appears'; text(typecast(double(0), 'double'),min(min(averagePupildiamter)),txt1)

    print(p6,['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results\Pupil\Graph\Graph_PostCertainty\',...
        files(i_Folder).name(1:end-4), 'PostCertainty_NoBaseline'],'-dpng','-r0') %dpng
    
    
     p7 =figure(7);
    for i_Postcert = 1:length(PostCertaintyLevels)
        plot(  Time_PostCert(:,i_Postcert) , averagePupildiamter_Corrected_PostCer(:,i_Postcert) , Plotcolor_D{i_Postcert},'MarkerEdge',Plotcolor_D2(i_Postcert,:) , 'MarkerFace',Plotcolor_D2(i_Postcert,:) , 'MarkerSize',6,'DisplayName',num2str(PostCertaintyLevels(i_Postcert)))   ;
        PostCert{i_Postcert} = num2str(PostCertaintyLevels(i_Postcert));
        legend(PostCert, 'FontSize',12);hold on   
    end
    line([0 0], get(gca,'YLim'),'Color','black','LineStyle','--');
    txt1 = 'PostCert appears'; text(typecast(double(0), 'double'),min(min(averagePupildiamter_Corrected_PostCer)),txt1)

    print(p7,['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Perceptual_Certainty\Results\Pupil\Graph\Graph_PostCertainty\',...
        files(i_Folder).name(1:end-4), 'PostCertainty_AverageBaselineBeforePostCerRatingAppeared'],'-dpng','-r0') %dpng
   
    close all;
   % end %all trials per subject  
end % loop through participants

% plot all trials alligned to the  a specified Baseline
% p6 =figure(6);
% 
% plot((hPupilData(NrTrial).leftPupil_ValidSamples.signal.t),hPupilData(NrTrial).leftPupil_ValidSamples.signal.pupilDiameter_NoRemovedGap ,'r.')
% 
% 
% aim: percent signal change... sample(t1)/baseline as mean of a time period*100
% table with different Baseline-Periods
% 
% if preprocessing.analysis
%     great_averaged = zeros(offset_reference_index - onset_reference_index + 1,1);
%     all = 0;
%     if analysis.perGroup
%         for currentGroup = 1:groups.Count
%             
%             [indexes_current_group] = find(ismember(abs([trialData.refAngle_Sample] - [trialData.refAngle_M2S]),(groups(currentGroup))));
%             [indexes_current_group] = intersect(indexes_current_group,succesfull_indexes);
%             group_size = numel(indexes_current_group);
%             
%             if analisys.PlottingPupilSize
%                 temporal = arrayfun(@(x) doNothing(x.save_Pupil_w(onset_reference_index:offset_reference_index)), trialData(indexes_current_group), 'UniformOutput', false)';
%                 medianTimeSeries = 0;
%                 
%                 if ~isempty(temporal)
%                     converted2matrix = zeros(numel(temporal), offset_reference_index - onset_reference_index + 1);
%                                                 suma = zeros(numel(temporal{1}),1);
%                     all = all + 1;
%                     for t = 1: numel(temporal)
%                         converted2matrix(t,:) = temporal{t,1};
%                                                         suma = suma + temporal{t,1};
%                     end
%                     medianTimeSeries = median(converted2matrix);
%                     plotting(current_time_reference(onset_reference_index:offset_reference_index),median(converted2matrix),0,0, [groups(currentGroup) group_size], 'combined',0);
%                 end
%                 
%                 if analysis.total_average
%                     great_averaged = great_averaged + medianTimeSeries;
%                     if currentGroup == groups.Count
%                         great_averaged = great_averaged/all;
%                         plotting(current_time_reference(onset_reference_index:offset_reference_index),great_averaged,0,0, ['average' 'all'], 'combined',1);
%                     end
%                 end
%             end
%             if analysis.statistics
%                 if ~isempty(indexes_current_group)
%                     statistics(field_names, 101, [66 72], [onset_reference_index offset_reference_index], indexes_current_group, [groups(currentGroup) group_size]);
%                 end
%             end
%             
%         end
%         line([trialData(reference_index).Time_matchtosample_appeared_run trialData(reference_index).Time_matchtosample_appeared_run], get(gca,'YLim'),'Color','black','LineStyle','--');
%         
%     else
%         statistics(field_names, 101, [66 72], [onset_reference_index offset_reference_index], succesfull_indexes);
%         
%     end
%     
% end
% 


%% Analyze Processed Data:
% After processing the valid samples, the pupil size data must be analyzed
% per segment. These segments are defined in the hPupilData instance, see
% the RawFileModel documentationn for details. The analyzeSegments method
% analyzes either the sole available pupil, or both pupils and their mean,
% and horizontally concatenates the results, together with the segmentsData
% table, per PupilDataModel instance. The results are returned in a cell
% array.

% Get the analysis results, vertically concatenate them, and view the
% resulting master table:
results      = hPupilData.analyzeSegments();
PupilResults = vertcat(results{:});

% The code below open the PupilResults results in the variable viewer:
openvar PupilResults

% All pupil size data are now saved in the results table. If baseline
% correction is desired, the baseline segments can easily be subtracted
% from their concerning response epoch(s). To simplify this, it is handy to
% add a metadata column to the segmentData table containing an identifier
% labeling the epoch as a baseline. A baseline table can then be generated
% by extracting the rows featuring the baseline identifier from the results
% table. Similarly, a response table can be generated by extracting the
% desired response epochs. These baseline and response tables can then be
% combined using the 'join' table operation, producing a table with
% matched baseline and response columns.


%% Visualizing the Data:
% The PupilDataModel class can visualize its data via the plotData method,
% which can be called on a PupilDataModel array. Note that if the
% keepFilterData flag was set to true, the intermediate filering steps will
% be plotted as well, which may be slow when plotting multiple objects.
%
% The segments are visualized as stacked non-overlapping rectangles in the
% top axes. Their colors are random.
%
% *** You can click the segment rectangle to view its info ***.

% Plotting:
hPupilData.plotData;