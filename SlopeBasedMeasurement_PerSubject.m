function SlopeBasedMeasurement_PerSubject




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
		meta               = meta_conf(:) + meta_err(:);
		meta_pre           = meta_conf_pre + meta_err_pre;

%% Slop control

sum_ntrials_suc_pre_F1   = sum(Y_matrix_ntrials_suc_pre_F1,2);
sum_ntrials_suc_pre_F2   = sum(Y_matrix_ntrials_suc_pre_F2,2);
sum_ntrials_suc_pre_F3   = sum(Y_matrix_ntrials_suc_pre_F3,2);
sum_ntrials_unsuc_pre_F1 = sum(Y_matrix_ntrials_unsuc_pre_F1,2);
sum_ntrials_unsuc_pre_F2 = sum(Y_matrix_ntrials_unsuc_pre_F2,2);
sum_ntrials_unsuc_pre_F3 = sum(Y_matrix_ntrials_unsuc_pre_F3,2);

%% CORRECT AND INCORRECT FOR EACH DIFFICULTY LEVEL POST

Trials_post       = flipud(NumberTrials_postwagering);
Trials_suc_post   = flipud(NumberTrials_correctTrials_postwagering);
Trials_unsuc_post = flipud(NumberTrials_ErrorTrials_postwagering);

n_trials_suc_post_diff(:,:,1)   = flipud(NumberTrials_correctTrials_postwagering_diff1);
n_trials_unsuc_post_diff(:,:,1) = flipud(NumberTrials_ErrorTrials_postwagering_diff1);
n_trials_post_diff(:,:,1)       = n_trials_suc_post_diff(:,:,1)+n_trials_unsuc_post_diff(:,:,1);

n_trials_suc_post_diff(:,:,2)   = flipud(NumberTrials_correctTrials_postwagering_diff2);
n_trials_unsuc_post_diff(:,:,2) = flipud(NumberTrials_ErrorTrials_postwagering_diff2);
n_trials_post_diff(:,:,2)       = n_trials_suc_post_diff(:,:,2)+n_trials_unsuc_post_diff(:,:,2);

n_trials_suc_post_diff(:,:,3)   = flipud(NumberTrials_correctTrials_postwagering_diff3);
n_trials_unsuc_post_diff(:,:,3) = flipud(NumberTrials_ErrorTrials_postwagering_diff3);
n_trials_post_diff(:,:,3)       = n_trials_suc_post_diff(:,:,3)+n_trials_unsuc_post_diff(:,:,3);

n_trials_suc_post_diff(:,:,4)   = flipud(NumberTrials_correctTrials_postwagering_diff4);
n_trials_unsuc_post_diff(:,:,4) = flipud(NumberTrials_ErrorTrials_postwagering_diff4);
n_trials_post_diff(:,:,4)       = n_trials_suc_post_diff(:,:,4)+n_trials_unsuc_post_diff(:,:,4);

n_trials_suc_post_diff(:,:,5)   = flipud(NumberTrials_correctTrials_postwagering_diff5);
n_trials_unsuc_post_diff(:,:,5) = flipud(NumberTrials_ErrorTrials_postwagering_diff5);
n_trials_post_diff(:,:,5)       = n_trials_suc_post_diff(:,:,5)+n_trials_unsuc_post_diff(:,:,5);

for s = 1:length(D)
	for d = 1:difficulty_degrees
		p_trials_suc_post_diff(s,:,d)   = (n_trials_suc_post_diff(s,:,d)/sum(n_trials_suc_post_diff(s,:,d),2))*100;
		p_trials_unsuc_post_diff(s,:,d) = (n_trials_unsuc_post_diff(s,:,d)/sum(n_trials_unsuc_post_diff(s,:,d),2))*100;
		p_trials_post_diff(s,:,d)       = (n_trials_post_diff(s,:,d)/sum(n_trials_post_diff(s,:,d),2))*100;
	end
end

p_trials_unsuc_diff1_NaN = find(isnan(p_trials_unsuc_post_diff(:,:,1)));
p_trials_unsuc_post_diff(p_trials_unsuc_diff1_NaN) = 0;









%% CORRECT AND INCORRECT FOR EACH DIFFICULTY LEVEL PRE

Trials_pre       = flipud(NumberTrials_prewagering);
Trials_suc_pre   = flipud(NumberTrials_correctTrials_prewagering);
Trials_unsuc_pre = flipud(NumberTrials_ErrorTrials_prewagering);

n_trials_suc_pre_diff(:,:,1)   = flipud(NumberTrials_correctTrials_prewagering_diff1);
n_trials_unsuc_pre_diff(:,:,1) = flipud(NumberTrials_ErrorTrials_prewagering_diff1);
n_trials_pre_diff(:,:,1)       = n_trials_suc_pre_diff(:,:,1)+n_trials_unsuc_pre_diff(:,:,1);

n_trials_suc_pre_diff(:,:,2)   = flipud(NumberTrials_correctTrials_prewagering_diff2);
n_trials_unsuc_pre_diff(:,:,2) = flipud(NumberTrials_ErrorTrials_prewagering_diff2);
n_trials_pre_diff(:,:,2)       = n_trials_suc_pre_diff(:,:,2)+n_trials_unsuc_pre_diff(:,:,2);

n_trials_suc_pre_diff(:,:,3)   = flipud(NumberTrials_correctTrials_prewagering_diff3);
n_trials_unsuc_pre_diff(:,:,3) = flipud(NumberTrials_ErrorTrials_prewagering_diff3);
n_trials_pre_diff(:,:,3)       = n_trials_suc_pre_diff(:,:,3)+n_trials_unsuc_pre_diff(:,:,3);

n_trials_suc_pre_diff(:,:,4)   = flipud(NumberTrials_correctTrials_prewagering_diff4);
n_trials_unsuc_pre_diff(:,:,4) = flipud(NumberTrials_ErrorTrials_prewagering_diff4);
n_trials_pre_diff(:,:,4)       = n_trials_suc_pre_diff(:,:,4)+n_trials_unsuc_pre_diff(:,:,4);

n_trials_suc_pre_diff(:,:,5)   = flipud(NumberTrials_correctTrials_prewagering_diff5);
n_trials_unsuc_pre_diff(:,:,5) = flipud(NumberTrials_ErrorTrials_prewagering_diff5);
n_trials_pre_diff(:,:,5)       = n_trials_suc_pre_diff(:,:,5)+n_trials_unsuc_pre_diff(:,:,5);

for s = 1:length(D)
	for d = 1:difficulty_degrees
		p_trials_suc_pre_diff(s,:,d)   = (n_trials_suc_pre_diff(s,:,d)/sum(n_trials_suc_pre_diff(s,:,d),2))*100;
		p_trials_unsuc_pre_diff(s,:,d) = (n_trials_unsuc_pre_diff(s,:,d)/sum(n_trials_unsuc_pre_diff(s,:,d),2))*100;
		p_trials_pre_diff(s,:,d)       = (n_trials_pre_diff(s,:,d)/sum(n_trials_pre_diff(s,:,d),2))*100;
	end
end

p_trials_unsuc_diff1_NaN       = find(isnan(p_trials_unsuc_pre_diff(:,:,1)));
p_trials_unsuc_pre_diff(p_trials_unsuc_diff1_NaN) = 0;

%% Reaction time

for s = 1:length(D)
	for d = 1:number_wagers
		RT2_post((length(D)+1)-s,d,:)   = cell2mat(D(s).Mean_wagerRT_postwager_difficulty{d});
		RT1_post((length(D)+1)-s,d,:)   = cell2mat(D(s).Mean_wagerRT1_postwager_difficulty{d});
	end
end

%% Index for Confidence and Error Detection

pt_suc_post   = p_trials_suc_post_diff;
pt_unsuc_post = p_trials_unsuc_post_diff;
pt_post       = p_trials_post_diff;

pt_suc_pre    = p_trials_suc_pre_diff;
pt_unsuc_pre  = p_trials_unsuc_pre_diff;
pt_pre        = p_trials_pre_diff;

for s = 1:length(D)
	for d = 1:difficulty_degrees
		p_conf_post(s,:,d)      = polyfit(1:n_wagers,pt_suc_post(s,:,d),1);
		r_conf_post(s,:,d)      = p_conf_post(s,1,d) .* [1:n_wagers] + p_conf_post(s,2,d);
		p_err_post(s,:,d)       = polyfit(1:n_wagers,pt_unsuc_post(s,:,d),1);
		r_err_post(s,:,d)       = p_err_post(s,1,d) .* [1:n_wagers] + p_err_post(s,2,d);
		p_conf_pre(s,:,d)       = polyfit(1:n_wagers,pt_suc_pre(s,:,d),1);
		r_conf_pre(s,:,d)       = p_conf_pre(s,1,d) .* [1:n_wagers] + p_conf_pre(s,2,d);
		p_err_pre(s,:,d)        = polyfit(1:n_wagers,pt_unsuc_pre(s,:,d),1);
		r_err_pre(s,:,d)        = p_err_pre(s,1,d) .* [1:n_wagers] + p_err_pre(s,2,d);
		
		p_post(s,:,d)           = polyfit(1:n_wagers,pt_post(s,:,d),1);
		r_post(s,:,d)           = p_post(s,1,d) .* [1:n_wagers] + p_post(s,2,d);
		p_pre(s,:,d)            = polyfit(1:n_wagers,pt_pre(s,:,d),1);
		r_pre(s,:,d)            = p_pre(s,1,d) .* [1:n_wagers] + p_pre(s,2,d);
		
		r_conf_corrected(s,:,d) = r_conf_post(s,:,d) - r_conf_pre(s,:,d);
		r_err_corrected(s,:,d)  = r_err_post(s,:,d)  - r_err_pre(s,:,d);
		
		meta_conf(s,:,d)          = p_conf_post(s,1,d) - p_conf_pre(s,1,d);
		meta_err(s,:,d)           = -p_err_post(s,1,d) + p_err_pre(s,1,d);
		meta_conf_pre(s,:,d)      = p_conf_post(s,1,d) - p_conf_pre(s,1,d);
		meta_err_pre(s,:,d)       = -p_err_post(s,1,d) + p_err_pre(s,1,d);
		meta(s,:,d)               = meta_conf(s,:,d) + meta_err(s,:,d);
		meta_pre(s,:,d)           = meta_conf_pre(s,:,d) + meta_err_pre(s,:,d);
	end
end