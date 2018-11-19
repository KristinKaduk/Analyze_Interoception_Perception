function  out = metaD_PerSubject(n_wagers, wagering_or_controll_wagering_post, success,selected_Left, wager_choosen_post,selected_Right)
%trial.unsuccess   0: completed trials 
%trial.M2SPeriod_SampleOrNonSample  1:Sample and Match same, 2: diffrent
%trial.N_Selected    Response:1 Sample and Match are same
%trial.V_selected    Response:1 Sample and Match are different
%trial.wagering_or_controll_wagering_post  2: wagering post
%trial.wagering_or_controll_wagering_pre  1: wagering pre
%trial.wager_choosen  value of chosen wager 



%%meta_d
% input transformation for meta d calculations (all difficulty levels)
  %% % example: Inputs for calculating meta d'
%     %for d=1:number_wagers
      %D(x).idx_amount_postwager_success_S2_right{d}= find(D(x).wagering_or_controll_wagering1 == 1 & D(x).success== 1 & D(x).ms_side== 2 & D(x).wager_choosen == d);
      %end
%     
%     D(x).nR_S2 =[length(D(x).idx_amount_postwager_success_S1_left{6}) length(D(x).idx_amount_postwager_success_S1_left{5}) length(D(x).idx_amount_postwager_success_S1_left{4}) length(D(x).idx_amount_postwager_success_S1_left{3}) length(D(x).idx_amount_postwager_success_S1_left{2}) length(D(x).idx_amount_postwager_success_S1_left{1}) length(D(x).idx_amount_postwager_unsuccess_S1_left{1}) length(D(x).idx_amount_postwager_unsuccess_S1_left{2}) length(D(x).idx_amount_postwager_unsuccess_S1_left{3}) length(D(x).idx_amount_postwager_unsuccess_S1_left{4}) length(D(x).idx_amount_postwager_unsuccess_S1_left{5}) length(D(x).idx_amount_postwager_unsuccess_S1_left{6})];
%     D(x).nR_S1 =[length(D(x).idx_amount_postwager_success_S2_right{6}) length(D(x).idx_amount_postwager_success_S2_right{5}) length(D(x).idx_amount_postwager_success_S2_right{4}) length(D(x).idx_amount_postwager_success_S2_right{3}) length(D(x).idx_amount_postwager_success_S2_right{2}) length(D(x).idx_amount_postwager_success_S2_right{1}) length(D(x).idx_amount_postwager_unsuccess_S2_right{1}) length(D(x).idx_amount_postwager_unsuccess_S2_right{2}) length(D(x).idx_amount_postwager_unsuccess_S2_right{3}) length(D(x).idx_amount_postwager_unsuccess_S2_right{4}) length(D(x).idx_amount_postwager_unsuccess_S2_right{5}) length(D(x).idx_amount_postwager_unsuccess_S2_right{6})];
%     D(x).nR_S2 =fliplr(D(x).nR_S2);
%     nR_S1      =D(x).nR_S1;
%     nR_S2      =D(x).nR_S2;
%     out        = type2_SDT_SSE(nR_S1, nR_S2);
%     D(x).out   = out;
%     D(x).type2_fit.obs_FAR2_rS1                              =D(x).out.type2_fit.obs_FAR2_rS1;
%     D(x).type2_fit.obs_FAR2_rS2                              =D(x).out.type2_fit.obs_FAR2_rS2;
%     

%%  Nr.of Trials per Wager seprated for Successfull or unsuccessful trials and Stimulus either left or right
      for d=1:n_wagers
        idx_amount_postwager_success_S1Left_respondedLeft{d}= find( (wagering_or_controll_wagering_post == 3  |   wagering_or_controll_wagering_post == 2) & success== 1 & selected_Left == 1 & wager_choosen_post == d);
        idx_amount_postwager_unsuccess_S1Left_respondedRight{d}= find((wagering_or_controll_wagering_post == 3  |   wagering_or_controll_wagering_post == 2) & success== 0 & selected_Left == 1 & wager_choosen_post == d);
        idx_amount_postwager_success_S2Right_respondedRight{d}= find((wagering_or_controll_wagering_post == 3  |   wagering_or_controll_wagering_post == 2) & success== 1 & selected_Right == 1 & wager_choosen_post == d);
        idx_amount_postwager_unsuccess_S2Right_respondedLeft{d}= find((wagering_or_controll_wagering_post == 3  |   wagering_or_controll_wagering_post == 2) & success== 0 & selected_Right == 1 & wager_choosen_post == d);
      end    
nR_S1=zeros(1,n_wagers*2);
nR_S2=zeros(1,n_wagers*2);
% why is S1 and s2 exchanged?
% successful trials unsuccessful trials
%nR_S1 = [100 2 100 100 50 50 150 20 30 1  1 1]
%nR_S2 = [1 1 1 30 20 150 0 50 100 100 2 0]

    nR_S1 =[length(idx_amount_postwager_success_S1Left_respondedLeft{6}) length(idx_amount_postwager_success_S1Left_respondedLeft{5}) length(idx_amount_postwager_success_S1Left_respondedLeft{4}) length(idx_amount_postwager_success_S1Left_respondedLeft{3}) length(idx_amount_postwager_success_S1Left_respondedLeft{2}) length(idx_amount_postwager_success_S1Left_respondedLeft{1}) length(idx_amount_postwager_unsuccess_S1Left_respondedRight{1}) length(idx_amount_postwager_unsuccess_S1Left_respondedRight{2}) length(idx_amount_postwager_unsuccess_S1Left_respondedRight{3}) length(idx_amount_postwager_unsuccess_S1Left_respondedRight{4}) length(idx_amount_postwager_unsuccess_S1Left_respondedRight{5}) length(idx_amount_postwager_unsuccess_S1Left_respondedRight{6})];
    nR_S2 =[length(idx_amount_postwager_unsuccess_S2Right_respondedLeft{6}) length(idx_amount_postwager_unsuccess_S2Right_respondedLeft{5}) length(idx_amount_postwager_unsuccess_S2Right_respondedLeft{4}) length(idx_amount_postwager_unsuccess_S2Right_respondedLeft{3}) length(idx_amount_postwager_unsuccess_S2Right_respondedLeft{2}) length(idx_amount_postwager_unsuccess_S2Right_respondedLeft{1}) length(idx_amount_postwager_success_S2Right_respondedRight{1}) length(idx_amount_postwager_success_S2Right_respondedRight{2}) length(idx_amount_postwager_success_S2Right_respondedRight{3}) length(idx_amount_postwager_success_S2Right_respondedRight{4}) length(idx_amount_postwager_success_S2Right_respondedRight{5}) length(idx_amount_postwager_success_S2Right_respondedRight{6})];
    suc = nR_S2(7:12) + flip(nR_S1(1:6));
    nR_S2(7:12)= suc;
    nR_S1(1:6) = flip(suc) ;
    unsuc = nR_S2(1:6) + flip(nR_S1(7:12));
    nR_S2(1:6) = unsuc;
    nR_S1(7:12) = flip(unsuc);
    out        = type2_SDT_SSE(nR_S1, nR_S2);
    type2_fit.obs_FAR2_rS1                              =out.type2_fit.obs_FAR2_rS1;
    type2_fit.obs_FAR2_rS2                              =out.type2_fit.obs_FAR2_rS2;  

% 
% 
% IndS1S1=find((unsuccess==0)&(control_post==2)&(Stimuli_LeftOrRight==2)&(selected_Right==1)); %index for Match and Sample: different Response: different
% IndS1S2=find((unsuccess==0)&(control_post==2)&(Stimuli_LeftOrRight==2)&(selected_Left==1)); %index for Match and Sample: differnt  Response: same
% IndS2S1=find((unsuccess==0)&(control_post==2)&(Stimuli_LeftOrRight==1)&(selected_Right==1)); %index for Match and Sample: same      Response: different
% IndS2S2=find((unsuccess==0)&(control_post==2)&(Stimuli_LeftOrRight==1)&(selected_Left==1)); %index for Match and Sample: same      Response: same
% 
% wagers_S1S1=zeros(1,n_wagers);
% wagers_S1S2=zeros(1,n_wagers);
% wagers_S2S1=zeros(1,n_wagers);
% wagers_S2S2=zeros(1,n_wagers);
% 
% for i=1:SETTINGS.NrOfWagers
%     wagers_S1S1(i)= sum(wager_choosen_post(IndS1S1)==i);
%     wagers_S1S2(i)= sum(wager_choosen_post(IndS1S2)==i);
%     wagers_S2S1(i)= sum(wager_choosen_post(IndS2S1)==i);
%     wagers_S2S2(i)= sum(wager_choosen_post(IndS2S2)==i);
% end
% nR_S1= cat(2,fliplr(wagers_S1S1),wagers_S1S2);
% nR_S2= cat(2,fliplr(wagers_S2S1),wagers_S2S2);
% 
% %meta-d function Maniscalo 
% out=[];
% out = type2_SDT_SSE(nR_S1, nR_S2);
% 
% %%slope-based measurements (all difficulty levels)
% 
% 
% %wager proportions ordered from lowest to highest wager
% %pre: 
% Ind_pre_correct=find((unsuccess==0)&(control_pre==1)&(((M2SPeriod_SampleOrNonSample==2)&(selected_Diff==1))|((M2SPeriod_SampleOrNonSample==1)&(selected_Same==1))));
% Ind_pre_incorrect=find((unsuccess==0)&(control_pre==1)&(((M2SPeriod_SampleOrNonSample==2)&(selected_Same==1))|((M2SPeriod_SampleOrNonSample==1)&(selected_Diff==1))));
% for i=1:n_wagers
%     wager_correct_pre(i)=sum(wager_choosen_pre(Ind_pre_correct)==i);
%     wager_incorrect_pre(i)=sum(wager_choosen_pre(Ind_pre_incorrect)==i);
% end
% wager_prop_correct_pre=(wager_correct_pre./sum(wager_correct_pre))*100;
% wager_prop_incorrect_pre=(wager_incorrect_pre./sum(wager_incorrect_pre))*100;
% % post: 
% wager_prop_correct_post=((wagers_S1S1+wagers_S2S2)./sum(wagers_S1S1+wagers_S2S2))*100;
% wager_prop_incorrect_post=((wagers_S1S2+wagers_S2S1)./sum(wagers_S1S2+wagers_S2S1))*100;
% 
% %linear fit
% wagers=[1 2 3 4 5 6];
% %pre
% b_correct_pre = wagers'\wager_prop_correct_pre';
% b_incorrect_pre=wagers'\wager_prop_incorrect_pre';
% %post
% b_correct_post = wagers'\wager_prop_correct_post';
% b_incorrect_post=wagers'\wager_prop_incorrect_post';
% 
% %readouts
% readout_correct=(atan(b_correct_post)-atan(b_correct_pre))*180/pi;
% readout_incorrect=-((atan(b_incorrect_post)-atan(b_incorrect_pre))*180/pi);
% slope_metacognitive_ability=readout_correct+readout_incorrect;
% 
% %%Meessen calculation (all difficulty levels)
% n_trials_post=length(find((unsuccess==0)&(control_post==2)));
% 
% percentage_correct=((sum(wagers_S1S1)+sum(wagers_S2S2))/n_trials_post)*100;
% performance=percentage_correct*ones(1,n_trials_post);
% confidence=wager_choosen_post(sort(cat(2,IndS1S1,IndS2S2,IndS1S2,IndS2S1))); 
% 
% meessen_metacognitive_ability=(sum((performance-rescalewag(confidence,1,6,0,100)).^2))/n_trials_post;
% 
% 
% 
% 


%% SUBFUNCTIONS
function out = type2_SDT_SSE(nR_S1, nR_S2)

% out = type2_SDT(nR_S1, nR_S2)
%
% Given data from an experiment where an observer discriminates between two
% stimulus alternatives on every trial and provides confidence ratings,
% provides a type 2 SDT analysis of the data.
%
% The function does a standard type 1 SDT analysis on the raw behavioral
% data and then does a type 2 SDT analysis using the function fit_meta_d
% with d_min = -5, d_grain = .01, d_max = 5
%
% INPUTS
%
% nR_S1, nR_S2
%
% nR_S1 and nR_S2 are vectors containing the total number of responses in
% each response category, conditional on presentation of S1 and S2.
% size of each array is 2*nRatings, where each element corresponds to a
% count of responses in each response category. Response categories are
% ordered as follows:
% highest conf "S1" ... lowest conf "S1", lowest conf "S2", ... highest conf "S2"
%
% e.g. if nR_S1 = [100 50 20 10 5 1], then when stimulus S1 was
% presented, the subject had the following response counts:
% responded S1, rating=3 : 100 times correct
% responded S1, rating=2 : 50 times correct
% responded S1, rating=1 : 20 times correct
% responded S2, rating=1 : 10 times incorrect
% responded S2, rating=2 : 5 times incorrect
% responded S2, rating=3 : 1 time incorrect
%
% The ordering of response / rating counts for S2 should be the same as it
% is for S1. e.g. if nR_S2 = [3 7 8 12 27 89], then when stimulus S2 was
% presented, the subject had the following response counts:
% responded S1, rating=3 : 3 times incorrect
% responded S1, rating=2 : 7 times
% responded S1, rating=1 : 8 times
% responded S2, rating=1 : 12 times correct
% responded S2, rating=2 : 27 times
% responded S2, rating=3 : 89 times
%
%
% OUTPUTS
%
% out.d_a       : d_a for input data. If s=1, d_a = d'
% out.meta_d_a  : meta_d_a for input data
% out.M_ratio   : meta_d_a / d_a; measure of metacognitive efficiency
% out.M_diff    : meta_d_a - d_a; measure of metacognitive efficiency
% out.s         : ratio of evidence distribution standard deviations assumed for the analysis.
% out.type2_fit : output of fit_meta_d_SSE for the type 2 SDT fit.

% 10/14/14 - removed support for unequal variance SDT model
%          - removed support for input data format describing
%            trial-by-trial outcomes (function "trials2counts" is now
%            provided on http://www.columbia.edu/~bsm2105/type2sdt/ to
%            convert trial data to response count data)
% 9/7/10 - bm - wrote it

%% parse inputs

% check for valid inputs
if ~( length(nR_S1) == length(nR_S2) && mod(length(nR_S2),2) == 0 )
    error('nR_S1 and nR_S2 must be the same length and have an even number of elements')
end

nRatings = length(nR_S1) / 2;


%% standard SDT analysis

  
%(?-KK) Why only one part of the input is used?
HR1  = sum(nR_S2(nRatings+1:end)) / sum(nR_S2);
FAR1 = sum(nR_S1(nRatings+1:end)) / sum(nR_S1);

for i=2:2*nRatings
    ratingHRs(i-1)  = sum(nR_S2(i:end)) / sum(nR_S2);
    ratingFARs(i-1) = sum(nR_S1(i:end)) / sum(nR_S1);
end

% if equalVariance
%     s = 1;
% else
%     p = polyfit(norminv(ratingFARs), norminv(ratingHRs), 1);
%     s = p(1);
% end

s = 1;

% 1 or 0 for FAR1 or HR1 -> 
if isnan(FAR1)==1; %~pFA ||
    FAR1 = 0.0001;
elseif FAR1==0
    FAR1= 0.0001;
elseif FAR1==1
    FAR1 = 0.9999;
end
if isnan(HR1)==1
    HR1= 0.0001;
elseif HR1==0; %~pHit ||
    HR1 = 0.0001;
elseif HR1==1
    HR1 = 0.9999;
end

% d' and c in terms of S1 distribution standard deviation units
d_1 = (1/s)*norminv(HR1) - norminv(FAR1);
c_1 = (-1/(1+s)) * (norminv(HR1)+norminv(FAR1));
cprime = c_1 / d_1;


%% type 2 SDT analysis

% get type 2 HR and FAR for S1 responses
for i = 1 : nRatings-1
    HR2_rS1(i)  = sum(nR_S1(1:i)) / sum(nR_S1(1:nRatings));
    FAR2_rS1(i) = sum(nR_S2(1:i)) / sum(nR_S2(1:nRatings));
    % FA_befor1 = FAR2_rS1(i)
    if isnan(FAR2_rS1(i));
        FAR2_rS1(i) = 0.0001;
    elseif FAR2_rS1(i)<0.0001;
        FAR2_rS1(i)= 0.0001;
    elseif FAR2_rS1(i)>0.9999;
        FAR2_rS1(i) = 0.9999;
    end
    % FA_after1 = FAR2_rS1(i)
end

% get type 2 HR and FAR for S2 responses
for i = nRatings+2 : 2*nRatings
    HR2_rS2(i - (nRatings+2) + 1)  = sum(nR_S2(i:end)) / sum(nR_S2(nRatings+1:end));
    FAR2_rS2(i - (nRatings+2) + 1) = sum(nR_S1(i:end)) / sum(nR_S1(nRatings+1:end));
    % FA_before2 = FAR2_rS2(i - (nRatings+2) + 1)
    if isnan(FAR2_rS2(i - (nRatings+2) + 1));
        FAR2_rS2(i - (nRatings+2) + 1)= 0.0001;
    elseif FAR2_rS2(i - (nRatings+2) + 1)<0.0001;
        FAR2_rS2(i - (nRatings+2) + 1) = 0.0001;
    elseif FAR2_rS2(i - (nRatings+2) + 1)>0.9999;
        FAR2_rS2(i - (nRatings+2) + 1) = 0.9999;
    end
    % FA_after2 = FAR2_rS2(i - (nRatings+2) + 1)
end

d_min = -5;
d_grain = .01;
d_max = 5;

fit = fit_meta_d_SSE(HR2_rS1, FAR2_rS1, HR2_rS2, FAR2_rS2, cprime, s, d_min, d_max, d_grain);

%% package output
out.da        = d_1 * s * sqrt(2/(1+s^2));
out.meta_da   = fit.meta_d1 * s * sqrt(2/(1+s^2));
out.M_ratio   = out.meta_da / out.da;
out.M_diff    = out.meta_da - out.da;
out.s         = s;
out.type2_fit = fit;
out.Typ1_criterion = c_1;
out.Typ1_cprime = cprime;

function fit = fit_meta_d_SSE(obs_HR2_rS1, obs_FAR2_rS1, obs_HR2_rS2, obs_FAR2_rS2, cprime, s, d_min, d_max, d_grain)
% fit = fit_metad_SSE(obs_HR2_rS1, obs_FAR2_rS1, obs_HR2_rS2, obs_FAR2_rS2,
%                     cprime, s, d_min, d_max, d_grain)
%
% Given response-conditional type 2 hit rates and type 2 false alarm rates,
% as well as the empirically estimated relative criterion cprime = c / d',
% use a signal detection theory model to estimate meta-d',
% the value of d' that would have been expected to generate the observed
% type 2 data.
%
% Estimation is done by testing different values for meta-d' to see what
% value gives the best fit to the observed data (best fit = minimizes sum
% squared error between observed and expected type 2 HR and type 2 FAR).
%
%
% required inputs
% ---------------
% obs_HR2_rS1, obs_FAR2_rS1 :
% Arrays of observed type 2 hit rate (p(high confidence|correct "S1" response))
% and type 2 false alarm rate (p(high confidence|incorrect "S1" response)).
% The size of each array is N-1, where N is the number of options for
% rating confidence. So for instance, if you use a confidence scale with
% 4 levels of confidence, these arrays should contain 3 elements each.
% The i_th element corresponds to the type 2 HR and type 2 FAR found by
% considering all confidence ratings of i+1 or higher to be "high
% confidence".
% The i_th element of obs_HR2_rS1 must correspond to the i_th element of
% obs_FAR2_rS1. Otherwise, ordering of the data is not important.
%
% obs_HR2_rS2, obs_FAR2_rS2 : same as above, for "S2" responses
%
% cprime :
% The relative type 1 criterion.
%
% c' = c / d', where
%
% d' = z(type 1 HR) - z(type 1 FAR)
% c = -.5 * [z(type 1 HR) + z(type 1 FAR)]
%
% and z is the inverse of the normal cumulative distribution function.
%
% If s != 1, specify c' in units of the S1 distribution, as follows.
%
% d' = (1/s)*z(type 1 HR) - z(type 1 FAR)
% c = [ -1 / (1+s) ] * [z(type 1 HR) + z(type 1 FAR)]
%
% optional inputs
% ---------------
% s :
% The ratio of the standard deviations of the evidence distributions for
% stimulus classes S1 and S2. Can be estimated from rating data.
% If unspecified, s = 1
%
% d_min :
% The minimimum value for meta-d' that will be tested.
% If unspecified, d_min = -5
%
% d_max :
% The maximum value for meta-d' that will be tested.
% If unspecified, d_max = 5
%
% d_grain :
% The step size used in testing values of meta-d'.
% If unspecified, d_grain = .01
%
%
% output
% ------
% fit.meta_d :
% meta_d' value that minimizes the SSE between observed and expected type 2
% data. If s != 1, meta_d' is specified in units of the S1 distribution.
%
% fit.meta_c :
% The value of type 1 criterion c used in conjunction with meta_d'.
% meta_c / meta_d = cprime, the constant type 1 criterion specified in the
% input. If s != 1, meta-c is specified in units of the S1 distribution.
%
% fit.s :
% The value of s used in the type 2 data fitting, where s = sd(S1) / sd(S2)
%
% fit.t2c_rS1 :
% Values for the type 2 criteria that, along with meta-d' and c', provide
% the best fit for type 2 data for "S1" responses
%
% fit.t2c_rS2 :
% Likewise, for "S2" responses
%
% fit.SSE :
% Sum of squared errors between observed and expected type 2 data
%
% fit.est_HR2_rS1 :
% The type 2 hit rates for "S1" responses expected from meta_d, meta_c, s,
% and t2c_rS1
%
% fit.obs_HR2_rS1 :
% Empirically observed type 2 hit rates for "S1" responses
%
% fit.est_FAR2_rS1, fit.obs_FAR2_rS1, fit.est_HR2_rS2, ...
% Likewise as above, for expected and observed type 2 FAR for "S1"
% responses and type 2 HR and type 2 FAR for "S2" responses

% 9/7/10 - bm - wrote it


%% initialize optional inputs
if ~exist('s','var') || isempty(s), s=1; end
if ~exist('d_min','var') || isempty(d_min), d_min = -5; end
if ~exist('d_max','var') || isempty(d_max), d_max = 5; end
if ~exist('d_grain','var')  || isempty(d_grain), d_grain = .01; end



%% initialize analysis
nRatings = length(obs_HR2_rS1);

ds = d_min : d_grain : d_max;

SSEmin = Inf;
meta_d = [];
meta_c = [];
t2c_rS1 = [];
t2c_rS2 = [];
est_HR2_rS1 = [];
est_FAR2_rS1 = [];
est_HR2_rS2 = [];
est_FAR2_rS2 = [];


%% search for meta-d' that minimizes type 2 SSE
for i=1:length(ds)
    
    % initialize parameters for current level of meta-d'
    d = ds(i);
    c = cprime * d;
    
    S1mu = -d/2;
    S2mu = d/2;
    S1sd = 1;
    S2sd = 1/s;
    
    x = S1mu - 5*max([S1sd S2sd]) : .001 : S2mu + 5*max([S1sd S2sd]);
    [diff c_ind] = min(abs(x-c));
    
    HRs = 1-normcdf(x,S2mu,S2sd);
    FARs = 1-normcdf(x,S1mu,S1sd);
    
    % fit type 2 data for S1 responses
    est_HR2s_rS1 = (1-FARs(1:c_ind)) / (1-FARs(c_ind));
    est_FAR2s_rS1 = (1-HRs(1:c_ind)) / (1-HRs(c_ind));
    
    for n=1:nRatings
        SSE = (est_HR2s_rS1 - obs_HR2_rS1(n)).^2 + (est_FAR2s_rS1 - obs_FAR2_rS1(n)).^2;
        [SSE_rS1(n) rS1_ind(n)] = min(SSE);
    end
    
    % fit type 2 data for S2 responses
    est_HR2s_rS2 = HRs(c_ind:end) / HRs(c_ind);
    est_FAR2s_rS2 = FARs(c_ind:end) / FARs(c_ind);
    
    for n=1:nRatings
        SSE = (est_HR2s_rS2 - obs_HR2_rS2(n)).^2 + (est_FAR2s_rS2 - obs_FAR2_rS2(n)).^2;
        [SSE_rS2(n) rS2_ind(n)] = min(SSE);
    end
    
    % update analysis
    SSEtot = sum(SSE_rS1) + sum(SSE_rS2);
    if SSEtot < SSEmin
        SSEmin = SSEtot;
        meta_d = d;
        meta_c = c;
        t2c_rS1 = x(rS1_ind);
        t2c_rS2 = x(c_ind + rS2_ind - 1);
        est_HR2_rS1  = est_HR2s_rS1(rS1_ind);
        est_FAR2_rS1 = est_FAR2s_rS1(rS1_ind);
        est_HR2_rS2  = est_HR2s_rS2(rS2_ind);
        est_FAR2_rS2 = est_FAR2s_rS2(rS2_ind);
    end
    
end

%% package output
fit.meta_d1  = meta_d;
fit.meta_c1  = meta_c;
fit.s        = s;
fit.t2c1_rS1 = t2c_rS1;
fit.t2c1_rS2 = t2c_rS2;

fit.SSE = SSEmin;

fit.est_HR2_rS1  = est_HR2_rS1;
fit.obs_HR2_rS1  = obs_HR2_rS1;

fit.est_FAR2_rS1 = est_FAR2_rS1;
fit.obs_FAR2_rS1 = obs_FAR2_rS1;

fit.est_HR2_rS2  = est_HR2_rS2;
fit.obs_HR2_rS2  = obs_HR2_rS2;

fit.est_FAR2_rS2 = est_FAR2_rS2;
fit.obs_FAR2_rS2 = obs_FAR2_rS2;

function [dpri,ccrit] = dprime(pHit,pFA,nTarget,nDistract)
% DPRIME  --  Calculate sensitivity index from signal-detection theory
%
%  USE:
%  dvalue = dprime(pHit,pFA) returns the sensitivity index, using a
%  standard value for correction of rates of 0 and 1.
%  dvalue = dprime(pHit,pFA,nTarget,nDistract) uses the given numbers of
%  targets and distractor trials value for correction.
%  [dvalue,cvalue] = dprime(pHit,pFA) returns also the individual bias
%  criterion.
%
%  Coding by Martin B?ckmann-Barthel, Otto-von-Guericke-Universit?t Magdeburg
%  The sensitivity index d' measures the distance between the signal and
%  the noise means in standard deviation units. c is the distance of the
%  bias criterion %  from the point where neither response is favored, also
%  in standard units. Positive c values indicate a bias towards high hit
%  rates and false alarm rates
% 
%  Obligatory input variables: 
%  pHit, pFA - hit rate and false alarm rate (max: 1.0, min 0.0)
%  Optional variables: nTarget, nDistract - number of targets and
%  distractor trials (needed for correction of perfect responses)
%  Perfect hit rates (pHit = 1) and FA rates (pFA = 0) are corrected 
%  by -1/(2 nTarget) and +1/(2 nDistract), respectively, if provided
%  cf. Stanislaw H, Todorov N, Behav Res Meth (1999) 31, 137-149, "1/2N rule"
%
%--------------------------------------------------------------------------
%
%  $Revision: 2.0 $  $Date: 2014-07-29 $
%  Handling of perfect response
%  $Revision: 2.01 $  $Date: 2016-08-31 $
%  Help edit
%  $Revision: 3.0 $  $Date: 2017-02-15 $
%  Criterium c (response bias) added
%  $Revision: 3.01 $  $Date: 2017-12-09 $
%  Help edit, input error check
%-- Replace rates equalling zero or one
if nargin < 4 % number of distractor presentations
    nDistract = 1e8; % if not specified, take a very high number
end
if nargin < 3 % number of target presentations
    nTarget = 1e8; % if not specified, take a very high number
end
if pHit > 1 | pFA > 1 
    error('Meaningless probabilities. (Do NOT enter percentage values!)');
end % if
if pHit < 0 | pFA < 0 
    error('Meaningless negative probabilities.');
end % if
if pHit > 0
    pHit = min(pHit,1-.5/nTarget);
else pHit = .5/nTarget;
end % if pHit=0
if pFA < 1
    pFA = max(pFA,.5/nTarget);
else pFA=1-.5/nTarget;
end % if pFA=0
%-- Convert to Z scores, no error checking
zHit = -sqrt(2).*erfcinv(2*pHit);
zFA = -sqrt(2).*erfcinv(2*pFA);
%-- Calculate d-prime
dpri = zHit - zFA ;
%-- Calculate d-prime and bias
dpri = zHit - zFA ;
if nargout > 1
    ccrit = -.5*(zHit + zFA);
end % if nargout
%  Return DPRIME

