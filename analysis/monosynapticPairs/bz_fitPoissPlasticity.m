function bz_fitPoissPlasticity(basepath, preID,postID,varargin)
% This function takes a pair of neurons and fits a Poisson GLM to account
% for long term monosynaptic coupling. This framework exploits a bin-by-bin baseline
% offset that reflects the postsynaptic spiking binned at a timescale
% slower than the fine timescale spike-spike coupling.
%
% INPUT
% basepath = directory where at least one of two files are located (1 requires) a buzcode
%            file [basename '.spikes.cellinfo.mat'] and (2 optional) a file
%            with a struct 'stim' with field 'ts' which gives a cell array
%            of stimulation times per location ([start stop])
% preID    = index of presynaptic neuron, e.g. spikes.times{preID}
% postID   = index of presynaptic neuron, e.g. spikes.times{postID}
%
% OUTPUT
% resaves inputs: 'basepath','preID', 'postID'
% synParamsT       = spk trans prob basis weights
% test_nll         = negative log likelihood on test set during cross validation
% timescales       = basis separation for spike trans to vary over (must be
%                    enough spikes for an estimate
% ltc_t            = spike transmission timeseries
%
%    'basepath', 'preID', 'postID', 'pre', 'post', 'synParamsT', 'test_nll', 'timescales', 'down_dt', 'ts_lastelement', ...
%    'pre_binned_inds', 'post_binned_inds', 'ltc_t', 'lambda','argmax', 'delta', 'binsize', 'dt', 'fitting_duration'
% 
% NOTES
%  - this function does not verify whether the passed neuron IDs are
%    monosynaptically connected
%    
% Roman Huszar & Sam McKenzie, August 2019

%% STEP 0 PROCESS INPUTS

tic

% Parse / save user inputs
p = inputParser;

addParameter(p, 'outpath', fullfile(basepath,'monosyn_plast'), @isstr)
addParameter(p, 'splitsize', 100, @isnumeric);
addParameter(p, 'delta', 0.015, @isnumeric);
addParameter(p, 'binsize', 0.015, @isnumeric);
addParameter(p, 'num_cv_iter', 10, @isnumeric);
addParameter(p, 'cv_fit_iter', 50, @isnumeric);
addParameter(p, 'timescales', [400 600 800 1000], @isvector);
addParameter(p, 'fit_iter', 1000, @isnumeric);
addParameter(p, 'dt', 0.0008, @isnumeric);
addParameter(p, 'down_dt', 0.1, @isnumeric);

parse(p,varargin{:})

outpath = p.Results.outpath;
splitsize = p.Results.splitsize;            % for crossvalidating
delta = p.Results.delta;                    % ( bz_SpktToSpkmat parameter )
binsize = p.Results.binsize;                % ( bz_SpktToSpkmat parameter )
num_cv_iter = p.Results.num_cv_iter;        % number of crossvalidations
cv_fit_iter = p.Results.cv_fit_iter;        % number of crosssvalidation fit iterations
timescales = p.Results.timescales;          % basis timescales
fit_iter = p.Results.fit_iter;              % number of fit iterations
dt = p.Results.dt;                          % binsize for spike train binning
down_dt = p.Results.down_dt;                % for downsampling excess synchrony estimates                                                                                                

% Create path to output
if ~exist(outpath, 'dir')
    mkdir(outpath)
end

basename = bz_BasenameFromBasepath(basepath);
pairID = ['_p' num2str(preID) '_' num2str(postID)];
outfile = fullfile(outpath,[basename pairID]);
% Give user feedback about what is being fit
fprintf('%s, %d to %d \n', basepath, preID, postID)

%% STEP 1 PROCESS SPIKES
    
% Load spikes
spikes = bz_GetSpikes('basepath', basepath);

% Pull out the full spike trains
post = spikes.times{postID}; 
pre = spikes.times{preID};

% Bin spikes
latest_spk = max(cellfun(@max, spikes.times));
ts = [0:dt:latest_spk]';
post_binned = [histcounts(post, [ts ; ts(end)+dt])]';
post_binned(post_binned > 1) = 1;
post = find(post_binned)*dt;
pre_binned = [histcounts(pre, [ts ; ts(end)+dt])]';
pre_binned(pre_binned > 1) = 1;
pre = find(pre_binned)*dt;

%% STEP 2 GENERATE COUPLING COVARIATE

% The latency to postsynaptic spike is estimated from the full CCG - we
% take the short-latency positive lag bin with the largest number of spikes
spk_times = [pre ; post];
groups = [ones(size(pre)) ; 2*ones(size(post))];

[spk_times,b] = sort(spk_times);
groups = groups(b);
[ccg,t] = CCG(spk_times,groups, 'binsize', dt, 'duration', 0.2);
ccg = ccg(:,1,2); 
ccg([1:126 135:end]) = 0;  % Ignore counts outside the relevant window
[ ~,ind ] = max(ccg);
argmax = t(ind) / dt;
% Generate the coupling covariate - an indicator function that is 1 at the
% monosynaptic transmission time lag following each presynaptic spike
pre_spk_inds = find(pre_binned);
couple_indicator = zeros(size(post_binned));
couple_indicator(pre_spk_inds+argmax) = 1;
if length(couple_indicator) > length(post_binned)
    couple_indicator(length(post_binned)+1:end) = [];
end
% Keep track of the presynaptic spike indices associated with bins where 
% the coupling indicator is 1
couple_indicator_inds = couple_indicator;
couple_indicator_inds(couple_indicator_inds == 1) = 1:sum(couple_indicator);

%% STEP 3 PREPARE ALL THE REGRESSORS

% Generate even/odd splits indicator for crossvalidation
% We take every other 100 second chunk of data for training the model. 
% The remaining data chunks are used for testing.
step = round(splitsize/dt);
splits_indicator = zeros(size(post_binned));
for ii = 0:ceil(length(splits_indicator)/step)-1
    if mod(ii,2) == 0
        if (ii+1)*step <= length(splits_indicator)
            splits_indicator(ii*step+1:(ii+1)*step) = 1;
        else
            splits_indicator(ii*step+1:end) = 1;
        end
    else
        continue;
    end
    
end
splits_indicator = logical(splits_indicator);

% Bin the postsynaptic spike train to get the slow rate
post_rate = bz_SpktToSpkmat({post}, 'binsize', binsize, 'dt', delta);
% Interpolate at the dt resolution to get an estimate of the coarse
% baseline at every time step
coarsened_rate = interp1(post_rate.timestamps, post_rate.data / post_rate.binsize, ts,'linear','extrap');
coarsened_rate(coarsened_rate <= 0) = 1e-6;
% We take the log of the coarsened rate, so that when it's passed through an
% exponential in the GLM, you get 'coarsened_rate' as your predicted Poisson rate
X_cPost = log(double(coarsened_rate));
clear smoothed_rate post_rate

% Preparing options for constrained mininimization
opts.dt = dt;
options = optimoptions('fmincon',...
    'Algorithm','interior-point',...    
    'OptimalityTolerance',1e-5,...     
    'StepTolerance',1e-8,...               
    'ConstraintTolerance',1e-5,...          
    'SpecifyObjectiveGradient', true);


%% STEP 4 CROSS-VALIDATION
%  Find the best spline basis timescale

fprintf('CV to find the right timescale...\n');
iter = 1;
% Store the test log likelihood for each timescale iteration
test_nll = {};
for timescale = timescales

    fprintf('Timescale %d\n', timescale);
    hyper_params_tmp.XTimeSpan = timescale;
    % Obtain basis - center them on presynaptic spike times
    X_couple = calculate_X({pre},hyper_params_tmp);
    if size(X_couple,1) > sum(couple_indicator)
        X_couple(sum(couple_indicator)+1:end,:) = [];
    end
    % Prepare training data regressors 
    coupling_bins_train = couple_indicator & splits_indicator;   % Coupling bins within training set
    opts.rate = X_cPost(coupling_bins_train,:);
    X_odd_temp =  X_couple( couple_indicator_inds(coupling_bins_train),: );
    nparams = size(X_odd_temp,2);

    % Handle for the cost function
    costfun = @(x) nll_poissPlasticity(x,X_odd_temp,post_binned(coupling_bins_train),opts);
    clear X_odd_temp coupling_bins_train

    % For each timescale, we fit model and compute the test likelihood
    for kk = 1:num_cv_iter

        % Set the number of iterations for initialization
        options.MaxIterations = cv_fit_iter;
        options.Display = 'off';
        
        % Initialize parameters
        bta_init = rand(nparams, 1);
        % Constrained optimization problem - we enforce nonnegative basis parameters
        [synParamsT,~,~,~] = fmincon(costfun, bta_init, [],[],[],[],...
                zeros(1,nparams),... 
                inf(1,nparams),...
                [],options);

        % Get negative log likelihood of the test set
        coupling_bins_test = couple_indicator & ~splits_indicator;   % Coupling bins within indicator
        opts.rate = X_cPost(coupling_bins_test,:);
        X_even_temp = X_couple( couple_indicator_inds(coupling_bins_test),: );
        [ f_test,~,~ ] = nll_poissPlasticity(synParamsT,X_even_temp,post_binned(coupling_bins_test),opts);
        clear X_even_temp coupling_bins_test

        % Store test likelihood for given this iteration 
        test_nll{ iter }(kk) = f_test;

    end
    
    iter = iter+1;
end
clear splits_indicator

%% STEP 5 FIT GLM

fprintf('Fit full model...\n');
[~,arg] = min(cellfun(@mean, test_nll));
timescale = timescales(arg);

hyper_params_tmp.XTimeSpan = timescale;
% Obtain basis
X_couple = calculate_X({pre},hyper_params_tmp);
if size(X_couple,1) > sum(couple_indicator)
    X_couple(sum(couple_indicator)+1:end,:) = [];
end

% Generate regressor matrix for all data
coupling_bins = logical(couple_indicator);
opts.rate = X_cPost(coupling_bins);
nparams = size(X_couple,2);

% Handle for the cost function
costfun = @(x) nll_poissPlasticity(x,X_couple,post_binned(coupling_bins),opts);

% Update some of the options
options.MaxIterations = fit_iter;
options.Display = 'iter-detailed';

% Initialize
bta_init = rand(nparams, 1);
% Optimize
[synParamsT,~,~,~] = fmincon(costfun, bta_init, [],[],[],[],...
        zeros(1,nparams),... 
        inf(1,nparams),...
        [],options);

% Get dt-resolution coupling term
cc = couple_indicator; cc(cc == 1) = X_couple*synParamsT;

% Compute excess synchrony - timeseries measure of spike transmission in units of Hz
lambda_full = exp(X_cPost + cc);
lambda_baseline = exp(X_cPost);
lambda_diff = lambda_full - lambda_baseline;

% Also store lambda as P(post = 1 | parameters)
lambda = lambda_full .* dt;

% Prepare kernel
gauss_window = 180;      % 3 minute size window (in seconds)
gauss_SD = 120;          % 2 minute standard deviation
grid = -gauss_window:dt:gauss_window;
gk = Gauss(grid, 0, gauss_SD); gk = gk(:);

gk = gausskernel(180/dt, 120/dt);
gk = gk ./ dt;

% Convolve excess synchrony and presynaptic spike train.
ldiff_conv = conv2fft(lambda_diff,gk,'same');
pre_binned_conv = conv2fft(pre_binned,gk,'same');
% Normalize to get excess synchrony per presynaptic spike
ltc_t = ldiff_conv ./ pre_binned_conv;
% Downsample
ltc_t = ltc_t(1:round(down_dt/dt):end);

%% STEP 6 SAVE

% We store the bare minimum in order to be able to post-hoc reconstruct all
% of the variables that were used to fit the model. 
ts_lastelement = ts(end);
pre_binned_inds = find(pre_binned);
post_binned_inds = find(post_binned);

fitting_duration = toc;

save(outfile, 'basepath', 'preID', 'postID', 'pre', 'post', 'synParamsT', 'test_nll', 'timescales', 'down_dt', 'ts_lastelement', ...
 'pre_binned_inds', 'post_binned_inds', 'ltc_t', 'lambda','argmax', 'delta', 'binsize', 'dt', 'fitting_duration', '-v7.3')
 

end
