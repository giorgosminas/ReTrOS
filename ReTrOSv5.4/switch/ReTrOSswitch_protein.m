function [profile, samples, posteriors, params] = ReTrOSswitch_protein(geneLocus,rawData,timescale,replicates,iterations,varargin)
%analyses a single timecourse dataset (with replicates), generating a 'transcriptional switch profile' and estimates of degradation rate
%do stuff with varargin - we have a number of optional parameters:
% 'showplot', <logical scalar> - if we want to view the raw output after analysis (default: true)
% 'figHandle', <graphics handle scalar> - handle to the figure to use for plotting
% 'plotSamples', <logical scalar> - if we want to plot all sample data, or just the posteriors (default: true) NOTE: only used if 'showplot' is true
% 'plotSwitchFitting', <logical scalar> - if we want to plot various stages of the getSwitchFit function (default: false) NOTE: only used if 'showplot' is true
% 'plotDegradationDist', <logical scalar> - if we want to plot estimated degradation rate distribution against the prior distribution (default: false) NOTE: only used if 'showplot' is true
% 'minSwitchStrength', <real scalar> - sets the minimum switch strength in the returned switch distribution (default: 0.2)
% 'switchEvidence', <real scalar> - sets the Bayes factor evidence threshold suggesting an N-switch model compared to a 0-switch model (default: 0)
% 'forceOnOffSwitches', <logical scalar> - if we want to enforce that the difference in birth rates before/after a switch point changes sign with each subsequent switch (default: false)
% 'minimumRateChange', <real scalar> - sets the level at which the minimum difference between two adjacent birth rates is allowed (impacts on overfitting) (default: 0)
% 'mRNADegradationRatePrior', <real scalar> - sets the alpha and beta parameters of the Gamma prior distribution for the mRNA degradation rate (default: [1.4828 0.1823])
% 'proteinDegradationRatePrior', <real scalar> - sets the alpha and beta parameters of the Gamma prior distribution for the protein degradation rate (default: [1.4828 0.1823])
% 'dataDetrend', <string> - sets the data detrending method to use (default: 'none')
% 'regressionMethod', <string> - sets the least squares regression method to use (default: 'leastsquares')
% 'useInitialSwitchSuggest', <logical scalar> - if we want to use an initial guess of switch points (default: true)
% 'initialSwitchThreshold', <real scalar> - sets the threshold in gradient space for suggesting initial switch points (default: 0.1)
% 'maxSamplesToUse', <integer scalar> - sets the maximum number of MCMC samples to estimate posterior distributions (default: # of iterations)
% 'maxSamplesToPlot', <integer scalar> - sets the maximum number of accepted MCMC samples plotted in figures (default: 10000)
% 'maxSwitches', <integer scalar> - sets the maximum number of switch points allowed by the algorithm (default: # of timepoints / 3)
% 'switchNumberPrior', <integer scalar> - sets the mean number of switch points to expect (default: 1)
% 'switchDelay', <real scalar> - sets the minimum delay between switch points (default: minimum time between timepoints)
% 'switchBaselineFactor', <real scalar> - sets the scaling factor of the switch density baseline modifier (default: 1)
% 'maxSwitchDeviationFactor' <real scalar> - sets the maximum switch deviation as a factor of time point sampling (default: 1)
% 'confidence', <integer scalar> - sets the confidence (number of standard deviations) at which to generate switch distributions (only used in plots) (default: 3)
% 'burnIn', <real scalar> - sets the proportion of samples from chain to remove as the 'burn in' period from beginning of chain (default: 0.25)
% 'percentiles', <real vector - sets the percentiles to return sigma and degradation rate estimates (default: [0.5 1 2.5 5 95 97.5 99 9.5])
% 'timepoints', <integer vector> - DEPRECATED sets which timepoints to use from the specified experimental dataset (default: all)
% 'replicates', <integer vecotr> - DEPRECATED sets which replicates to use from the specified experimental dataset (default: all)
% 'averageReplicates', <logical scalar> - DEPRECATED determines if we consider all replicates at each timepoint simulatenously (recommended), or if we average replicates generating a single value at each timepoint (default: false)
% 'maskSamples', <integer matrix> - allows the user to specify non-rectangular input data samples (replaces and supercedes 'timepoints' and 'replicates' parameters)

% disp('Required changes:');
% disp('Check that treatment name is valid for experiment - DONE');
% disp('Update "getSwitchFit" function to take max switches and force on/off parameter into account');
% disp('Switch posterior to take co-occurrence into account');

%we want to use all treatments in this experiment together
%they must have same number of timepoints and replicates
trueTimescale = timescale;
timepoints = length(timescale);


%parameter defaults - can be changed with function arguments
isPlot = true;
figHandle = [];
plotSamples = true;
plotSwitchFitting = false;
plotDegDist = false;
forceOnOff = false;
brRatioTol = 0;
maxSwitches = -1;%floor( (timepoints-2) / 2);
burnIn = 0.25;
confidence = 1.96;
confidences = [1 1.96 3];
percentiles = [0.5 1 2.5 5 95 97.5 99 99.5];
strengthThreshold = 0.2;
dimChangeThreshold = 0;
fullPlot = true;
estimationPoints = 1000;
useLogExplore = true;
userLambda = false;
switchDelay = -1;
maskedSamples = [];
switchBaselineFactor = 1;
maxSwitchDeviationFactor = 1;
explorePrec = true;
exploreAlpha = false;
suggestInitialSwitches = true;
initialSwTimes = [];
initialSwThreshold = 0.1;
detrend = 'none';
maxSamplesToPlot = 10000;
maxSamplesToUse = iterations;
regressionWeightPhiMax = 1;

regressionTypes = {'leastsquares','weightedleastsquares'};
regressionMethod = 'leastsquares';

% flatish prior
% priorAlpha_deg = 1;
% priorBeta_deg = 10;

%prior from literature estimates
% priorAlpha_deg = 1.5;
% priorBeta_deg = 0.2;

% fast prior
% priorAlpha_deg = 5;
% priorBeta_deg = 0.2;

%prior_deg = [10^0 10^3]+eps;
prior_deg_m = [1.4828 0.1823]; %Plant Cell. 2007 Nov;19(11):3418-36%[24.9519 0.0915];%luc mRNA deg measured in ROBUST project
prior_deg_p = [3.3765 0.0822]; %Mol Cell Proteomics. 2012 June; 11(6): M111.010025.%[155.0504 0.000835];%luc protein deg measured in ROBUST project
%[1 0.75];
prior_alpha = [1 1];

%%% Parse and check input parameters
nvarargin = length(varargin);
currArg = 1;
while currArg <= nvarargin
    if strcmpi(varargin{currArg},'showPlot')
        %we expect a single following argument, a logical, specifying whether we display plots of raw results after running
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"showPlot" requires a logical (true/false) scalar');
        else
            if islogical(varargin{currArg+1}) && isscalar(varargin{currArg+1})
                isPlot = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"showPlot" requires a logical (true/false) scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'plotSamples')
        %we expect a single following argument, a logical, specifying whether we display sample data, or just posteriors
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"plotSamples" requires a logical (true/false) scalar');
        else
            if islogical(varargin{currArg+1}) && isscalar(varargin{currArg+1})
                plotSamples = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"plotSamples" requires a logical (true/false) scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'plotSwitchFitting')
        %we expect a single following argument, a logical, specifying whether we plot various stages of the getSwitchFit function
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"plotSwitchFitting" requires a logical (true/false) scalar');
        else
            if islogical(varargin{currArg+1}) && isscalar(varargin{currArg+1})
                plotSwitchFitting = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"plotSwitchFitting" requires a logical (true/false) scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'plotDegradationDist')
        %we expect a single following argument, a logical, specifying whether we plot the estimated degradation rate distribution against the prior
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"plotDegradationDist" requires a logical (true/false) scalar');
        else
            if islogical(varargin{currArg+1}) && isscalar(varargin{currArg+1})
                plotDegDist = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"plotDegradationDist" requires a logical (true/false) scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'figHandle')
        %we expect a single argument, a non-negative integer, specifying the figure handle to use
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"figHandle" requies a graphics handle scalar');
        else
            if isscalar(varargin{currArg+1}) && ishghandle(varargin{currArg+1})
                figHandle = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"figHandle" requies a graphics handle scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'mRNADegradationRatePrior')
        %we expect a single following argument, either a string, or a 2-element non-negative float vector (mean and s.dev)
        if (currArg + 1) > nvarargin
            error('switchTool:varargcheck','"mRNADegradationRatePrior" requires a 2-element non-negative real vector');
        else
            if isreal(varargin{currArg+1})
                priorVal = varargin{currArg+1};
                if length(priorVal) == 2
                    if priorVal > 0
                        %                         alpha = priorVal(1) ^ 2 / priorVal(2) ^ 2;
                        %                         beta = priorVal(2) ^ 2 / priorVal(1);
                        %                         prior_deg_m = [alpha beta]+eps;
                        prior_deg_m = priorVal+eps;
                        currArg = currArg + 2;
                    else
                        error('switchTool:varargCheck','"mRNADegradationRatePrior" requires a 2-element non-negative real vector');
                    end
                else
                    error('switchTool:varargCheck','"mRNADegradationRatePrior" requires a 2-element non-negative real vector');
                end
            else
                error('switchTool:varargcheck','"mRNADegradationRatePrior" requires a 2-element non-negative real vector');
            end
        end
    elseif strcmpi(varargin{currArg},'proteinDegradationRatePrior')
        %we expect a single following argument, either a string, or a 2-element non-negative float vector (mean and s.dev)
        if (currArg + 1) > nvarargin
            error('switchTool:varargcheck','"proteinDegradationRatePrior" requires a 2-element non-negative real vector');
        else
            if isreal(varargin{currArg+1})
                priorVal = varargin{currArg+1};
                if length(priorVal) == 2
                    if priorVal > 0
                        %                         alpha = priorVal(1) ^ 2 / priorVal(2) ^ 2;
                        %                         beta = priorVal(2) ^ 2 / priorVal(1);
                        %                         prior_deg_p = [alpha beta]+eps;
                        prior_deg_p = priorVal+eps;
                        currArg = currArg + 2;
                    else
                        error('switchTool:varargCheck','"proteinDegradationRatePrior" requires a 2-element non-negative real vector');
                    end
                else
                    error('switchTool:varargCheck','"proteinDegradationRatePrior" requires a 2-element non-negative real vector');
                end
            else
                error('switchTool:varargcheck','"proteinDegradationRatePrior" requires a 2-element non-negative real vector');
            end
        end
    elseif strcmpi(varargin{currArg},'forceOnOffSwitches')
        %we expect a single following argument, a logical, specifying whether we force on-off structure of switches
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"forceOnOffSwitches" requires a logical (true/false) scalar');
        else
            if islogical(varargin{currArg+1}) && isscalar(varargin{currArg+1})
                forceOnOff = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"forceOnOffSwitches" requires a logical (true/false) scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'logExplore')
        %we expect a single following argument, a logical, specifying whether we explore log or linear parameter space
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"logExplore" requires a logical (true/false) scalar');
        else
            if islogical(varargin{currArg+1}) && isscalar(varargin{currArg+1})
                useLogExplore = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"logExplore" requires a logical (true/false) scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'samplePrecision')
        %we expect a single following argument, a logical, specifying whether we sample the precision, or fix an estimate
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"samplePrecision" requires a logical (true/false) scalar');
        else
            if islogical(varargin{currArg+1}) && isscalar(varargin{currArg+1})
                explorePrec = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"samplePrecision" requires a logical (true/false) scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'sampleAlpha')
        %we expect a single following argument, a logical, specifying whether we sample the alpha, or fix an estimate
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"sampleAlpha" requires a logical (true/false) scalar');
        else
            if islogical(varargin{currArg+1}) && isscalar(varargin{currArg+1})
                exploreAlpha = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"samplePrecision" requires a logical (true/false) scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'minimumRateChange')
        %we expect a single argument, a non-negative floating point value, specifying when we can accept changes in birth rates (help prevent over-fitting)
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"minimumRateChange" requires a non-negative real scalar');
        else
            if isreal(varargin{currArg+1}) && isscalar(varargin{currArg+1}) && varargin{currArg+1} >= 0
                brRatioTol = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"minimumRateChange" requires a non-negative real scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'minSwitchStrength')
        %we expect a single argument, a non-negative floating point value, specifying when we can accept switches based on strength
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"minSwitchStrength" requires a real scalar between 0 and 1');
        else
            if isreal(varargin{currArg+1}) && isscalar(varargin{currArg+1}) && varargin{currArg+1} >= 0 && varargin{currArg+1} <= 1
                strengthThreshold = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"minSwitchStrength" requires a real scalar between 0 and 1');
            end
        end
    elseif strcmpi(varargin{currArg},'switchEvidence')
        %we expect a single argument, a non-negative floating point value, specifying the Bayes factor evidence threshold suggesting an >0-switch model compared to a 0-switch model
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"switchEvidence" requires a non-negative real scalar');
        else
            if isreal(varargin{currArg+1}) && isscalar(varargin{currArg+1}) && varargin{currArg+1} >= 0
                dimChangeThreshold = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"switchEvidence" requires a non-negative real scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'maxSwitches')
        %we expect a single argument, a non-negative integer, specifying maximum number of switches (dimensions) the model can use
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"maxSwitches" requies a non-negative integer scalar');
        else
            if isint(varargin{currArg+1}) && isscalar(varargin{currArg+1}) && varargin{currArg+1} >= 0
                maxSwitches = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"maxSwitches" requies a non-negative integer scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'switchNumberPrior')
        %we expect a single argument, a non-negative value, specifying mean number of switches (dimensions) the model expects (uses a Poisson distribution)
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"switchNumberPrior" requies a non-negative real scalar');
        else
            if isreal(varargin{currArg+1}) && isscalar(varargin{currArg+1}) && varargin{currArg+1} >= 0
                lambda = varargin{currArg+1};
                userLambda = true;
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"switchNumberPrior" requies a non-negative real scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'switchDelay')
        %we expect a single argument, a non-negative value, specifying the minimum time delay between switches (in hours)
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"switchDelay" requies a non-negative real scalar');
        else
            if isreal(varargin{currArg+1}) && isscalar(varargin{currArg+1}) && varargin{currArg+1} >= 0
                switchDelay = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"switchDelay" requies a non-negative real scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'switchBaselineFactor')
        %we expect a single argument, a real scalar, specifying the scaling factor of the switch density baseline modifier
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"switchBaselineFactor" requies a real scalar');
        else
            if isreal(varargin{currArg+1}) && isscalar(varargin{currArg+1}) && varargin{currArg+1} >= 0
                switchBaselineFactor = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"switchBaselineFactor" requies a real scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'maxSwitchDeviationFactor')
        %we expect a single argument, a non-negative floating point value, specifying the scaling factor of time point sampling specifying maximum switch deviation
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"maxSwitchDeviationFactor" requires a non-negative real scalar');
        else
            if isreal(varargin{currArg+1}) && isscalar(varargin{currArg+1}) && varargin{currArg+1} >= 0
                maxSwitchDeviationFactor = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"maxSwitchDeviationFactor" requires a non-negative real scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'confidence')
        %we expect a single argument, a non-negative integer, specifying confidence (# SDs) of switch distributions
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"confidence" requies an integer scalar between 1 and 3');
        else
            if isint(varargin{currArg+1}) && isscalar(varargin{currArg+1}) && varargin{currArg+1} >= 1 && varargin{currArg+1} <= 3
                confidence = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"confidence" requies an integer scalar between 1 and 3');
            end
        end
    elseif strcmpi(varargin{currArg},'burnIn')
        %we expect a single argument, a non-negative, max 1, floating point value, specifying how much of the beginning of the chain to remove as burn in
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"burnIn" requires a real scalar between 0 and 1');
        else
            if isreal(varargin{currArg+1}) && isscalar(varargin{currArg+1}) && varargin{currArg+1} >= 0 && varargin{currArg+1} <= 1
                burnIn = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"burnIn" requires a real scalar between 0 and 1');
            end
        end
    elseif strcmpi(varargin{currArg},'percentiles')
        %we expect a single argument, a non-negative floating point vector, specifying lower and upper percentiles for sigma and degradation rate estimates
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"percentiles" requires a non-negative, 2 element real vector');
        else
            if isreal(varargin{currArg + 1}) && length( varargin{currArg+1} ) == 2 && sum(varargin{currArg +1} >= 0) == 2 && sum(varargin{currArg +1} <= 100) == 2
                percentiles = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"percentiles" requires a non-negative, 2 element real vector');
            end
        end
    elseif strcmpi(varargin{currArg},'maskSamples')
        %we expect a single argument, a non-negative integer vector, specifying samples to ignore.
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"maskSamples" requires a non-negative, non-zero integer vector');
        else
            if isint(varargin{currArg+1}) && (length(varargin{currArg+1}) < (timepoints*replicates)) && (sum(varargin{currArg +1} >= 1) == length(varargin{currArg+1})) && (sum(varargin{currArg +1} <= (timepoints*replicates)) == length(varargin{currArg+1}))
                maskedSamples = unique(sort(varargin{currArg+1},'ascend'));
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"maskSamples" requires a non-negative, non-zero integer vector');
            end
        end
    elseif strcmpi(varargin{currArg},'initialSwitchThreshold')
        %we expect a single argument, a non-negative, max 1, floating point value, specifying threshold in gradient space for initial switch selection
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"initialSwitchThreshold" requires a real scalar between 0 and 1');
        else
            if isreal(varargin{currArg+1}) && isscalar(varargin{currArg+1}) && varargin{currArg+1} >= 0 && varargin{currArg+1} <= 1
                initialSwThreshold = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"initialSwitchThreshold" requires a real scalar between 0 and 1');
            end
        end
    elseif strcmpi(varargin{currArg},'detrendData')
        %we expect a single argument, a string
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"detrendData" requires a string');
        else
            if ischar(varargin{currArg+1})
                detrend = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"detrendData" requires a string');
            end
        end
    elseif strcmpi(varargin{currArg},'regressionMethod')
        %we expect a single argument, a string
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"regressionMethod" requires a string');
        else
            if ischar(varargin{currArg+1})
                regressionMethod = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"regressionMethod" requires a string');
            end
        end
    elseif strcmpi(varargin{currArg},'initialSwitches')
        %we expect a single argument, a non-negative integer vector, specifying initial switch points.
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"initialSwitches" requires a non-negative, non-zero integer vector');
        else
            if isnumeric(varargin{currArg+1}) && (length(varargin{currArg+1}) < timepoints) && (min(varargin{currArg +1}) > timescale(1)) && (max(varargin{currArg +1}) < timescale(end))
                initialSwTimes = unique(varargin{currArg+1});
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"initialSwitches" requires a non-negative, non-zero integer vector');
            end
        end
    elseif strcmpi(varargin{currArg},'useInitialSwitchSuggest')
        %we expect a single following argument, a logical, specifying whether we want to start from a 0-switch model, or from a suggested set of switches generated from data
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"useInitialSwitchSuggest" requires a logical (true/false) scalar');
        else
            if islogical(varargin{currArg+1}) && isscalar(varargin{currArg+1})
                suggestInitialSwitches = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"useInitialSwitchSuggest" requires a logical (true/false) scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'maxSamplesToUse')
        %we expect a single argument, a non-negative integer, specifying maximum number of samples to use for posteriors (thinning)
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"maxSamplesToUse" requies a non-negative integer scalar');
        else
            if isint(varargin{currArg+1}) && isscalar(varargin{currArg+1}) && varargin{currArg+1} >= 0
                maxSamplesToUse = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"maxSamplesToUse" requies a non-negative integer scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'maxSamplesToPlot')
        %we expect a single argument, a non-negative integer, specifying maximum number of samples to plot in figure
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"maxSamplesToPlot" requies a non-negative integer scalar');
        else
            if isint(varargin{currArg+1}) && isscalar(varargin{currArg+1}) && varargin{currArg+1} >= 0
                maxSamplesToPlot = varargin{currArg+1};
                currArg = currArg + 2;
            else
                error('switchTool:varargCheck','"maxSamplesToPlot" requies a non-negative integer scalar');
            end
        end
    elseif strcmpi(varargin{currArg},'timepoints')
        %we expect a single argument, a non-negative integer vector, specifying which timepoints to use
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"timepoints" requires a non-negative integer vector NOTE: parameter deprecated');
        else
            currArg = currArg + 2;
            disp('NOTE: "timepoints" parameter deprecated, ignoring');
        end
    elseif strcmpi(varargin{currArg},'replicates')
        %we expect a single argument, a non-negative integer vector, specifying which replicates to use
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"replicates" requires a non-negative integer vector NOTE: parameter deprecated');
        else
            currArg = currArg + 2;
            disp('NOTE: "replicates" parameter deprecated, ignoring');
        end
    elseif strcmpi(varargin{currArg},'averageReplicates')
        %we expect a single argument, a logical, specifying if we want to average replicates (not recommended)
        if (currArg + 1) > nvarargin
            error('switchTool:varargCheck','"averageReplicates" requires a logical (true/false) scalar NOTE: parameter deprecated');
        else
            currArg = currArg + 2;
            disp('NOTE: "averageReplicates" parameter deprecated, ignoring');
        end
    else
        error('switchTool:varargCheck', ['Invalid parameter (' int2str(currArg) ') - "' varargin{currArg} '"']);
    end
end

if ~isPlot
    plotSwitchFitting = false;
    plotSamples = false;
end

%if our timescales dont start at 0 the algorithm fails
%therefore we must shift the timescale to accomodate this, and inform the user
% shiftedTimescale = false;
timescaleShift = 0;
if timescale(1) > 0
    %     shiftedTimescale = true;
    timescaleShift = timescale(1);
    timescale = timescale - timescaleShift;
end
%timescale

if switchDelay > 0
    if switchDelay >= (timescale(end)/2)
        error('switchTool:switchDelay','"switchDelay" is too large for timescale');
    end
else
    deltaTime = zeros(1,timepoints-1);
    for x=1:length(deltaTime)
        deltaTime(x) = timescale(x+1) - timescale(x);
    end
    switchDelay = min(deltaTime)*0.75;
end

maxSwitchesDelay = floor((timescale(end)- 2*switchDelay) / switchDelay);
maxSwitchesDelay = min(maxSwitchesDelay,timepoints-2);
if maxSwitchesDelay < maxSwitches
    disp(['Current maximum switch number (' int2str(maxSwitches) ') too large for switch delay']);
    disp(['Using new maximum switch number of ' int2str(maxSwitchesDelay)]);
    maxSwitches = maxSwitchesDelay;
end
if maxSwitches == -1
    maxSwitches = maxSwitchesDelay;
elseif maxSwitches == 0
    maxSwitches = 1;
end
if ~userLambda
    lambda = maxSwitches/2;
end
if lambda > maxSwitches
    disp(['Switch number lambda (' int2str(lambda) ') larger than max switches, setting to max switches (' int2str(maxSwitches) ')']);
    lambda = maxSwitches;
elseif lambda == 0
    lambda = 1;
end

maxswTimes = timescale(end);%timescale(end) - switchDelay;
minswTimes = 0;%switchDelay;
maxDegRate_mRNA = 10^1;    %min half-life: log(2) / 10  < 5 minutes
minDegRate_mRNA = 10^-2; %max half-life: log(2) / 0.01 > 24 hours
maxDegRate_protein = 10^1;    %min protein: log(2) / 10  < 5 minutes
minDegRate_protein = 10^-2; %max protein: log(2) / 0.01 > 24 hours
maxAlpha = 10^2;
minAlpha = 10^-2;


%%% set up priors
%lambda = 4;     %suggested expected number of switches
priorSwitchNum = calculateSwitchNumPrior(maxSwitches,lambda);

changeDimConst = 0.4;


data = rawData;
useDatapoint = data >= 0;
useDatapoint(maskedSamples) = 0;
data(useDatapoint==0) = NaN;

d = reshape(data,timepoints,replicates);
%we have several ways to extract data, depending on input arguments

% we need to do some additional work to check for missing values

Y = reshape(d,timepoints * replicates,1);
data = d';

switch detrend
    %     case 'log-linear'
    %         [rawData trend] = detrendLinear(log(rawData),timescale,replicates);
    %         rawData = exp(rawData)';
    %         trend = exp(trend)';
    
    case 'log-linear'
        disp('Applying log-linear data detrend');
        newData = nan(replicates,timepoints);
        trendData = nan(replicates,timepoints);
        betas = nan(replicates,2);
        for z = 1:replicates
            [detrendData, trend, beta] = detrendLinear(log(data( z,: )),timescale-timescale(1),1);
            newData(z,:) = exp(detrendData)';
            trendData(z,:) = exp(trend)';
            %newData(((z-1)*timepoints+1):(z*timepoints)) = newData(((z-1)*timepoints+1):(z*timepoints)) +  trendData((z-1)*timepoints+1);
            betas(z,:) = beta;
        end
        data = newData;
        Y = reshape(data',timepoints * replicates,1);
    case 'linear'
        disp('Applying linear data detrend');
        newData = nan(replicates,timepoints);
        trendData = nan(replicates,timepoints);
        betas = nan(replicates,2);
        for z = 1:replicates
            [detrendData, trend, beta] = detrendLinear(data(z,:),timescale-timescale(1),1);
            newData(z,:) = detrendData';
            trendData(z,:) = trend';
            %newData(((z-1)*timepoints+1):(z*timepoints)) = newData(((z-1)*timepoints+1):(z*timepoints)) +  trendData((z-1)*timepoints+1);
            betas(z,:) = beta;
        end
        data = newData;
        Y = reshape(data',timepoints * replicates,1);
end


%constants used to calculate precision
TC = sum(Y>=0);
prec_prior_r = 10^1;%0.001;
prec_prior_lambda = 10^1;%0.001;
p_1 = prec_prior_r + 0.5 * TC;

maxData = max(max(data));

minData = min(min(data));


maxP_0 = maxData*2; %should estimate these from data!
minP_0 = max(0,minData*0.5);  %should estimate these from data!


degRate_mRNA = fastgamrnd(prior_deg_m(1),prior_deg_m(2),1);%fastnormrnd(0.75,0.25,1);   %degradation rate of mRNA
degRate_mRNA = max(min(degRate_mRNA,maxDegRate_mRNA),minDegRate_mRNA);

degRate_protein = fastgamrnd(prior_deg_p(1),prior_deg_p(2),1);%fastnormrnd(0.75,0.25,1);   %degradation rate of mRNA
degRate_protein = max(min(degRate_protein,maxDegRate_protein),minDegRate_protein);

if exploreAlpha
    alpha = fastgamrnd(prior_alpha(1),prior_alpha(2),1);%fastnormrnd(1,0.25,1);
    alpha = max(min(alpha,maxAlpha),minAlpha);
else
    alpha = 1;
end

prec = 1;

eHatMat = zeros(replicates,timepoints);
e_hat = ones(replicates*timepoints,1);
Xmatrix = cell(maxSwitches+1,1);
for x = 1:(maxSwitches+1)
    Xmatrix{x} = zeros(timepoints * replicates,x+2);
end

lLikelihoodHistory = zeros(iterations,1);
degRate_mRNAHistory = zeros(iterations,1);
degRate_proteinHistory = zeros(iterations,1);
alphaHistory = zeros(iterations,1);
%P_0History = zeros(iterations,1);
%M_0History = zeros(iterations,1);
swTimesHistory = cell(iterations,1);
swTimesHistory{iterations} = [];
%tausHistory = cell(iterations,1);
%tausHistory{iterations} = [];
sigma2History = zeros(iterations,1);
%mRNAFitHistory = zeros(iterations,timepoints);
%proteinFitHistory = zeros(iterations,timepoints);

swTimes = initialSwTimes - timescaleShift;
if suggestInitialSwitches
    if isempty(initialSwTimes) %user hasnt supplied switch times and wants suggestions from data
        %         %need to update if we have replicates
        %         if replicates > 1
        %             d = median(data);
        %         else
        %             d = data;
        %         end
        %
        %         %simply check for changes in sign of first derivative
        %         g1 = gradient(d)
        %         d = (g1>0) + (g1<0).*-1
        %         swTimes = timescale(diff(d) ~= 0)
        %lets be a little smarter
        f = feval(fit(repmat(timescale',replicates,1),Y,'smoothingspline'),timescale);
        g1 = gradient(f)./max(max(gradient(f)),abs(min(gradient(f))));
        
        maskPos = ((g1 > 0) .* (g1 <= initialSwThreshold)) == 1;
        maskNeg = ((g1 < 0) .* (g1 >= -initialSwThreshold)) == 1;
        
        g1(maskPos) = 0;
        g1(maskNeg) = 0;
        
        d = (g1>0) + (g1<0).*-1;
        
        idx = find(d,1);
        if idx > 0
            curr = d(idx);
            swTimes = timescale(idx);
            for x = 1:length(d)
                if curr == -d(x)
                    swTimes = [swTimes timescale(x)];
                    curr = d(x);
                end
            end
            swTimes = swTimes(swTimes > timescale(1));
            swTimes = swTimes(swTimes < timescale(end));
        end
    else
        suggestInitialSwitches = false;
        disp('Using user supplied initial switch positions');
    end
end

if length(swTimes) > maxSwitches
    swTimes = [];
end

% set up regression methods
regression.type = regressionMethod;
if strcmpi(regression.type,'weightedleastsquares')
    regression.weightsPhi = [];
    regression.weights = [];
    regression.Q = [];
    regression.weights = feval(fit(repmat(timescale',replicates,1),Y,'smoothingspline'),timescale);
    regressionWeightsPhiHistory = nan(iterations,1) ;
end

temperature = 1;

dimChangeUp = nan(iterations,maxSwitches+1);
dimChangeDown = nan(iterations,maxSwitches+1);
%acceptRate = 0;
n = 1;
earlyConverged = false; %this refers to if we pass the burn-in period
burnInIterations = round(iterations * burnIn);
while n <= iterations && ~earlyConverged
    
    if strcmp(regression.type,'weightedleastsquares')
        yphi = log(prec*e_hat.^2) ;
        xphi = -repmat(log(regression.weights),replicates,1) ;
        phi_hat = ((xphi'*xphi)\xphi')*yphi;
        phi_hat = max(min(phi_hat,regressionWeightPhiMax),0);
        regression.weightsPhi = phi_hat;
        q = repmat(regression.weights.^(-regression.weightsPhi),replicates,1);
        regression.Q = diag(q);
    end
    
    %calculate baseline model
    [P_0 M_0 taus] = calculateModel_protein(timescale,swTimes,degRate_mRNA,degRate_protein,alpha,Y,timepoints,replicates,useDatapoint,Xmatrix,regression);
    [mRNA protein] = nStateSwitchODE_protein(timescale,M_0,taus,degRate_mRNA,swTimes,P_0,alpha,degRate_protein);
    %calculate SSE of profile to data
    e_hat = calculateEHat(protein,data,eHatMat,replicates,useDatapoint);
    
    if explorePrec
        if strcmpi(regressionMethod,'weightedleastsquares')
            sse = e_hat'*(regression.Q\e_hat);
        else
            sse = sum(e_hat .^ 2);
        end
        
        p_2 = prec_prior_lambda + 0.5 * sse;
        
        %calculate precision (used for sigma_sq) - Gibbs
        prec = (1/p_2) * randg(p_1);
    end
    lLikelihood = logLikelihood(e_hat,prec,regression) * temperature;
    
    
    %update delta_p - MH
    %degRateProposed = exp( log(degRate) + (rand(1)/2 - 0.25) );
    if useLogExplore
        degRateProposed = exp(log(degRate_protein) + fastnormrnd(0,0.25,1));
        while degRateProposed < minDegRate_protein || maxDegRate_protein < degRateProposed
            degRateProposed = exp(log(degRate_protein) + fastnormrnd(0,0.25,1));
        end
        proposalRatio = degRateProposed / degRate_protein;
    else
        degRateProposed = degRate_protein + fastnormrnd(0,0.25,1);
        while degRateProposed < minDegRate_protein || maxDegRate_protein < degRateProposed
            degRateProposed = degRate_protein + fastnormrnd(0,0.25,1);
        end
        proposalRatio = 1;
    end
    %degRateProposed = min(degRateProposed,maxDegRate);
    %degRateProposed = max(degRateProposed,minDegRate);
    
    %     if degRateProposed >= minDegRate && maxDegRate > degRateProposed
    %calculate model with proposed degradation rate
    [P_0Proposed M_0Proposed tausProposed] = calculateModel_protein(timescale,swTimes,degRate_mRNA,degRateProposed,alpha,Y,timepoints,replicates,useDatapoint,Xmatrix,regression);
    [mRNAProposed proteinProposed] = nStateSwitchODE_protein(timescale,M_0Proposed,tausProposed,degRate_mRNA,swTimes,P_0Proposed,alpha,degRateProposed);
    e_hat = calculateEHat(proteinProposed,data,eHatMat,replicates,useDatapoint);
    lLikelihoodProposed = logLikelihood(e_hat,prec,regression) * temperature;
    
    lLikelihoodRatio = exp(lLikelihoodProposed - lLikelihood);
    priorRatio = gampdf(degRateProposed,prior_deg_p(1),prior_deg_p(2)) /  gampdf(degRate_protein,prior_deg_p(1),prior_deg_p(2));
    
    %acceptance: likelihood ratio * prior proposal * proposal
    acceptanceDegRate = min(1, lLikelihoodRatio * priorRatio * proposalRatio);
    
    if sum(tausProposed < 0) > 0
        noFail = 0;
    else
        noFail = 1;
    end
    acceptanceDegRate = acceptanceDegRate * (M_0Proposed >= 0) * (P_0Proposed >= minP_0) * (P_0 <= maxP_0) * noFail;
    %     else
    %         acceptanceDegRate = 0;
    %     end
    
    if rand(1) <= acceptanceDegRate
        degRate_protein = degRateProposed;
        lLikelihood = lLikelihoodProposed;
        %acceptRate = acceptRate + 1;
    end
    
    %update delta_m - MH
    %degRateProposed = exp( log(degRate) + (rand(1)/2 - 0.25) );
    if useLogExplore
        degRateProposed = exp(log(degRate_mRNA) + fastnormrnd(0,0.25,1));
        while degRateProposed < minDegRate_mRNA || maxDegRate_mRNA < degRateProposed
            degRateProposed = exp(log(degRate_mRNA) + fastnormrnd(0,0.25,1));
        end
        proposalRatio = degRateProposed / degRate_mRNA;
    else
        degRateProposed = degRate_mRNA + fastnormrnd(0,0.25,1);
        while degRateProposed < minDegRate_mRNA || maxDegRate_mRNA < degRateProposed
            degRateProposed = degRate_mRNA + fastnormrnd(0,0.25,1);
        end
        proposalRatio = 1;
    end
    %degRateProposed = min(degRateProposed,maxDegRate);
    %degRateProposed = max(degRateProposed,minDegRate);
    
    %     if degRateProposed >= minDegRate && maxDegRate > degRateProposed
    %calculate model with proposed degradation rate
    [P_0Proposed M_0Proposed tausProposed] = calculateModel_protein(timescale,swTimes,degRateProposed,degRate_protein,alpha,Y,timepoints,replicates,useDatapoint,Xmatrix,regression);
    [mRNAProposed proteinProposed] = nStateSwitchODE_protein(timescale,M_0Proposed,tausProposed,degRateProposed,swTimes,P_0Proposed,alpha,degRate_protein);
    e_hat = calculateEHat(proteinProposed,data,eHatMat,replicates,useDatapoint);
    lLikelihoodProposed = logLikelihood(e_hat,prec,regression) * temperature;
    
    lLikelihoodRatio = exp(lLikelihoodProposed - lLikelihood);
    priorRatio = gampdf(degRateProposed,prior_deg_m(1),prior_deg_m(2)) /  gampdf(degRate_mRNA,prior_deg_m(1),prior_deg_m(2));
    
    %acceptance: likelihood ratio * prior proposal * proposal
    acceptanceDegRate = min(1, lLikelihoodRatio * priorRatio * proposalRatio);
    
    if sum(tausProposed < 0) > 0
        noFail = 0;
    else
        noFail = 1;
    end
    acceptanceDegRate = acceptanceDegRate * (M_0Proposed >= 0) * (P_0Proposed >= minP_0) * (P_0 <= maxP_0) * noFail;
    %     else
    %         acceptanceDegRate = 0;
    %     end
    
    if rand(1) <= acceptanceDegRate
        degRate_mRNA = degRateProposed;
        lLikelihood = lLikelihoodProposed;
        %acceptRate = acceptRate + 1;
    end
    
    %update alpha - MH
    if exploreAlpha
        if useLogExplore
            alphaProposed = exp(log(alpha) + fastnormrnd(0,0.5,1));
            while alphaProposed < minAlpha || maxAlpha < alphaProposed
                alphaProposed = exp(log(alpha) + fastnormrnd(0,0.5,1));
            end
            proposalRatio = alphaProposed / alpha;
        else
            alphaProposed = alpha + fastnormrnd(0,0.25,1);
            while alphaProposed < minAlpha || maxAlpha < alphaProposed
                alphaProposed = alpha + fastnormrnd(0,0.25,1);
            end
            proposalRatio = 1;
        end
        %degRateProposed = min(degRateProposed,maxDegRate);
        %degRateProposed = max(degRateProposed,minDegRate);
        
        %     if degRateProposed >= minDegRate && maxDegRate > degRateProposed
        %calculate model with proposed alpha
        [P_0Proposed M_0Proposed tausProposed] = calculateModel_protein(timescale,swTimes,degRate_mRNA,degRate_protein,alphaProposed,Y,timepoints,replicates,useDatapoint,Xmatrix,regression);
        [mRNAProposed proteinProposed] = nStateSwitchODE_protein(timescale,M_0Proposed,tausProposed,degRate_mRNA,swTimes,P_0Proposed,alphaProposed,degRate_protein);
        e_hat = calculateEHat(proteinProposed,data,eHatMat,replicates,useDatapoint);
        lLikelihoodProposed = logLikelihood(e_hat,prec,regression) * temperature;
        
        lLikelihoodRatio = exp(lLikelihoodProposed - lLikelihood);
        priorRatio = gampdf(alphaProposed,prior_alpha(1),prior_alpha(2)) /  gampdf(alpha,prior_alpha(1),prior_alpha(2));
        
        %acceptance: likelihood ratio * prior proposal * proposal
        acceptanceAlpha = min(1, lLikelihoodRatio * priorRatio * proposalRatio);
        
        if sum(tausProposed < 0) > 0
            noFail = 0;
        else
            noFail = 1;
        end
        acceptanceAlpha = acceptanceAlpha * (M_0Proposed >= 0) * (P_0Proposed >= minP_0) * (P_0 <= maxP_0) * noFail;
        %     else
        %         acceptanceDegRate = 0;
        %     end
        
        if rand(1) <= acceptanceAlpha
            alpha = alphaProposed;
            lLikelihood = lLikelihoodProposed;
            %acceptRate = acceptRate + 1;
        end
    end
    
    %update switch position/number
    switchNum = length(swTimes);
    
    %we only allow single dimension jumps, so only have three moves
    if switchNum == 0
        pBirthOne = 0;
        if maxSwitches > 0
            pBirthOne = 1;
        end
        pDeathOne = 0;
        pMove = 0;
    else
        if switchNum < maxSwitches
            pBirthOne = changeDimConst * min(1,priorSwitchNum(switchNum + 2) / priorSwitchNum(switchNum + 1));
        else
            pBirthOne = 0;
        end
        if switchNum > 0
            pDeathOne = changeDimConst * min(1,priorSwitchNum(switchNum) / priorSwitchNum(switchNum +1));
        else
            pDeathOne = 0;
        end
        
        pMove = 1 - pBirthOne - pDeathOne;
    end
    
    r = rand(1);
    acceptanceSwitch = 0;
    if r <= pMove
        %Move a single switch point
        
        swTimesProposed = swTimes;
        %pick a switch uniformly
        switchToMove = randi(switchNum);
        %we move the switch (again uniformly) between the previous and next switch
        if switchToMove == 1
            minRange = minswTimes;
            if switchNum == 1
                maxRange = maxswTimes;
            else
                maxRange = swTimes(2) - switchDelay;
            end
        else
            minRange = swTimes(switchToMove-1) + switchDelay;
            if switchToMove == switchNum
                maxRange = maxswTimes;
            else
                maxRange = swTimes(switchToMove + 1) - switchDelay;
            end
        end
        range = maxRange - minRange;
        if range > 0
            move = range * rand(1);
            swTimesProposed(switchToMove) = minRange + move;
            
            %calculate model with proposed swTimes
            [P_0Proposed M_0Proposed tausProposed] = calculateModel_protein(timescale,swTimesProposed,degRate_mRNA,degRate_protein,alpha,Y,timepoints,replicates,useDatapoint,Xmatrix,regression);
            [mRNAProposed proteinProposed] = nStateSwitchODE_protein(timescale,M_0Proposed,tausProposed,degRate_mRNA,swTimesProposed,P_0Proposed,alpha,degRate_protein);
            e_hat = calculateEHat(proteinProposed,data,eHatMat,replicates,useDatapoint);
            lLikelihoodProposed = logLikelihood(e_hat,prec,regression) * temperature;
            
            lLikelihoodRatio = exp(lLikelihoodProposed - lLikelihood);
            
            distanceToSwitch = [maxRange - swTimesProposed(switchToMove), swTimesProposed(switchToMove) - minRange, maxRange - swTimes(switchToMove), swTimes(switchToMove) - minRange];
            acceptanceSwitch = min(1,lLikelihoodRatio * ( (distanceToSwitch(1) * distanceToSwitch(2)) / (distanceToSwitch(3) * distanceToSwitch(4)) ) );
        end
    elseif r <= (pMove + pBirthOne)
        %Add a single switch point
        
        %new switch proposed uniformly over entire timescale
        maxRange = maxswTimes;
        minRange = minswTimes;
        range = maxRange - minRange;
        newSwitch = (range * rand(1)) + minRange;
        %need to find interval in which this new switch falls
        if switchNum > 0
            foundSwitch = false;
            %check beginning and end separately
            if minRange <= newSwitch && newSwitch <= swTimes(1)
                insertSwitch = 1;
                maxRange = swTimes(1) - switchDelay;
            else
                currSwitch = 1;
                while ~foundSwitch && currSwitch < switchNum
                    if newSwitch < swTimes(currSwitch+1)
                        insertSwitch = currSwitch+1;
                        minRange = swTimes(currSwitch) + switchDelay;
                        maxRange = swTimes(currSwitch+1) - switchDelay;
                        foundSwitch = true;
                    end
                    currSwitch = currSwitch + 1;
                end
                %if its not found, then we know its at the end
                if ~foundSwitch
                    minRange = swTimes(switchNum) + switchDelay;
                    maxRange = maxswTimes;
                    insertSwitch = switchNum + 1;
                end
            end
        else
            insertSwitch = 1;
        end
        swTimesProposed = zeros(1,switchNum + 1);
        for x = 1:(insertSwitch-1)
            swTimesProposed(x) = swTimes(x);
        end
        swTimesProposed(insertSwitch) = newSwitch;
        for x = (insertSwitch+1):(switchNum+1)
            swTimesProposed(x) = swTimes(x-1);
        end
        range = maxRange - minRange;
        if range > 0
            %calculate model with proposed swTimes
            [P_0Proposed M_0Proposed tausProposed] = calculateModel_protein(timescale,swTimesProposed,degRate_mRNA,degRate_protein,alpha,Y,timepoints,replicates,useDatapoint,Xmatrix,regression);
            [mRNAProposed proteinProposed] = nStateSwitchODE_protein(timescale,M_0Proposed,tausProposed,degRate_mRNA,swTimesProposed,P_0Proposed,alpha,degRate_protein);
            e_hat = calculateEHat(proteinProposed,data,eHatMat,replicates,useDatapoint);
            lLikelihoodProposed = logLikelihood(e_hat,prec,regression) * temperature;
            
            lLikelihoodRatio = exp(lLikelihoodProposed - lLikelihood);
            
            %we know thw range is positive (including the switch delay on either side, so we just need to calculate how many switchDelay instances
            %our time space for switches is reduced by enforcing a switch delay between switches
            %invalidLength = (length(swTimes) + 1 + 2) * switchDelay;
            invalidLength = switchDelay * 2;
            for m = 1:length(swTimes)
                if m == 1
                    left = minswTimes;
                else
                    left = swTimes(m-1);
                end
                invalidLength = invalidLength + min(swTimes(m) - left,switchDelay);
                
                if m == length(swTimes)
                    invalidLength = invalidLength + min(maxswTimes-swTimes(m),switchDelay);
                end
            end
            
            %priorRatio = (priorSwitchNum(switchNum+2) / priorSwitchNum(switchNum +1)) * ((2 * (switchNum + 1) * (2 * switchNum + 3)) / ((maxswTimes-minswTimes) ^ 2)) * (((newSwitch - minRange) * (maxRange - newSwitch)) / (maxRange - minRange));
            %add in invalid length term
            priorRatio = (priorSwitchNum(switchNum+2) / priorSwitchNum(switchNum +1)) * (((2 * switchNum + 1) * (2 * switchNum + 3)) / ((timescale(end) - invalidLength)^2)) * (((newSwitch - minRange) * (maxRange - newSwitch)) / (maxRange - minRange));
            
            pDeath1 = changeDimConst * min(1,priorSwitchNum(switchNum + 1) / priorSwitchNum(switchNum + 2));
            %proposalRatio = (pDeath1 * (maxswTimes * minswTimes)) / (pBirthOne * (switchNum + 2)); %spotted by Kirsty, should be maxswTimes-minSwTimes
            %proposalRatio = (pDeath1 * (maxswTimes - minswTimes)) / (pBirthOne * (switchNum + 2));
            %add in invalid length term
            proposalRatio = (pDeath1 * ((timescale(end) - invalidLength))) / (pBirthOne * (switchNum + 2));
            acceptanceSwitch = min(1,lLikelihoodRatio * priorRatio * proposalRatio );
            
            %record log likelihood ratio
            dimChangeUp(n,switchNum+1) = lLikelihoodProposed - lLikelihood;
        else
            swTimesProposed = swTimes;
        end
    elseif r <= (pMove + pBirthOne + pDeathOne)
        %Remove a single switch point
        
        switchToRemove = randi(switchNum);
        if switchToRemove == 1
            minRange = minswTimes;
            if switchNum == 1
                maxRange = maxswTimes;
            else
                maxRange = swTimes(2) - switchDelay;
            end
        else
            minRange = swTimes(switchToRemove-1) + switchDelay;
            if switchNum == switchToRemove
                maxRange = maxswTimes;
            else
                maxRange = swTimes(switchToRemove+1) - switchDelay;
            end
        end
        swTimesProposed = [swTimes(1:(switchToRemove-1)) swTimes((switchToRemove+1):switchNum)];
        oldSwitch = swTimes(switchToRemove);
        %calculate model with proposed swTimes
        [P_0Proposed M_0Proposed tausProposed] = calculateModel_protein(timescale,swTimesProposed,degRate_mRNA,degRate_protein,alpha,Y,timepoints,replicates,useDatapoint,Xmatrix,regression);
        [mRNAProposed proteinProposed] = nStateSwitchODE_protein(timescale,M_0Proposed,tausProposed,degRate_mRNA,swTimesProposed,P_0Proposed,alpha,degRate_protein);
        e_hat = calculateEHat(proteinProposed,data,eHatMat,replicates,useDatapoint);
        lLikelihoodProposed = logLikelihood(e_hat,prec,regression) * temperature;
        lLikelihoodRatio = exp(lLikelihoodProposed - lLikelihood);
        
        invalidLength = switchDelay * 2;
        for m = 1:length(swTimes)
            if m == 1
                left = minswTimes;
            else
                left = swTimes(m-1);
            end
            invalidLength = invalidLength + min(swTimes(m) - left,switchDelay);
            
            if m == length(swTimes)
                invalidLength = invalidLength + min(maxswTimes-swTimes(m),switchDelay);
            end
        end
        
        %length to top, change big interval over smaller intervals
        %priorRatio = (priorSwitchNum(switchNum) / priorSwitchNum(switchNum + 1)) * (( (maxswTimes-minswTimes) ^ 2) / (2 * (switchNum + 1) * (2 * switchNum + 3) )) * ((maxRange - minRange) / ((oldSwitch - minRange) * (maxRange - oldSwitch)));
        priorRatio = (priorSwitchNum(switchNum) / priorSwitchNum(switchNum + 1)) * (( (timescale(end)- invalidLength)^2) / ((2 * switchNum + 1) * (2 * switchNum + 3) )) * ((maxRange - minRange) / ((oldSwitch - minRange) * (maxRange - oldSwitch)));
        pBirth1 = changeDimConst * min(1,priorSwitchNum(switchNum+1) / priorSwitchNum(switchNum));
        %birth on top + relabel
        %proposalRatio = ((pBirth1 * (switchNum) / pDeathOne * (maxswTimes-minswTimes)));
        proposalRatio = ((pBirth1 * (switchNum) / pDeathOne * (timescale(end) - invalidLength)));
        acceptanceSwitch = min(1,lLikelihoodRatio * priorRatio * proposalRatio );
        
        %record log likelihood ratio
        dimChangeDown(n,switchNum+1) = lLikelihoodProposed - lLikelihood;
    end
    
    
    if acceptanceSwitch == 0 || sum(tausProposed < 0) > 0
        noFail = 0;
    else
        noFail = 1;
    end
    
    %ensure on-off structure
    if acceptanceSwitch > 0 && noFail
        switchesRemoved = 0;
        if length(swTimesProposed) > 1
            if forceOnOff && brRatioTol > 0
                switchesRemoved = removeSameTypeOrRateSwitches(birthRatesProposed,degRateProposed,brRatioTol);
            elseif ~forceOnOff && brRatioTol > 0
                switchesRemoved = removeSameRateSwitches(birthRatesProposed,degRateProposed,brRatioTol);
            elseif forceOnOff && brRatioTol <= 0
                switchesRemoved = removeSameTypeSwitches(birthRatesProposed);
            end
        elseif length(swTimesProposed) == 1 && brRatioTol > 0
            switchesRemoved = removeSameRateSwitches(birthRatesProposed,degRateProposed,brRatioTol);
        end
        if switchesRemoved
            acceptanceSwitch = 0;
        end
    end
    
    %further checking of initial expression value and negative birth rates (result of previous checking - noFail)
    if acceptanceSwitch > 0
        acceptanceSwitch = acceptanceSwitch * (M_0Proposed >= 0) * (P_0Proposed >= minP_0) * (P_0 <= maxP_0) * noFail;
    end
    
    if rand(1) <= acceptanceSwitch
        swTimes = swTimesProposed;
        lLikelihood = lLikelihoodProposed;
    end
    
    [P_0 M_0 taus] = calculateModel_protein(timescale,swTimes,degRate_mRNA,degRate_protein,alpha,Y,timepoints,replicates,useDatapoint,Xmatrix,regression);
    [mRNA protein] = nStateSwitchODE_protein(timescale,M_0,taus,degRate_mRNA,swTimes,P_0,alpha,degRate_protein);
    %     mRNAFitHistory(n,:) = mRNA;
    %     proteinFitHistory(n,:) = protein;
    lLikelihoodHistory(n) = lLikelihood;
    %     M_0History(n) = M_0;
    %     P_0History(n) = P_0;
    swTimesHistory{n} = swTimes + timescaleShift;
    %     tausHistory{n} = taus;
    
    sigma2History(n) = 1 / prec;
    
    degRate_mRNAHistory(n) = degRate_mRNA;
    degRate_proteinHistory(n) = degRate_protein;
    alphaHistory(n) = alpha;
    
    if strcmp(regression.type,'weightedleastsquares')
        regressionWeightsPhiHistory(n) = regression.weightsPhi;
    end
    
    %check to see if we have hit the end of the burn in period
    if n == burnInIterations
        changeUp = cell(maxSwitches+1,1);
        changeUpBurnIn = zeros(maxSwitches+1,1);
        for x = 1:(maxSwitches+1)
            changeUp{x} = dimChangeUp( isnan(dimChangeUp(:,x)) == 0,x);
            validUp = dimChangeUp(1:burnInIterations,x);
            changeUpBurnIn(x) = mean( validUp(isnan(validUp) == 0) );
        end
        %changeUpBurnIn
        if dimChangeThreshold > 0
            earlyConverged = sum(changeUpBurnIn >= dimChangeThreshold) == 0;
        end
    end
    n = n + 1;
end

usableIterations = iterations - round(iterations * burnIn);
startIteration = iterations - usableIterations + 1;
%     samples = [swTimesHistory{ startIteration:end }];

profile = struct;
profile.profileName = geneLocus;
profile.modelType = 'protein';
profile.timescale = trueTimescale;
profile.timepoints = timepoints;
profile.replicates = replicates;
profile.data = Y;
profile.useDatapoint = useDatapoint;
profile.dataDetrend = detrend;

%construct parameters structure
params = struct;
params.iterations = iterations;
params.forceOnOff = forceOnOff;
params.maxSwitches = maxSwitches;
params.brRatioTolerance = brRatioTol;
params.burnIn = burnIn;
params.minSwitchStrength = strengthThreshold;
params.switchEvidence = dimChangeThreshold;
params.useInitialSwitchSuggest = suggestInitialSwitches;
params.initialSwitches = initialSwTimes;
params.initialSwitchThreshold = initialSwThreshold;
params.startIteration = startIteration;
params.percentiles = percentiles;
params.switchConfidences = confidences;
params.confidence = confidence;
params.estimationPoints = estimationPoints;
params.mRNADegradationRatePrior = prior_deg_m;
params.proteinDegradationRatePrior = prior_deg_p;
params.translationRatePrior = prior_alpha;
params.regressionMethod = regression.type;
if strcmpi(regression.type,'weightedleastsquares')
    params.regressionWeights = regression.weights;
end
params.logExplore = useLogExplore;
params.samplePrecision = explorePrec;
params.sampleAlpha = exploreAlpha;
params.switchNumberPrior = lambda;
params.switchDelay = switchDelay;
params.switchBaselineFactor = switchBaselineFactor;
params.maxSwitchDeviationFactor = maxSwitchDeviationFactor;
params.maskedSamples = maskedSamples;
params.totalSamplesUsed = sum(useDatapoint);
params.maxSamplesToUse = maxSamplesToUse;
params.maxSamplesToPlot = maxSamplesToPlot;

samples = struct;
samples.degradationRate_mRNAHistory = degRate_mRNAHistory;
samples.degradationRate_proteinHistory = degRate_proteinHistory;
samples.swTimesHistory = swTimesHistory;
% samples.tausHistory = tausHistory;
% samples.M_0History = M_0History;
% samples.P_0History = P_0History;
samples.logLikelihoodHistory = lLikelihoodHistory;
samples.sigma2History = sigma2History;
samples.alphaHistory = alphaHistory;
% samples.fit_mRNAHistory = mRNAFitHistory;
% samples.fit_proteinHistory = proteinFitHistory;
samples.lLikelihoodHistory = lLikelihoodHistory;
if strcmpi(regression.type,'weightedleastsquares')
    samples.regressionWeightsPhiHistory = regressionWeightsPhiHistory;
end

%construct posterior distributions
posteriors = extractPosteriorDistributions(profile,samples,params);

if isPlot
    if plotSamples
        plotSwitch(profile,posteriors,params,figHandle,samples);
    else
        plotSwitch(profile,posteriors,params,figHandle);
    end
    if plotDegDist
        figure('Name','mRNA degradation rate distributions','NumberTitle','off','WindowStyle','docked');
        t = 0.001:0.01:3;
        
        hold on;
        plot(t,gampdf(t,prior_deg_m(1), prior_deg_m(2)),'--k','linewidth',2);
        %if fit.posteriors.switchDistribution.switchNumber>0
        m = fit.posteriors.degradationRate_mRNA.mean;
        v = fit.posteriors.degradationRate_mRNA.standardDeviation^2;
        plot(t,gampdf(t,m*m / v, v / m),'b','linewidth',2);
        %end
        hold off;
        figure('Name','protein degradation rate distributions','NumberTitle','off','WindowStyle','docked');
        t = 0.001:0.01:3;
        
        hold on;
        plot(t,gampdf(t,prior_deg_p(1), prior_deg_p(2)),'--k','linewidth',2);
        %if fit.posteriors.switchDistribution.switchNumber>0
        m = fit.posteriors.degradationRate_protein.mean;
        v = fit.posteriors.degradationRate_protein.standardDeviation^2;
        plot(t,gampdf(t,m*m / v, v / m),'b','linewidth',2);
        %end
        hold off;
    end
end


pause(1);
end

% function [initialExp birthRates] = calculateProteinModelMissingData(timescale,swTimes,degRate,Y,timepoints,replicates,useSamples,Xmatrix)
% X = nStateSwitchModel(timescale,swTimes,degRate);
% %X = repmat(X,replicates,1);
% X1 = Xmatrix{length(swTimes)+1};
% for z = 1:replicates
%     X1( (timepoints * (z-1) + 1):(timepoints * z),:) = X;
% end
% X = X1(useSamples,:);
% Y = Y(useSamples);
% %least squares regression to obtain model initial exp and birth rates
% %alpha = inv(X' * X) * X' * Y;
% alpha = (X' * X) \ X' * Y;
%
% initialExp = alpha(1);
% birthRates = alpha(2:length(alpha)) * degRate;
% for m = 2:length(birthRates)
%     birthRates(m) = birthRates(m) + birthRates(m-1);
% end
%
% end

function [taus] = getTausHistory(means,sDevs,tausHistory,sTimesHistory)

taus = nan(length(tausHistory),length(means)+1);

swNum = length(means);

count = 1;
for x = 1:length(tausHistory)
    tau = tausHistory{x}';
    sTimes = sTimesHistory{x};
    if length(sTimes) == swNum
        if swNum == 0
            taus(count,:) = tau;
            count = count + 1;
        elseif(length(means) == length(sTimes)) && (length(means) == length(sDevs))
            %we know we have the same number of switches, now do those switches match with our fitted switches
            if sum( (sTimes >= (means - (1.96 .* sDevs))) .* (sTimes <= (means + (1.96 .* sDevs))) ) == swNum
                %we have the correct number of matching switches...include this set of birthRates
                taus(count,:) = tau;
                count = count + 1;
            end
        else
            disp('Incorrect mean/sDevs/sTimes vector');
            means
            sDevs
            sTimes
            swNum
            tau
        end
    end
end

taus = taus(1:(count-1),:);

end

function [switchesRemoved] = removeSameTypeOrRateSwitches(birthRates,degRate,tolerance)
nR = length(birthRates);
switchesRemoved = 0;

%check rates
adjustedRates = birthRates ./ degRate;
if sum(abs(adjustedRates(2:nR) - adjustedRates(1:(nR-1))) < tolerance) > 0
    switchesRemoved = 1;
    return;
end

%check if we have 'on-off' structure
%can we simplify and speed this up???
x = 2;
while ~switchesRemoved && x <= (nR-1)
    prevRate = birthRates(x-1);
    currRate = birthRates(x);
    nextRate = birthRates(x+1);
    if ((prevRate < currRate) && (currRate < nextRate)) || ((prevRate > currRate) && (currRate > nextRate))
        switchesRemoved = 1;
        return;
    end
    
    x = x + 1;
end

end

function [switchesRemoved] = removeSameTypeSwitches(birthRates)
nR = length(birthRates);
switchesRemoved = 0;
x = 2;
while ~switchesRemoved && x <= (nR-1)
    prevRate = birthRates(x-1);
    currRate = birthRates(x);
    nextRate = birthRates(x+1);
    %check if we have 'on-off' structure
    if ((prevRate < currRate) && (currRate < nextRate)) || ((prevRate > currRate) && (currRate > nextRate))
        switchesRemoved = 1;
        return;
    end
    
    x = x + 1;
end

end

function [switchesRemoved] = removeSameRateSwitches(birthRates,degRate,tolerance)
nR = length(birthRates);
switchesRemoved = 0;
adjustedRates = birthRates ./ degRate;
if sum(abs(adjustedRates(2:nR) - adjustedRates(1:(nR-1))) < tolerance) > 0
    switchesRemoved = 1;
    return;
end

end

function [y trend beta] = detrendLinear(data,timescale,replicates)

t = timescale ./ timescale(end);
t=t';
data=data';

%least squares regression
X  = [repmat(t,1,replicates) ones(length(t) * replicates,1)];
beta = (X' * X) \ X' * data;
trend = (t.*beta(1)) + beta(2);
y = data - trend;

end

function prior = calculateSwitchNumPrior(maxSwitches,lambda)
% prior on switch num is a truncated poisson distribution (capped at max switches)
prior = zeros(maxSwitches+1,1);

cumProb = 0;
for n = 0:maxSwitches
    prior(n+1) = exp(-lambda) * (lambda^n / factorial(n));
    cumProb = cumProb + prior(n+1);
end
remainder = 1 - cumProb;

%prior(end) = prior(end) + remainder;

% propToAdd = remainder / (maxSwitches + 1);
% for n = 0:maxSwitches
%     prior(n+1) = prior(n+1) + propToAdd;
% end
end

function [e_hat] = calculateEHat(profile,data,eHatMat,nReps,useDatapoint)
%[nReps timepoints] = size(data);
for x = 1:nReps
    eHatMat(x,:) = profile;
end
% e_hat = repmat(profile,nReps,1) - data;
e_hat = eHatMat - data;
e_hatT = e_hat';
e_hat = e_hatT(:);
e_hat = e_hat(useDatapoint);
end

function r = fastnormrnd(mu,sigma,numRands)
r = randn(1,numRands) .* sigma + mu;
end

function r = fastgamrnd(alpha,beta,numRands)
r = beta .* randg(alpha,1,numRands);
end

function i = isint(x)
i = (isnumeric(x) && isscalar(x)) .* (sum(round(x) == x) == length(x));
end

function [value] = getPercentile(data,percentile)
n = length(data);
%we know data is fine and sorted
%now get the percentiles
value = data(floor((n / 100) * percentile + 0.5));
end