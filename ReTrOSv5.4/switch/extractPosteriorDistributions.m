function posteriors = extractPosteriorDistributions(profileInfo,samples,params)
plotSwitchFitting = false;
posteriors = struct;

%29.07.2014
%TEMPORARY REMOVE
% if earlyConverged %|| isempty(samples)
%     disp('Fitting zero switches');
%
%     switchDist = struct;
%     switchDist.switchNumber = 0;
%     posteriors.switchDistribution = switchDist;
%
% else
%get out the params/chains needed
iterations = params.iterations;
startIteration = params.startIteration;
percentiles = params.percentiles;
confidences = params.switchConfidences;
maxSwitches = params.maxSwitches;
estimationPoints = params.estimationPoints;
switchDelay = params.switchDelay;
switchBaselineFactor = params.switchBaselineFactor;
strengthThreshold = params.minSwitchStrength;
regressionMethod = params.regressionMethod;
if strcmpi(regressionMethod,'weightedleastsquares')
    regressionWeights = params.regressionWeights;
    phiHistory = samples.regressionWeightsPhiHistory;
end
maxSwitchDeviationFactor = params.maxSwitchDeviationFactor;

if (iterations-startIteration+1) <= params.maxSamplesToUse
    sampleIndexToUse = 1:(iterations-startIteration+1);
else
    sampleIndexToUse = unique(round(linspace(1,iterations-startIteration+1,params.maxSamplesToUse)));
end
numSamplesUsedForPosteriors = length(sampleIndexToUse);

sigma2History = samples.sigma2History;
degRateHistory_mRNA = samples.degradationRate_mRNAHistory;
swTimesHistory = samples.swTimesHistory;
%fit_mRNAHistory = samples.fit_mRNAHistory;
lLikelihoodHistory = samples.lLikelihoodHistory;
%tausHistory = samples.tausHistory;

timescale = profileInfo.timescale;
deltaTimescale = max(timescale(2:end) - timescale(1:(end-1)));
timescaleShift = profileInfo.timescale(1);
timepoints = length(timescale);
replicates = profileInfo.replicates;
Y = profileInfo.data;
useDatapoint = profileInfo.useDatapoint;

%construct the default posteriors

%get the central X% of sigma (default: 90)
sigma2Struct = struct;
%sigmaHistory = sqrt(sigmaHistory);
sigma2HistoryToUse = sigma2History(startIteration:iterations);
sigma2HistoryToUse = sigma2HistoryToUse(sampleIndexToUse);
sortedSigma2 = sort(sigma2HistoryToUse);
sigma2Percentiles = zeros(length(percentiles),1);
for m = 1:length(percentiles)
    sigma2Percentiles(m) = getPercentile(sortedSigma2,percentiles(m));
end
m = getPercentile(sortedSigma2,50);
l = getPercentile(sortedSigma2,25);
u = getPercentile(sortedSigma2,75);
sigma2Struct.median = m;
sigma2Struct.lowerQuartile = l;
sigma2Struct.upperQuartile = u;
sigma2Struct.percentiles = sigma2Percentiles;
sigma2Struct.mean = mean(sortedSigma2);
sigma2Struct.standardDeviation = std(sortedSigma2);
posteriors.sigma2 = sigma2Struct;

if strfind(profileInfo.modelType,'protein')
    %get the central X% of alpha (default: 90)
    alphaStruct = struct;
    alphaHistoryToUse = samples.alphaHistory(startIteration:iterations);
    alphaHistoryToUse = alphaHistoryToUse(sampleIndexToUse);
    sortedAlpha = sort(alphaHistoryToUse);
    alphaPercentiles = zeros(length(percentiles),1);
    for m = 1:length(percentiles)
        alphaPercentiles(m) = getPercentile(sortedAlpha,percentiles(m));
    end
    m = getPercentile(sortedAlpha,50);
    l = getPercentile(sortedAlpha,25);
    u = getPercentile(sortedAlpha,75);
    alphaStruct.median = m;
    alphaStruct.lowerQuartile = l;
    alphaStruct.upperQuartile = u;
    alphaStruct.percentiles = alphaPercentiles;
    alphaStruct.mean = mean(sortedAlpha);
    alphaStruct.standardDeviation = std(sortedAlpha);
    posteriors.alpha = alphaStruct;
end

%get the central X% of degradation rates (mRNA) (default: 90)
degStruct = struct;
degRateHistoryToUse = degRateHistory_mRNA(startIteration:iterations);
degRateHistoryToUse = degRateHistoryToUse(sampleIndexToUse);
sortedDeg = sort(degRateHistoryToUse);
degPercentiles = zeros(length(percentiles),1);
for m = 1:length(percentiles)
    degPercentiles(m) = getPercentile(sortedDeg,percentiles(m));
end
m = getPercentile(sortedDeg,50);
l = getPercentile(sortedDeg,25);
u = getPercentile(sortedDeg,75);
degStruct.median = m;
degStruct.lowerQuartile = l;
degStruct.upperQuartile = u;
degStruct.percentiles = degPercentiles;
degStruct.mean = mean(sortedDeg);
degStruct.standardDeviation = std(sortedDeg);
p = gamfit(sortedDeg);
degStruct.alpha = p(1);
degStruct.beta = p(2);
posteriors.degradationRate_mRNA = degStruct;

%get the log of the degradation too
logDegStruct = struct;
sortedLogDeg = sort(log(degRateHistoryToUse));
logDegPercentiles = zeros(length(percentiles),1);
for m = 1:length(percentiles)
    logDegPercentiles(m) = getPercentile(sortedLogDeg,percentiles(m));
end
m = getPercentile(sortedLogDeg,50);
l = getPercentile(sortedLogDeg,25);
u = getPercentile(sortedLogDeg,75);
logDegStruct.median = m;
logDegStruct.lowerQuartile = l;
logDegStruct.upperQuartile = u;
logDegStruct.percentiles = logDegPercentiles;
logDegStruct.mean = mean(sortedLogDeg);
logDegStruct.standardDeviation = std(sortedLogDeg);
posteriors.logDegradationRate_mRNA = logDegStruct;

%     degpost = [l;m;u;mean(sortedLogDeg)-1.96*std(sortedLogDeg);mean(sortedLogDeg);mean(sortedLogDeg)+1.96*std(sortedLogDeg)]
%     figure;
%     plot(log(degRateHistory));
%     hold on;
%     plot([(iterations * burnIn)+1 iterations],[ posteriors.logDegradationRate.median posteriors.logDegradationRate.median], 'k');
%     plot([(iterations * burnIn)+1 iterations],[ posteriors.logDegradationRate.lowerQuartile posteriors.logDegradationRate.lowerQuartile], '--k');
%     plot([(iterations * burnIn)+1 iterations],[ posteriors.logDegradationRate.upperQuartile posteriors.logDegradationRate.upperQuartile], '--k');
%     %mean and standard deviation
%     plot([(iterations * burnIn)+1 iterations],[posteriors.logDegradationRate.mean posteriors.logDegradationRate.mean],'r');
%     plot([(iterations * burnIn)+1 iterations],[(posteriors.logDegradationRate.mean - confidence * posteriors.logDegradationRate.standardDeviation) (posteriors.logDegradationRate.mean - confidence * posteriors.logDegradationRate.standardDeviation)],'--r');
%     plot([(iterations * burnIn)+1 iterations],[(posteriors.logDegradationRate.mean + confidence * posteriors.logDegradationRate.standardDeviation) (posteriors.logDegradationRate.mean + confidence * posteriors.logDegradationRate.standardDeviation)],'--r');
%     hold off;


halfLifeStruct = struct;
sortedHL = sort(log(2) ./ degRateHistoryToUse);
hlPercentiles = zeros(length(percentiles),1);
for m = 1:length(percentiles)
    hlPercentiles(m) = getPercentile(sortedHL,percentiles(m));
end
m = getPercentile(sortedHL,50);
l = getPercentile(sortedHL,25);
u = getPercentile(sortedHL,75);
halfLifeStruct.median = m;
halfLifeStruct.lowerQuartile = l;
halfLifeStruct.upperQuartile = u;
halfLifeStruct.percentiles = hlPercentiles;
halfLifeStruct.mean = mean(sortedHL);
halfLifeStruct.standardDeviation = std(sortedHL);
p = gamfit(sortedHL);
halfLifeStruct.alpha = p(1);
halfLifeStruct.beta = p(2);
posteriors.halfLife_mRNA = halfLifeStruct;
if strfind(profileInfo.modelType,'protein')
    %get the central X% of degradation rates (default: 90)
    degStruct = struct;
    degRateHistoryToUse = samples.degradationRate_proteinHistory(startIteration:iterations);
    degRateHistoryToUse = degRateHistoryToUse(sampleIndexToUse);
    sortedDeg = sort(degRateHistoryToUse);
    degPercentiles = zeros(length(percentiles),1);
    for m = 1:length(percentiles)
        degPercentiles(m) = getPercentile(sortedDeg,percentiles(m));
    end
    m = getPercentile(sortedDeg,50);
    l = getPercentile(sortedDeg,25);
    u = getPercentile(sortedDeg,75);
    degStruct.median = m;
    degStruct.lowerQuartile = l;
    degStruct.upperQuartile = u;
    degStruct.percentiles = degPercentiles;
    degStruct.mean = mean(sortedDeg);
    degStruct.standardDeviation = std(sortedDeg);
    p = gamfit(sortedDeg);
    degStruct.alpha = p(1);
    degStruct.beta = p(2);
    posteriors.degradationRate_protein = degStruct;
    
    %get the log of the degradation too
    logDegStruct = struct;
    sortedLogDeg = sort(log(degRateHistoryToUse));
    logDegPercentiles = zeros(length(percentiles),1);
    for m = 1:length(percentiles)
        logDegPercentiles(m) = getPercentile(sortedLogDeg,percentiles(m));
    end
    m = getPercentile(sortedLogDeg,50);
    l = getPercentile(sortedLogDeg,25);
    u = getPercentile(sortedLogDeg,75);
    logDegStruct.median = m;
    logDegStruct.lowerQuartile = l;
    logDegStruct.upperQuartile = u;
    logDegStruct.percentiles = logDegPercentiles;
    logDegStruct.mean = mean(sortedLogDeg);
    logDegStruct.standardDeviation = std(sortedLogDeg);
    posteriors.logDegradationRate_protein = logDegStruct;
    
    halfLifeStruct = struct;
    sortedHL = sort(log(2) ./ degRateHistoryToUse);
    hlPercentiles = zeros(length(percentiles),1);
    for m = 1:length(percentiles)
        hlPercentiles(m) = getPercentile(sortedHL,percentiles(m));
    end
    m = getPercentile(sortedHL,50);
    l = getPercentile(sortedHL,25);
    u = getPercentile(sortedHL,75);
    halfLifeStruct.median = m;
    halfLifeStruct.lowerQuartile = l;
    halfLifeStruct.upperQuartile = u;
    halfLifeStruct.percentiles = hlPercentiles;
    halfLifeStruct.mean = mean(sortedHL);
    halfLifeStruct.standardDeviation = std(sortedHL);
    p = gamfit(sortedHL);
    halfLifeStruct.alpha = p(1);
    halfLifeStruct.beta = p(2);
    posteriors.halfLife_protein = halfLifeStruct;
end

%switch frequencies (non-parametric)
swTimesHistoryToUse = swTimesHistory(startIteration:iterations);
swTimesHistoryToUse = swTimesHistoryToUse(sampleIndexToUse);


switchCounts = struct;
swCounts = zeros(maxSwitches+1,1);
swNum = zeros(numSamplesUsedForPosteriors,1);
for x = 1:numSamplesUsedForPosteriors
    currS = swTimesHistoryToUse{x};
    swCounts(length(currS) +1) = swCounts( length(currS) + 1) + 1;
    swNum(x) = length(currS);
end
%     switchCounts.frequencies = swCounts ./ usableIterations;
switchCounts.frequencies = swCounts ./ numSamplesUsedForPosteriors;
switchCounts.mean = mean(swNum);
switchCounts.standardDeviation = std(swNum);
swNum = sort(swNum);
swNumPercentiles = zeros(length(percentiles),1);
for m = 1:length(percentiles)
    swNumPercentiles(m) = getPercentile(swNum,percentiles(m));
end
m = getPercentile(swNum,50);
l = getPercentile(swNum,25);
u = getPercentile(swNum,75);
switchCounts.median = m;
switchCounts.lowerQuartile = l;
switchCounts.upperQuartile = u;
switchCounts.percentiles = swNumPercentiles;
posteriors.switchFrequency = switchCounts;

%construct the (marginal) switch distributions
%we'll fit both a gaussian mixture model and individual recombined gaussians

%GMM
%[f xi baseline mus sigmas heights] = getSwitchFitNonLinear(swTimesHistory,startIteration,trueTimescale,estimationPoints,switchDelay,true,plotSwitchFitting);


[f, xi, h, baseline, mus, sigmas weights] = getSwitchFit(swTimesHistoryToUse,timescale,estimationPoints,switchDelay,switchBaselineFactor,true,plotSwitchFitting);

switchDist = struct;
%     switchDist.function = f;
%     switchDist.estimationPoints = xi;
%     switchDist.baseline = baseline;
%need to generate switch strength (the percentage of iterations which sampled each switch with the given confidence
% actually we want to automatically compute these for each of the 3 sd levels
sampled = zeros(length(mus),length(confidences));
for z = 1:length(confidences)
    swTimesL = mus - confidences(z) * sigmas;
    swTimesU = mus + confidences(z) * sigmas;
    %for x = startIteration:iterations
    for x = 1:numSamplesUsedForPosteriors
        currS = swTimesHistoryToUse{x};
        for y = 1:length(mus)
            sampled(y,z) = sampled(y,z) + (sum((swTimesL(y) <= currS) .* (swTimesU(y) >= currS)) > 0);
        end
    end
end
%     strengths = sampled ./ usableIterations;
strengths = sampled ./ numSamplesUsedForPosteriors;
%once we have the strengths we want to threshold on these
%use the 1.96SD strengths (~95%)
keep = ((strengths(:,2) >= strengthThreshold) .* (sigmas <= (maxSwitchDeviationFactor * deltaTimescale))) == 1;
keepSwitchNum = sum(keep);
switchDist.switchNumber = keepSwitchNum;
switchDist.means = mus(keep);
switchDist.standardDeviations = sigmas(keep);
switchDist.strengths = strengths(keep,:);
switchDist.weights = weights(keep);
switchDist.function = f;
switchDist.bandwidth = h;
posteriors.switchDistribution = switchDist;

%construct dimension-specific switch models
modelSwitchDist = struct;
modelSwitchDist.models = getModelsBasic(switchCounts.frequencies,swTimesHistoryToUse,timescale, estimationPoints,switchDelay,switchBaselineFactor);
posteriors.switchModelsBasic = modelSwitchDist;

%         %individual
%         [f xi baseline mus sigmas heights] = getSwitchFitDelay(swTimesHistory,startIteration,trueTimescale,estimationPoints,switchDelay,false,false);
%         switchDist = struct;
%         %     switchDist.function = f;
%         %     switchDist.estimationPoints = xi;
%         %     switchDist.baseline = baseline;
%         %need to generate switch strength (the percentage of iterations which sampled each switch with the given confidence
%         % actually we want to automatically compute these for each of the 3 sd levels
%         sampled = zeros(length(mus),length(confidences));
%         for z = 1:length(confidences)
%             swTimesL = mus - confidences(z) * sigmas;
%             swTimesU = mus + confidences(z) * sigmas;
%             for x = startIteration:iterations
%                 currS = swTimesHistory{x};
%                 for y = 1:length(mus)
%                     sampled(y,z) = sampled(y,z) + (sum((swTimesL(y) <= currS) .* (swTimesU(y) >= currS)) > 0);
%                 end
%             end
%         end
%         strengths = sampled ./ usableIterations;
%         %once we have the strengths we want to threshold on these
%         %use the 1.96SD strengths (~95%)
%         keep = strengths(:,2) >= strengthThreshold;
%         keepSwitchNum = sum(keep);
%         switchDist.switchNumber = keepSwitchNum;
%         switchDist.means = mus(keep);
%         switchDist.standardDeviations = sigmas(keep);
%         switchDist.strengths = strengths(keep,:);
%         switchDist.heights = heights(keep);
%         posteriors.switchDistributionInd = switchDist;



% %reconstruct fit and create distribution
% fitDistribution = struct;
%
% %get distributions
% fit_mRNAHistoryToUse = fit_mRNAHistory(startIteration:iterations,:);
% fit_mRNAHistoryToUse = fit_mRNAHistoryToUse(sampleIndexToUse,:);
% fitMean = mean(fit_mRNAHistoryToUse);
% fitSigma = std(fit_mRNAHistoryToUse);
% sFit = sort(fit_mRNAHistoryToUse);
% fitPercentiles = zeros(length(percentiles),timepoints);
% for x = 1:timepoints
%     for m = 1:length(percentiles)
%         fitPercentiles(m,x) = getPercentile(sFit(:,x),percentiles(m));
%     end
% end
% fitMedian = zeros(timepoints,1);
% fitUpperQ = zeros(timepoints,1);
% fitLowerQ = zeros(timepoints,1);
% for x = 1:timepoints
%     fitMedian(x) = getPercentile(sFit(:,x),50);
%     fitLowerQ(x) = getPercentile(sFit(:,x),25);
%     fitUpperQ(x) = getPercentile(sFit(:,x),75);
% end
% fitDistribution.mean = fitMean;
% fitDistribution.standardDeviation = fitSigma;
% fitDistribution.median = fitMedian;
% fitDistribution.lowerQuartile = fitLowerQ;
% fitDistribution.upperQuartile = fitUpperQ;
% fitDistribution.percentiles = fitPercentiles;
% posteriors.fit_mRNADistribution = fitDistribution;
%
% if strcmpi(gene.modelType,'protein')
%     %reconstruct fit and create distribution
%     fitDistribution = struct;
%
%     %get distributions
%     fit_proteinHistoryToUse = samples.fit_proteinHistory(startIteration:iterations,:);
%     fit_proteinHistoryToUse = fit_proteinHistoryToUse(sampleIndexToUse,:);
%     fitMean = mean(fit_proteinHistoryToUse);
%     fitSigma = std(fit_proteinHistoryToUse);
%     sFit = sort(fit_proteinHistoryToUse);
%     fitPercentiles = zeros(length(percentiles),timepoints);
%     for x = 1:timepoints
%         for m = 1:length(percentiles)
%             fitPercentiles(m,x) = getPercentile(sFit(:,x),percentiles(m));
%         end
%     end
%     fitMedian = zeros(timepoints,1);
%     fitUpperQ = zeros(timepoints,1);
%     fitLowerQ = zeros(timepoints,1);
%     for x = 1:timepoints
%         fitMedian(x) = getPercentile(sFit(:,x),50);
%         fitLowerQ(x) = getPercentile(sFit(:,x),25);
%         fitUpperQ(x) = getPercentile(sFit(:,x),75);
%     end
%     fitDistribution.mean = fitMean;
%     fitDistribution.standardDeviation = fitSigma;
%     fitDistribution.median = fitMedian;
%     fitDistribution.lowerQuartile = fitLowerQ;
%     fitDistribution.upperQuartile = fitUpperQ;
%     fitDistribution.percentiles = fitPercentiles;
%     posteriors.fit_proteinDistribution = fitDistribution;
% end

fitStruct = struct;

regression.type = regressionMethod;

if strcmpi(regressionMethod,'weightedleastsquares')
    phiStruct = struct;
    %sigmaHistory = sqrt(sigmaHistory);
    phiHistoryToUse = phiHistory(startIteration:iterations);
    phiHistoryToUse = phiHistoryToUse(sampleIndexToUse);
    sortedPhi = sort(phiHistoryToUse);
    phiPercentiles = zeros(length(percentiles),1);
    for m = 1:length(percentiles)
        phiPercentiles(m) = getPercentile(sortedPhi,percentiles(m));
    end
    m = getPercentile(sortedPhi,50);
    l = getPercentile(sortedPhi,25);
    u = getPercentile(sortedPhi,75);
    phiStruct.median = m;
    phiStruct.lowerQuartile = l;
    phiStruct.upperQuartile = u;
    phiStruct.percentiles = phiPercentiles;
    phiStruct.mean = mean(sortedPhi);
    phiStruct.standardDeviation = std(sortedPhi);
    posteriors.regressionWeightsPhi = phiStruct;
    
    %create Q from median of sampled distribution
    q = repmat(regressionWeights.^(-phiStruct.median),replicates,1);
    regression.Q = diag(q);
end

if strfind(profileInfo.modelType,'mRNA')
    %birth rates and initial exp
    maxComponents = length(posteriors.switchDistribution.means);
    Xmatrix = cell(maxComponents+1,1);
    for x = 1:(maxComponents+1)
        Xmatrix{x} = zeros(timepoints * replicates,x+1);
    end
    [initialExp, birthRates] = calculateModel_mRNA(timescale-timescaleShift,posteriors.switchDistribution.means-timescaleShift,posteriors.degradationRate_mRNA.mean,Y,timepoints,replicates,useDatapoint,Xmatrix,regression);
    fitStruct.M_0_meanModel = initialExp;
    fitStruct.taus_meanModel = birthRates;
    [initialExp, birthRates] = calculateModel_mRNA(timescale-timescaleShift,posteriors.switchDistribution.means-timescaleShift,posteriors.degradationRate_mRNA.median,Y,timepoints,replicates,useDatapoint,Xmatrix,regression);
    fitStruct.M_0_medianModel = initialExp;
    fitStruct.taus_medianModel = birthRates;
else
    %birth rates and initial exp
    maxComponents = length(posteriors.switchDistribution.means);
    Xmatrix = cell(maxComponents+1,1);
    for x = 1:(maxComponents+1)
        Xmatrix{x} = zeros(timepoints * replicates,x+2);
    end
    [P_0, M_0, birthRates] = calculateModel_protein(timescale-timescaleShift,posteriors.switchDistribution.means-timescaleShift,posteriors.degradationRate_mRNA.mean,posteriors.degradationRate_protein.mean,posteriors.alpha.mean,Y,timepoints,replicates,useDatapoint,Xmatrix,regression);
    fitStruct.M_0_meanModel = M_0;
    fitStruct.P_0_meanModel = P_0;
    fitStruct.taus_meanModel = birthRates;
    [P_0, M_0, birthRates] = calculateModel_protein(timescale-timescaleShift,posteriors.switchDistribution.means-timescaleShift,posteriors.degradationRate_mRNA.median,posteriors.degradationRate_protein.median,posteriors.alpha.median,Y,timepoints,replicates,useDatapoint,Xmatrix,regression);
    fitStruct.M_0_medianModel = M_0;
    fitStruct.P_0_medianModel = P_0;
    fitStruct.taus_medianModel = birthRates;
end

posteriors.modelFit = fitStruct;
% tausHistoryToUse = tausHistory(startIteration:iterations);
% tausHistoryToUse = tausHistoryToUse(sampleIndexToUse);
% specificTaus = getTausHistory(posteriors.switchDistribution.means',posteriors.switchDistribution.standardDeviations',tausHistoryToUse,swTimesHistoryToUse);
%
% brMeans = nan(posteriors.switchDistribution.switchNumber+1,1);
% brSigmas = nan(posteriors.switchDistribution.switchNumber+1,1);
% brPercentiles = nan(posteriors.switchDistribution.switchNumber+1,length(percentiles));
% brMedians = nan(posteriors.switchDistribution.switchNumber+1,1);
% brLowerQs = nan(posteriors.switchDistribution.switchNumber+1,1);
% brUpperQs = nan(posteriors.switchDistribution.switchNumber+1,1);
%
% if size(specificTaus,1) > 100
%     for x = 1:(posteriors.switchDistribution.switchNumber+1)
%         usableBirthRates = specificTaus(:,x);
%         brMeans(x) = mean(usableBirthRates);
%         brSigmas(x) = std(usableBirthRates);
%         sBR = sort(usableBirthRates);
%         for m = 1:length(percentiles)
%             brPercentiles(x,m) = getPercentile(sBR,percentiles(m));
%         end
%         brMedians(x) = getPercentile(sBR,50);
%         brLowerQs(x) = getPercentile(sBR,25);
%         brUpperQs(x) = getPercentile(sBR,75);
%     end
% end
% birthRatesStruct.mean = brMeans;
% birthRatesStruct.standardDeviation = brSigmas;
% birthRatesStruct.median = brMedians;
% birthRatesStruct.lowerQuartile = brLowerQs;
% birthRatesStruct.upperQuartile = brUpperQs;
% birthRatesStruct.percentile = brPercentiles;
%
% posteriors.taus = birthRatesStruct;
% posteriors.M_0 = M_0Struct;
% if strcmpi(gene.modelType,'protein');
%     posteriors.P_0 = P_0Struct;
% end

%log likelihood
likStruct = struct;
lLikelihoodHistoryToUse = lLikelihoodHistory(startIteration:iterations);
lLikelihoodHistoryToUse = lLikelihoodHistoryToUse(sampleIndexToUse);
sortedLik = sort(lLikelihoodHistoryToUse);
likPercentiles = zeros(length(percentiles),1);
for m = 1:length(percentiles)
    likPercentiles(m) = getPercentile(sortedLik,percentiles(m));
end
m = getPercentile(sortedLik,50);
l = getPercentile(sortedLik,25);
u = getPercentile(sortedLik,75);
likStruct.median = m;
likStruct.lowerQuartile = l;
likStruct.upperQuartile = u;
likStruct.percentiles = likPercentiles;
likStruct.mean = mean(sortedLik);
likStruct.standardDeviation = std(sortedLik);
posteriors.logLikelihood = likStruct;

%29.07.2014
%TEMPORARY REMOVE
%     changeUp = cell(maxSwitches+1,1);
%     changeDown = cell(maxSwitches+1,1);
%     changeUpBurnIn = zeros(maxSwitches+1,1);
%     changeDownBurnIn = zeros(maxSwitches+1,1);
%     changeUpMean = zeros(maxSwitches+1,1);
%     changeDownMean = zeros(maxSwitches+1,1);
%     for x = 1:(maxSwitches+1)
%         changeUp{x} = dimChangeUp( isnan(dimChangeUp(:,x)) == 0,x);
%         validUp = dimChangeUp(1:(startIteration-1),x);
%         changeUpBurnIn(x) = mean( validUp(isnan(validUp) == 0) );
%         validUp = dimChangeUp(startIteration:iterations,x);
%         changeUpMean(x) = mean( validUp(isnan(validUp) == 0) );
%         changeDown{x} = dimChangeDown( isnan(dimChangeDown(:,x)) == 0,x);
%         validDown = dimChangeDown(1:(startIteration-1),x);
%         changeDownBurnIn(x) = mean( validDown(isnan(validDown) == 0) );
%         validDown = dimChangeDown(startIteration:iterations,x);
%         changeDownMean(x) = mean( validDown(isnan(validDown) == 0) );
%     end
%     samples.dimChangeUp = changeUp;
%     samples.dimChangeDown = changeDown;
%     %     fit.dimChangeUpBurnIn = changeUpBurnIn;
%     %     fit.dimChangeDownBurnIn = changeDownBurnIn;
%     %     fit.dimChangeUpMean = changeUpMean;
%     %     fit.dimChangeDownMean = changeDownMean;
%     samples.dimChangeUpBurnIn = changeUpBurnIn;
%     samples.dimChangeDownBurnIn = changeDownBurnIn;
%     samples.dimChangeUpMean = changeUpMean;
%     samples.dimChangeDownMean = changeDownMean;

%end

end

function [value] = getPercentile(data,percentile)
n = length(data);
%we know data is fine and sorted
%now get the percentiles
value = data(floor((n / 100) * percentile + 0.5));
end

function [birthRates] = getTausHistory(means,sDevs,birthRatesHistory,sTimesHistory)

birthRates = nan(length(birthRatesHistory),length(means)+1);

swNum = length(means);

count = 1;
for x = 1:length(birthRatesHistory)
    br = birthRatesHistory{x}';
    sTimes = sTimesHistory{x};
    if length(sTimes) == swNum
        if swNum == 0
            birthRates(count,:) = br;
            count = count + 1;
        elseif(length(means) == length(sTimes)) && (length(means) == length(sDevs))
            %we know we have the same number of switches, now do those switches match with our fitted switches
            if sum( (sTimes >= (means - (1.96 .* sDevs))) .* (sTimes <= (means + (1.96 .* sDevs))) ) == swNum
                %we have the correct number of matching switches...include this set of birthRates
                birthRates(count,:) = br;
                count = count + 1;
            end
        else
            disp('Incorrect mean/sDevs/sTimes vector');
            means
            sDevs
            sTimes
            swNum
            br
        end
    end
end

birthRates = birthRates(1:(count-1),:);

end