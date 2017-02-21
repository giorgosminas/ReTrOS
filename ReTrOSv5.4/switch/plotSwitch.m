function [] = plotSwitch(varargin)

if length(varargin) == 4
    fullFit = false;
elseif length(varargin) == 5
    fullFit = true;
else
    error('plotSwitch:argCheck','Invalid input parameters');
end
profile = varargin{1};
posteriors = varargin{2};
params = varargin{3};
figHandle = varargin{4};

if isempty(figHandle)
    figHandle = figure;
end

% timepointsUsed = fit.algorithmParameters.timepointsUsed;
% replicatesUsed = fit.algorithmParameters.replicatesUsed;
name = profile.profileName;
timescale = profile.timescale;
timepoints = length(timescale);
replicates = profile.replicates;
maskedSamples = params.maskedSamples;
% averageReplicates = fit.algorithmParameters.averageReplicates;
maxSwitches = params.maxSwitches;
confidence = params.confidence;
iterations = params.iterations;
burnIn = params.burnIn;
startIteration = params.startIteration;
percentiles = params.percentiles;
estimationPoints = params.estimationPoints;
maxSamplesToPlot = params.maxSamplesToPlot;
modelType = profile.modelType;
Y = profile.data;
timeResolution = 100;

regression.type = params.regressionMethod;
if strcmpi(regression.type,'weightedleastsquares')
    regression.weights = params.regressionWeights;
    q = repmat(regression.weights.^(-posteriors.regressionWeightsPhi.median),replicates,1);
    regression.Q = diag(q);
end

%new parameter determining which specific samples to use (defaults to all, but allows the user to mark specific samples either as missing or outliers)
useDatapoint = profile.useDatapoint;

% useSample(maskedSamples) = 0;
% maskedD = nan(1,timepoints*replicates);
% maskedD(maskedSamples) = data(maskedSamples);
% maskedData = reshape(maskedD,timepoints,replicates)';

Y(useDatapoint==0) = NaN;

data = reshape(Y',timepoints,replicates)';

if replicates == 1
    medianProfile = data;
else
    %due to the missing values, we'll have to do it on a timepoint basis
    medianProfile = median(data);
    for m = find(isnan(medianProfile))
        medianProfile(m) = median(data(data(:,m)>=0,m));
    end
end

%get ylims
minData = min(min(data));
maxData = max(max(data));
yLims = [minData * 0.9, maxData * 1.1];

lowerData = zeros(1,timepoints);
upperData = zeros(1,timepoints);
for x = 1:timepoints
    lowerData(x) = min(data(:,x));
    upperData(x) = max(data(:,x));
end

%extract out data we need
swTimes = posteriors.switchDistribution.means;
swSigmas = posteriors.switchDistribution.standardDeviations;
swWeights = posteriors.switchDistribution.weights;
swFunction = posteriors.switchDistribution.function;
% swFunction = fit.posteriors.switchDistribution.function;
swEstimationPoints = linspace(timescale(1),timescale(end),estimationPoints);%fit.posteriors.switchDistribution.estimationPoints;
% swBaseline = fit.posteriors.switchDistribution.baseline;
% fitMean = fit.posteriors.fitDistribution.mean;
% fitStd = fit.posteriors.fitDistribution.standardDeviation;
%do we have the full fit (inc. chain samples) or just the posteriors


if strfind(modelType,'mRNA')
    if fullFit
        samples = varargin{5};
        swTimesHistory = samples.swTimesHistory;
        %         lLikelihood = posteriors.logLikelihood;
        %         lLikelihoodHistory = samples.logLikelihoodHistory;
        degRateHistory = samples.degradationRate_mRNAHistory;
        %initialExpHistory = fit.samples.initialExpHistory;
        sigma2History = samples.sigma2History;
        mainPlotCols = 3;
        cols = mainPlotCols + 1;
        rows = 3;
    else
        mainPlotCols = 1;
        cols = 1;
        rows = 2;
    end
    %adjustedTimescale = sort([timescale swTimes']);
    %adjustedProfile = nStateSwitchODE(adjustedTimescale,initialExp,birthRates,degRate,swTimes);
    
    %get the names of the profile
    
    %scrsz = get(0,'ScreenSize');
    %figure('Name',['ReTrOS-switch: ' name],'NumberTitle','off','Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)*0.5 scrsz(4)*0.5],'WindowStyle','docked');
    
    figure(figHandle);
    clf(figHandle);
    set(figHandle,'Name',['ReTrOS-switch: ' name],'NumberTitle','off');
    
    mainArea = 1:mainPlotCols;
    subplot(rows,cols,mainArea);
    dataCol = 'b';
    fillCol = [0.65 0.65 1];
    %fill function plots a single filled polygon, such that x1 is connected to xn (last point in vector)
    %therefore the last half of the x vector is reversed to do the top side of the polygon, and upper bound of data is reversed
    
    %check to see if we ignore all data for single timepoint
    hasData = lowerData>=0;
    if replicates > 1
        patch = fill( [timescale(hasData) fliplr(timescale(hasData))], [lowerData(hasData) fliplr(upperData(hasData))], fillCol); %final parameter is colour
        set(patch,'EdgeColor','none');%,'FaceAlpha',0.5); % makes the polygon transparent
    end
    hold on;
    %plot median profile
    plot(timescale,medianProfile,dataCol,'linewidth',2);
    for x = 1:replicates
        for y = 1:timepoints
            plot(timescale(y),data(x,y),[dataCol 'o'],'MarkerSize',4);
        end
    end
    
    %plot the ignored data
    % for x = 1:replicates
    %     for y = 1:timepoints
    %         plot(timescale(y),maskedData(x,y),'ko','MarkerSize',4);
    %     end
    % end
    
    title(name,'Interpreter','none');
    
    
    % %plot model profile
    
    %plot(sort([timescale swTimes']),nStateSwitchODE(sort([timescale swTimes'])-timescale(1), fit.posteriors.initialExpression.medianModel, fit.posteriors.birthRates.medianModel, fit.posteriors.degradationRate.median,swTimes-timescale(1)),'r','linewidth',2);
    Xmatrix = cell(length(swTimes)+1,1);
    for x = 1:(length(swTimes)+1)
        Xmatrix{x} = zeros(timepoints * replicates,x+1);
    end
    [M_0, taus] = calculateModel_mRNA(timescale-timescale(1),swTimes-timescale(1),posteriors.degradationRate_mRNA.median,Y,timepoints,replicates,useDatapoint,Xmatrix,regression);
    timeSw = unique([linspace(timescale(1),timescale(end),timeResolution) swTimes' timescale(end)]);
    plot(timeSw,nStateSwitchODE_mRNA(timeSw-timescale(1), M_0, taus, posteriors.degradationRate_mRNA.median,swTimes-timescale(1)),'k','linewidth',2);

    %plot switches
    for x = 1:length(swTimes)
        if taus(x+1) > taus(x)
            swCol = 'r';
        else
            swCol = 'g';
        end
        plot([swTimes(x) swTimes(x) ], [minData maxData],swCol,'LineWidth',2)
        plot([swTimes(x) - confidence*swSigmas(x) swTimes(x) - confidence*swSigmas(x) ], [minData maxData],['--' swCol],'linewidth',2)
        plot([swTimes(x) + confidence*swSigmas(x) swTimes(x) + confidence*swSigmas(x) ], [minData maxData],['--' swCol],'linewidth',2)
    end
    xlim([timescale(1) timescale(end)]);
    ylim(yLims);
    %xlabel('Timescale (hours)');
    ylabel('mRNA expression');
    hold off;
    
    
    %swiches PDF
    
    mainArea = (cols+1):(cols+mainPlotCols);
    subplot(rows,cols,mainArea);
    %plot the switch fits (gaussian)
    switchNum = length(swTimes);
    
    sFit = normFunc(swEstimationPoints,swTimes,swSigmas,swWeights);
    for z = 1:switchNum
        %colour the portion of the gaussian used to calculate switch size/strength
        l = swTimes(z) - confidence * swSigmas(z);
        u = swTimes(z) + confidence * swSigmas(z);
        idx = (swEstimationPoints >= l) .* (swEstimationPoints <= u);
        x1 = swEstimationPoints(idx == 1);
        x = [x1 fliplr(x1)];
        y1 = sFit(idx == 1)';
        y = [y1 zeros(1,length(x1))];
        
        if taus(z+1) > taus(z)
            swCol = [1 0.5 0.5];
        else
            swCol = [0.5 1 0.5];
        end
        p = fill(x, y,swCol);
        set(p,'edgecolor','none');%,'facealpha',0.5);
        hold on;
    end
    plot(swEstimationPoints,sFit,'k','linewidth',2);
    plot(swEstimationPoints,swFunction,'b','linewidth',2);
    hold off;
    ylabel('Switch PDF');
    xlim([timescale(1) timescale(end)]);
    %     xlabel('Timescale in hours');
    
    
    if fullFit
        if iterations <= maxSamplesToPlot
            %switch samples
            mainArea = (cols*2+1):(cols*2+mainPlotCols);
            subplot(rows,cols,mainArea);
            
            swData = cell(maxSwitches,1);
            for x = 1:maxSwitches
                swData{x} = zeros(iterations,1) * NaN;
            end
            
            for x = 1:iterations
                %plot( ones(length(swTimesHistory{x}),1) * x,swTimesHistory{x},'o','MarkerSize',1);
                currSwData = swTimesHistory{x};
                for y = 1:length(currSwData)
                    swData{y}(x) = currSwData(y);
                end
            end
            for x = 1:maxSwitches
                plot(swData{x},'ob','MarkerSize',1);
                hold on;
            end
            
            %NOTE: we plot everything the 'wrong way round' and then rotate for efficient speed of plotting
            plot([startIteration-1 startIteration-1],[timescale(1) timescale(end)],'--k','linewidth',2);
            if iterations > 1
                xlim([1 iterations]);
            end
            
            %plot X% for each switch (+/- confidence SDs)
            for x = 1:length(swTimes)
                if taus(x+1) > taus(x)
                    swCol = 'r';
                else
                    swCol = 'g';
                end
                plot([startIteration iterations],[swTimes(x) swTimes(x)],swCol,'linewidth',2);
                plot([startIteration iterations],[swTimes(x) - confidence*swSigmas(x) swTimes(x) - confidence*swSigmas(x)],['--' swCol],'linewidth',2);
                plot([startIteration iterations],[swTimes(x) + confidence*swSigmas(x) swTimes(x) + confidence*swSigmas(x)],['--' swCol],'linewidth',2);
            end
            view(-90,90);
            set(gca,'ydir','reverse');
            ylim([timescale(1) timescale(end)]);
            ylabel('Timescale (hours)');
            xlabel('Switch samples');
            hold off;
            
            
            %             subplot(rows,cols,cols);
            %             %just plot likelihood history
            %             plot(lLikelihoodHistory);
            %             hold on;
            %             maxL = max(lLikelihoodHistory(startIteration:iterations));
            %             minL = min(lLikelihoodHistory(startIteration:iterations));
            %             plot([startIteration-1 startIteration-1],[minL maxL],'--k','linewidth',2);
            %             %median + percentiles
            %             plot([startIteration iterations],[lLikelihood.median lLikelihood.median],'k','linewidth',2);
            %             plot([startIteration iterations],[lLikelihood.lowerQuartile lLikelihood.lowerQuartile],'--k','linewidth',2);
            %             plot([startIteration iterations],[lLikelihood.upperQuartile lLikelihood.upperQuartile],'--k','linewidth',2);
            %             %mean + standard deviation
            %             plot([startIteration iterations],[lLikelihood.mean lLikelihood.mean],'r','linewidth',2);
            %             plot([startIteration iterations],[(lLikelihood.mean - confidence * lLikelihood.standardDeviation) (lLikelihood.mean - confidence * lLikelihood.standardDeviation)],'--r','linewidth',2);
            %             plot([startIteration iterations],[(lLikelihood.mean + confidence * lLikelihood.standardDeviation) (lLikelihood.mean + confidence * lLikelihood.standardDeviation)],'--r','linewidth',2);
            %             hold off;
            %             title('Log Likelihood');
            %             %set the y limits based on burn in period
            %             minLMult = 0.9;
            %             if minL < 0
            %                 minLMult = 1.1;
            %             end
            %             maxLMult = 1.1;
            %             if maxL < 0
            %                 maxLMult = 0.9;
            %             end
            %             ylim([minL*minLMult maxL*maxLMult]);
            %             %set(gca,'OuterPosition',[mainPlotWidth 0.5 (1 / cols) 0.5]);
            
            %sigma history
            subplot(rows,cols,cols);
            if strcmpi(params.regressionMethod,'leastsquares')
                plot(sigma2History);
                hold on;
                maxS = max(sigma2History(startIteration:iterations));
                minS = min(sigma2History(startIteration:iterations));
                plot([startIteration-1 startIteration-1],[minS maxS],'--k','linewidth',2);
                %median + percentiles
                plot([startIteration iterations],[posteriors.sigma2.median posteriors.sigma2.median],'k','linewidth',2);
                plot([startIteration iterations],[posteriors.sigma2.lowerQuartile posteriors.sigma2.lowerQuartile],'--k','linewidth',2);
                plot([startIteration iterations],[posteriors.sigma2.upperQuartile posteriors.sigma2.upperQuartile],'--k','linewidth',2);
                %mean + standard deviation
                plot([startIteration iterations],[posteriors.sigma2.mean posteriors.sigma2.mean],'r','linewidth',2);
                plot([startIteration iterations],[(posteriors.sigma2.mean - confidence * posteriors.sigma2.standardDeviation) (posteriors.sigma2.mean - confidence * posteriors.sigma2.standardDeviation)],'--r','linewidth',2);
                plot([startIteration iterations],[(posteriors.sigma2.mean + confidence * posteriors.sigma2.standardDeviation) (posteriors.sigma2.mean + confidence * posteriors.sigma2.standardDeviation)],'--r','linewidth',2);
                hold off;
                ylim([minS*0.9 maxS*1.1]);
                title('sigma^2');
                %set(gca,'OuterPosition',[(mainPlotWidth + 1 / cols) 0.5 (1 / cols) 0.5]);
            elseif strcmpi(params.regressionMethod,'weightedleastsquares')
                [hAx,hLine1,hLine2] = plotyy(1:iterations,sigma2History,1:iterations,samples.regressionWeightsPhiHistory);
                hLine1.LineStyle = '-';
                hLine2.LineStyle = '-';
                title('sigma^2 and phi');
            end
            
            %degradation rate
            subplot(rows,cols,cols*2);
            plot(degRateHistory);
            hold on;
            maxD = max(degRateHistory((iterations * burnIn+1):iterations));
            minD = min(degRateHistory((iterations * burnIn+1):iterations));
            plot([iterations * burnIn iterations * burnIn],[ minD maxD], '--k','linewidth',2);
            %median and percentiles
            plot([(iterations * burnIn)+1 iterations],[ posteriors.degradationRate_mRNA.median posteriors.degradationRate_mRNA.median], 'k','linewidth',2);
            plot([(iterations * burnIn)+1 iterations],[ posteriors.degradationRate_mRNA.lowerQuartile posteriors.degradationRate_mRNA.lowerQuartile], '--k','linewidth',2);
            plot([(iterations * burnIn)+1 iterations],[ posteriors.degradationRate_mRNA.upperQuartile posteriors.degradationRate_mRNA.upperQuartile], '--k','linewidth',2);
            %mean and standard deviation
            plot([(iterations * burnIn)+1 iterations],[posteriors.degradationRate_mRNA.mean posteriors.degradationRate_mRNA.mean],'r','linewidth',2);
            plot([(iterations * burnIn)+1 iterations],[(posteriors.degradationRate_mRNA.mean - confidence * posteriors.degradationRate_mRNA.standardDeviation) (posteriors.degradationRate_mRNA.mean - confidence * posteriors.degradationRate_mRNA.standardDeviation)],'--r','linewidth',2);
            plot([(iterations * burnIn)+1 iterations],[(posteriors.degradationRate_mRNA.mean + confidence * posteriors.degradationRate_mRNA.standardDeviation) (posteriors.degradationRate_mRNA.mean + confidence * posteriors.degradationRate_mRNA.standardDeviation)],'--r','linewidth',2);
            hold off;
            title('Degradation rate');
            ylim([minD*0.9 maxD*1.1]);
            %set(gca,'OuterPosition',[mainPlotWidth 0 (1 / cols) 0.5]);
            
            %model size frequencies
            subplot(rows,cols,cols*3);
            firstNonZero = find(posteriors.switchFrequency.frequencies>0,1,'first');
            lastNonZero = find(posteriors.switchFrequency.frequencies>0,1,'last');
            r = lastNonZero - firstNonZero;
            bar(posteriors.switchFrequency.frequencies(firstNonZero:lastNonZero));
            if r < 7
                set(gca,'xtick',1:(r+1),'xticklabel',(firstNonZero-1):(lastNonZero-1));
            elseif r < 14
                set(gca,'xtick',1:2:(r+1),'xticklabel',(firstNonZero-1):2:(lastNonZero-1));
            else
                set(gca,'xtick',1:3:(r+1),'xticklabel',(firstNonZero-1):3:(lastNonZero-1));
            end
            xlim([0 r+2]);
            title('Model size frequency');
            xlabel('Number of switches');
        else
            %switch samples (thinned)
            
            samplesToPlot = round(linspace(1,iterations,params.maxSamplesToPlot));
            startIterationIdx = find(startIteration < samplesToPlot,1);
            
            mainArea = (cols*2+1):(cols*2+mainPlotCols);
            subplot(rows,cols,mainArea);
            
            swData = cell(maxSwitches,1);
            for x = 1:maxSwitches
                swData{x} = zeros(maxSamplesToPlot,1) * NaN;
            end
            
            count = 1;
            for x = samplesToPlot
                %plot( ones(length(swTimesHistory{x}),1) * x,swTimesHistory{x},'o','MarkerSize',1);
                currSwData = swTimesHistory{x};
                for y = 1:length(currSwData)
                    swData{y}(count) = currSwData(y);
                end
                count = count + 1;
            end
            for x = 1:maxSwitches
                plot(swData{x},'ob','MarkerSize',1);
                hold on;
            end
            
            %NOTE: we plot everything the 'wrong way round' and then rotate for efficient speed of plotting
            plot([startIterationIdx startIterationIdx],[timescale(1) timescale(end)],'--k','linewidth',2);
            if iterations > 1
                xlim([1 maxSamplesToPlot]);
            end
            
            %plot X% for each switch (+/- confidence SDs)
            for x = 1:length(swTimes)
                if taus(x+1) > taus(x)
                    swCol = 'r';
                else
                    swCol = 'g';
                end
                plot([startIterationIdx maxSamplesToPlot],[swTimes(x) swTimes(x)],swCol,'linewidth',2);
                plot([startIterationIdx maxSamplesToPlot],[swTimes(x) - confidence*swSigmas(x) swTimes(x) - confidence*swSigmas(x)],['--' swCol],'linewidth',2);
                plot([startIterationIdx maxSamplesToPlot],[swTimes(x) + confidence*swSigmas(x) swTimes(x) + confidence*swSigmas(x)],['--' swCol],'linewidth',2);
            end
            view(-90,90);
            set(gca,'ydir','reverse');
            ylim([timescale(1) timescale(end)]);
            ylabel('Timescale (hours)');
            xlabel('Switch samples (thinned)');
            hold off;
            
            %             subplot(rows,cols,cols);
            %             %just plot likelihood history
            %             plot(lLikelihoodHistory(samplesToPlot));
            %             hold on;
            %             maxL = max(lLikelihoodHistory(startIteration:iterations));
            %             minL = min(lLikelihoodHistory(startIteration:iterations));
            %             plot([startIterationIdx startIterationIdx],[minL maxL],'--k','linewidth',2);
            %             %median + percentiles
            %             plot([startIterationIdx maxSamplesToPlot],[lLikelihood.median lLikelihood.median],'k','linewidth',2);
            %             plot([startIterationIdx maxSamplesToPlot],[lLikelihood.lowerQuartile lLikelihood.lowerQuartile],'--k','linewidth',2);
            %             plot([startIterationIdx maxSamplesToPlot],[lLikelihood.upperQuartile lLikelihood.upperQuartile],'--k','linewidth',2);
            %             %mean + standard deviation
            %             plot([startIterationIdx maxSamplesToPlot],[lLikelihood.mean lLikelihood.mean],'r','linewidth',2);
            %             plot([startIterationIdx maxSamplesToPlot],[(lLikelihood.mean - confidence * lLikelihood.standardDeviation) (lLikelihood.mean - confidence * lLikelihood.standardDeviation)],'--r','linewidth',2);
            %             plot([startIterationIdx maxSamplesToPlot],[(lLikelihood.mean + confidence * lLikelihood.standardDeviation) (lLikelihood.mean + confidence * lLikelihood.standardDeviation)],'--r','linewidth',2);
            %             hold off;
            %             title('Log Likelihood (thinned)');
            %             %set the y limits based on burn in period
            %             minLMult = 0.9;
            %             if minL < 0
            %                 minLMult = 1.1;
            %             end
            %             maxLMult = 1.1;
            %             if maxL < 0
            %                 maxLMult = 0.9;
            %             end
            %             ylim([minL*minLMult maxL*maxLMult]);
            %             %set(gca,'OuterPosition',[mainPlotWidth 0.5 (1 / cols) 0.5]);
            
            %sigma history
            subplot(rows,cols,cols);
            if strcmpi(params.regressionMethod,'leastsquares')
                plot(sigma2History(samplesToPlot));
                hold on;
                maxS = max(sigma2History(startIteration:iterations));
                minS = min(sigma2History(startIteration:iterations));
                plot([startIterationIdx startIterationIdx],[minS maxS],'--k','linewidth',2);
                %median + percentiles
                plot([startIterationIdx maxSamplesToPlot],[posteriors.sigma2.median posteriors.sigma2.median],'k','linewidth',2);
                plot([startIterationIdx maxSamplesToPlot],[posteriors.sigma2.lowerQuartile posteriors.sigma2.lowerQuartile],'--k','linewidth',2);
                plot([startIterationIdx maxSamplesToPlot],[posteriors.sigma2.upperQuartile posteriors.sigma2.upperQuartile],'--k','linewidth',2);
                %mean + standard deviation
                plot([startIterationIdx maxSamplesToPlot],[posteriors.sigma2.mean posteriors.sigma2.mean],'r','linewidth',2);
                plot([startIterationIdx maxSamplesToPlot],[(posteriors.sigma2.mean - confidence * posteriors.sigma2.standardDeviation) (posteriors.sigma2.mean - confidence * posteriors.sigma2.standardDeviation)],'--r','linewidth',2);
                plot([startIterationIdx maxSamplesToPlot],[(posteriors.sigma2.mean + confidence * posteriors.sigma2.standardDeviation) (posteriors.sigma2.mean + confidence * posteriors.sigma2.standardDeviation)],'--r','linewidth',2);
                hold off;
                ylim([minS*0.9 maxS*1.1]);
                
                title('sigma^2 (thinned)');
            elseif strcmpi(params.regressionMethod,'weightedleastsquares')
                [hAx,hLine1,hLine2] = plotyy(1:length(samplesToPlot),sigma2History(samplesToPlot),1:length(samplesToPlot),samples.regressionWeightsPhiHistory(samplesToPlot));
                hLine1.LineStyle = '-';
                hLine2.LineStyle = '-';
                title('sigma^2 and phi (thinned)');
            end
            %set(gca,'OuterPosition',[(mainPlotWidth + 1 / cols) 0.5 (1 / cols) 0.5]);
            
            %degradation rate
            subplot(rows,cols,cols*2);
            plot(degRateHistory(samplesToPlot));
            hold on;
            maxD = max(degRateHistory((iterations * burnIn+1):iterations));
            minD = min(degRateHistory((iterations * burnIn+1):iterations));
            plot([startIterationIdx startIterationIdx],[ minD maxD], '--k','linewidth',2);
            %median and percentiles
            plot([startIterationIdx maxSamplesToPlot],[ posteriors.degradationRate_mRNA.median posteriors.degradationRate_mRNA.median], 'k','linewidth',2);
            plot([startIterationIdx maxSamplesToPlot],[ posteriors.degradationRate_mRNA.lowerQuartile posteriors.degradationRate_mRNA.lowerQuartile], '--k','linewidth',2);
            plot([startIterationIdx maxSamplesToPlot],[ posteriors.degradationRate_mRNA.upperQuartile posteriors.degradationRate_mRNA.upperQuartile], '--k','linewidth',2);
            %mean and standard deviation
            plot([startIterationIdx maxSamplesToPlot],[posteriors.degradationRate_mRNA.mean posteriors.degradationRate_mRNA.mean],'r','linewidth',2);
            plot([startIterationIdx maxSamplesToPlot],[(posteriors.degradationRate_mRNA.mean - confidence * posteriors.degradationRate_mRNA.standardDeviation) (posteriors.degradationRate_mRNA.mean - confidence * posteriors.degradationRate_mRNA.standardDeviation)],'--r','linewidth',2);
            plot([startIterationIdx maxSamplesToPlot],[(posteriors.degradationRate_mRNA.mean + confidence * posteriors.degradationRate_mRNA.standardDeviation) (posteriors.degradationRate_mRNA.mean + confidence * posteriors.degradationRate_mRNA.standardDeviation)],'--r','linewidth',2);
            hold off;
            title('Degradation rate (thinned)');
            ylim([minD*0.9 maxD*1.1]);
            
            %model size frequencies
            subplot(rows,cols,cols*3);
            firstNonZero = find(posteriors.switchFrequency.frequencies>0,1,'first');
            lastNonZero = find(posteriors.switchFrequency.frequencies>0,1,'last');
            r = lastNonZero - firstNonZero;
            bar(posteriors.switchFrequency.frequencies(firstNonZero:lastNonZero));
            if r < 7
                set(gca,'xtick',1:(r+1),'xticklabel',(firstNonZero-1):(lastNonZero-1));
            elseif r < 14
                set(gca,'xtick',1:2:(r+1),'xticklabel',(firstNonZero-1):2:(lastNonZero-1));
            else
                set(gca,'xtick',1:3:(r+1),'xticklabel',(firstNonZero-1):3:(lastNonZero-1));
            end
            xlim([0 r+2]);
            title('Model size frequency');
            xlabel('Number of switches');
        end
    else
        xlabel('Timescale (hours)');
    end
elseif strfind(modelType,'protein')
    if fullFit
        samples = varargin{5};
        swTimesHistory = samples.swTimesHistory;
        %     lLikelihood = fit.posteriors.logLikelihood;
        %     lLikelihoodHistory = fit.samples.logLikelihoodHistory;
        degRate_mRNAHistory = samples.degradationRate_mRNAHistory;
        degRate_proteinHistory = samples.degradationRate_proteinHistory;
        alphaHistory = samples.alphaHistory;
        %     M_0History = fit.samples.M_0History;
        %     P_0History = fit.samples.M_0History;
        sigma2History = samples.sigma2History;
        mainPlotCols = 3;
        cols = mainPlotCols + 1;
        rows = 4;
    else
        mainPlotCols = 1;
        cols = 1;
        rows = 3;
    end
    %adjustedTimescale = sort([timescale swTimes']);
    %adjustedProfile = nStateSwitchODE(adjustedTimescale,initialExp,birthRates,degRate,swTimes);
    
    %scrsz = get(0,'ScreenSize');
    %figure('Name',['ReTrOS-switch: ' name],'NumberTitle','off','Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)*0.5 scrsz(4)*0.5],'WindowStyle','docked');
    figure(figHandle);
    clf(figHandle);
    set(figHandle,'Name',['ReTrOS-switch: ' name],'NumberTitle','off');
    
    
    
    mainArea = 1:mainPlotCols;
    subplot(rows,cols,mainArea);
    dataCol = 'b';
    fillCol = [0.65 0.65 1];
    %fill function plots a single filled polygon, such that x1 is connected to xn (last point in vector)
    %therefore the last half of the x vector is reversed to do the top side of the polygon, and upper bound of data is reversed
    
    %check to see if we ignore all data for single timepoint
    hasData = lowerData>=0;
    if replicates > 1
        patch = fill( [timescale(hasData) fliplr(timescale(hasData))], [lowerData(hasData) fliplr(upperData(hasData))], fillCol); %final parameter is colour
        set(patch,'EdgeColor','none');%,'FaceAlpha',0.5); % makes the polygon transparent
    end
    hold on;
    %plot median profile
    plot(timescale,medianProfile,dataCol,'linewidth',2);
    for x = 1:replicates
        for y = 1:timepoints
            plot(timescale(y),data(x,y),[dataCol 'o'],'MarkerSize',4);
        end
    end
    
    % %plot the ignored data
    % for x = 1:replicates
    %     for y = 1:timepoints
    %         plot(timescale(y),maskedData(x,y),'ko','MarkerSize',4);
    %     end
    % end
    
    title(name,'Interpreter','none');
    
    
    % %plot model profile
    
    Xmatrix = cell(length(swTimes)+1,1);
    for x = 1:(length(swTimes)+1)
        Xmatrix{x} = zeros(timepoints * replicates,x+2);
    end
    
    [P_0, M_0, taus] = calculateModel_protein(timescale-timescale(1),swTimes-timescale(1),posteriors.degradationRate_mRNA.median,posteriors.degradationRate_protein.median,posteriors.alpha.median,Y,length(timescale),replicates,useDatapoint,Xmatrix,regression);
    
    timeSw = unique([linspace(timescale(1),timescale(end),timeResolution) swTimes' timescale(end)]);
    [mRNA, protein] = nStateSwitchODE_protein(timeSw-timescale(1),M_0,taus,posteriors.degradationRate_mRNA.median,swTimes-timescale(1),P_0,posteriors.alpha.median,posteriors.degradationRate_protein.median);
    plot(timeSw,protein,'k','linewidth',2);
    
    %plot switches
    for x = 1:length(swTimes)
        if taus(x+1) > taus(x)
            swCol = 'r';
        else
            swCol = 'g';
        end
        plot([swTimes(x) swTimes(x) ], [minData maxData],swCol,'LineWidth',2);
        plot([swTimes(x) - confidence*swSigmas(x) swTimes(x) - confidence*swSigmas(x) ], [minData maxData],['--' swCol],'linewidth',2);
        plot([swTimes(x) + confidence*swSigmas(x) swTimes(x) + confidence*swSigmas(x) ], [minData maxData],['--' swCol],'linewidth',2)
    end
    xlim([timescale(1) timescale(end)]);
    ylim([minData*0.9 maxData*1.1]);
    xlabel('Timescale in hours');
    ylabel('Protein expression');
    hold off;
    
    % set(gca,'OuterPosition',[0 0.5 mainPlotWidth 0.5]);
    
    %mRNA expression
    minmRNA = max(min(mRNA),0);
    maxmRNA = max(max(mRNA),0);
    
    mainArea = (cols+1):(cols+mainPlotCols);
    subplot(rows,cols,mainArea);
    plot(timeSw,mRNA,'k','linewidth',2);
    hold on;
    %plot switches
    for x = 1:length(swTimes)
        if taus(x+1) > taus(x)
            swCol = 'r';
        else
            swCol = 'g';
        end
        plot([swTimes(x) swTimes(x) ], [minmRNA maxmRNA],swCol,'LineWidth',2);
        plot([swTimes(x) - confidence*swSigmas(x) swTimes(x) - confidence*swSigmas(x) ], [minmRNA maxmRNA],['--' swCol],'linewidth',2);
        plot([swTimes(x) + confidence*swSigmas(x) swTimes(x) + confidence*swSigmas(x) ], [minmRNA maxmRNA],['--' swCol],'linewidth',2)
    end
    xlim([timescale(1) timescale(end)]);
    ylim([minmRNA * 0.9 (maxmRNA * 1.1)+eps]);
    xlabel('Timescale in hours');
    ylabel('mRNA expression');
    hold off;
    
    %switch distribution
    mainArea = (cols*2+1):(cols*2+mainPlotCols);
    subplot(rows,cols,mainArea);
    %plot the switch fits (gaussian)
    switchNum = length(swTimes);
    
    sFit = normFunc(swEstimationPoints,swTimes,swSigmas,swWeights);
    
    for z = 1:switchNum
        %colour the portion of the gaussian used to calculate switch size/strength
        l = swTimes(z) - confidence * swSigmas(z);
        u = swTimes(z) + confidence * swSigmas(z);
        idx = (swEstimationPoints >= l) .* (swEstimationPoints <= u);
        x1 = swEstimationPoints(idx == 1);
        x = [x1 fliplr(x1)];
        y1 = sFit(idx == 1)';
        y = [y1 zeros(1,length(x1))];
        
        if taus(z+1) > taus(z)
            swCol = [1 0.5 0.5];
        else
            swCol = [0.5 1 0.5];
        end
        p = fill(x, y,swCol);
        set(p,'edgecolor','none');%,'facealpha',0.5);
        hold on;
    end
    
    plot(swEstimationPoints,sFit,'k','linewidth',2);
    plot(swEstimationPoints,swFunction,'b','linewidth',2);
    hold off;
    ylabel('Switch PDF');
    xlim([timescale(1) timescale(end)]);
    
    % set(gca,'OuterPosition',[0 0 mainPlotWidth 0.5]);
    if fullFit
        if iterations <= maxSamplesToPlot
            %swiches samples
            mainPlotArea = (cols*3+1):(cols*3+mainPlotCols);
            subplot(rows,cols,mainPlotArea);
            hold on;
            swData = cell(maxSwitches,1);
            for x = 1:maxSwitches
                swData{x} = zeros(iterations,1) * NaN;
            end
            
            for x = 1:iterations
                %plot( ones(length(swTimesHistory{x}),1) * x,swTimesHistory{x},'o','MarkerSize',1);
                currSwData = swTimesHistory{x};
                for y = 1:length(currSwData)
                    swData{y}(x) = currSwData(y);
                end
            end
            for x = 1:maxSwitches
                plot(swData{x},'ob','MarkerSize',1);
                %    %plot(swTimesHistory(:,x),'o','MarkerSize',1);
            end
            plot([startIteration-1 startIteration-1],[timescale(1) timescale(end)],'--k','linewidth',2');
            if iterations > 1
                xlim([1 iterations]);
            end
            
            %plot X% for each switch (+/- confidence SDs)
            for x = 1:length(swTimes)
                if taus(x+1) > taus(x)
                    swCol = 'r';
                else
                    swCol = 'g';
                end
                plot([startIteration iterations], [swTimes(x) swTimes(x)], swCol,'linewidth',2);
                plot([startIteration iterations], [swTimes(x) - confidence*swSigmas(x) swTimes(x) - confidence*swSigmas(x)], ['--' swCol],'linewidth',2);
                plot([startIteration iterations], [swTimes(x) + confidence*swSigmas(x) swTimes(x) + confidence*swSigmas(x)], ['--' swCol],'linewidth',2);
            end
            view(-90,90);
            set(gca,'ydir','reverse');
            ylim([timescale(1) timescale(end)]);
            xlabel('Switch samples');
            ylabel('Timescale (hours)');
            %set(gca,'OuterPosition',[(mainPlotWidth + 2 / cols) 0 (2 / cols) 0.5]);
            %     pos = get(gca,'OuterPosition');
            %     pos(3) = 0.21;%65;
            %     set(gca,'OuterPosition',pos);
            hold off;
            
            %             %alpha history
            %             subplot(rows,cols,mainPlotCols + 1);
            %             plot(alphaHistory);
            %             hold on;
            %             maxA = max(alphaHistory(startIteration:iterations));
            %             minA = min(alphaHistory(startIteration:iterations));
            %             plot([startIteration-1 startIteration-1],[minA maxA],'--k','linewidth',2);
            %             %median + percentiles
            %             plot([startIteration iterations],[posteriors.alpha.median posteriors.alpha.median],'k','linewidth',2);
            %             plot([startIteration iterations],[posteriors.alpha.lowerQuartile posteriors.alpha.lowerQuartile],'--k','linewidth',2);
            %             plot([startIteration iterations],[posteriors.alpha.upperQuartile posteriors.alpha.upperQuartile],'--k','linewidth',2);
            %             %mean + standard deviation
            %             plot([startIteration iterations],[posteriors.alpha.mean posteriors.alpha.mean],'r','linewidth',2);
            %             plot([startIteration iterations],[(posteriors.alpha.mean - confidence * posteriors.alpha.standardDeviation) (posteriors.alpha.mean - confidence * posteriors.alpha.standardDeviation)],'--r','linewidth',2);
            %             plot([startIteration iterations],[(posteriors.alpha.mean + confidence * posteriors.alpha.standardDeviation) (posteriors.alpha.mean + confidence * posteriors.alpha.standardDeviation)],'--r','linewidth',2);
            %             hold off;
            %             ylim([minA*0.9 maxA*1.1]);
            %             title('alpha');
            
            
            
            %sigma history
            subplot(rows,cols,mainPlotCols + 1);
            if strcmpi(params.regressionMethod,'leastsquares')
                plot(sigma2History);
                hold on;
                maxS = max(sigma2History(startIteration:iterations));
                minS = min(sigma2History(startIteration:iterations));
                plot([startIteration-1 startIteration-1],[minS maxS],'--k','linewidth',2);
                %median + percentiles
                plot([startIteration iterations],[posteriors.sigma2.median posteriors.sigma2.median],'k','linewidth',2);
                plot([startIteration iterations],[posteriors.sigma2.lowerQuartile posteriors.sigma2.lowerQuartile],'--k','linewidth',2);
                plot([startIteration iterations],[posteriors.sigma2.upperQuartile posteriors.sigma2.upperQuartile],'--k','linewidth',2);
                %mean + standard deviation
                plot([startIteration iterations],[posteriors.sigma2.mean posteriors.sigma2.mean],'r','linewidth',2);
                plot([startIteration iterations],[(posteriors.sigma2.mean - confidence * posteriors.sigma2.standardDeviation) (posteriors.sigma2.mean - confidence * posteriors.sigma2.standardDeviation)],'--r','linewidth',2);
                plot([startIteration iterations],[(posteriors.sigma2.mean + confidence * posteriors.sigma2.standardDeviation) (posteriors.sigma2.mean + confidence * posteriors.sigma2.standardDeviation)],'--r','linewidth',2);
                hold off;
                ylim([minS*0.9 maxS*1.1]);
                title('sigma^2');
                %set(gca,'OuterPosition',[(mainPlotWidth + 1 / cols) 0.5 (1 / cols) 0.5]);
            elseif strcmpi(params.regressionMethod,'weightedleastsquares')
                [hAx,hLine1,hLine2] = plotyy(1:iterations,sigma2History,1:iterations,samples.regressionWeightsPhiHistory);
                hLine1.LineStyle = '-';
                hLine2.LineStyle = '-';
                title('sigma^2 and phi');
            end
            
            
            %degradation rates
            subplot(rows,cols,cols + mainPlotCols + 1);
            plot(degRate_mRNAHistory);
            hold on;
            maxD = max(degRate_mRNAHistory((iterations * burnIn+1):iterations));
            minD = min(degRate_mRNAHistory((iterations * burnIn+1):iterations));
            plot([iterations * burnIn iterations * burnIn],[ minD maxD], '--k','linewidth',2);
            %median and percentiles
            plot([(iterations * burnIn)+1 iterations],[ posteriors.degradationRate_mRNA.median posteriors.degradationRate_mRNA.median], 'k','linewidth',2);
            plot([(iterations * burnIn)+1 iterations],[ posteriors.degradationRate_mRNA.lowerQuartile posteriors.degradationRate_mRNA.lowerQuartile], '--k','linewidth',2);
            plot([(iterations * burnIn)+1 iterations],[ posteriors.degradationRate_mRNA.upperQuartile posteriors.degradationRate_mRNA.upperQuartile], '--k','linewidth',2);
            %mean and standard deviation
            plot([(iterations * burnIn)+1 iterations],[posteriors.degradationRate_mRNA.mean posteriors.degradationRate_mRNA.mean],'r','linewidth',2);
            plot([(iterations * burnIn)+1 iterations],[(posteriors.degradationRate_mRNA.mean - confidence * posteriors.degradationRate_mRNA.standardDeviation) (posteriors.degradationRate_mRNA.mean - confidence * posteriors.degradationRate_mRNA.standardDeviation)],'--r','linewidth',2);
            plot([(iterations * burnIn)+1 iterations],[(posteriors.degradationRate_mRNA.mean + confidence * posteriors.degradationRate_mRNA.standardDeviation) (posteriors.degradationRate_mRNA.mean + confidence * posteriors.degradationRate_mRNA.standardDeviation)],'--r','linewidth',2);
            hold off;
            title('mRNA degradation rate');
            ylim([minD*0.9 maxD*1.1]);
            %set(gca,'OuterPosition',[mainPlotWidth 0 (1 / cols) 0.5]);
            
            subplot(rows,cols,cols*2 + mainPlotCols+1);
            plot(degRate_proteinHistory);
            hold on;
            maxD = max(degRate_proteinHistory((iterations * burnIn+1):iterations));
            minD = min(degRate_proteinHistory((iterations * burnIn+1):iterations));
            plot([iterations * burnIn iterations * burnIn],[ minD maxD], '--k','linewidth',2);
            %median and percentiles
            plot([(iterations * burnIn)+1 iterations],[ posteriors.degradationRate_protein.median posteriors.degradationRate_protein.median], 'k','linewidth',2);
            plot([(iterations * burnIn)+1 iterations],[ posteriors.degradationRate_protein.lowerQuartile posteriors.degradationRate_protein.lowerQuartile], '--k','linewidth',2);
            plot([(iterations * burnIn)+1 iterations],[ posteriors.degradationRate_protein.upperQuartile posteriors.degradationRate_protein.upperQuartile], '--k','linewidth',2);
            %mean and standard deviation
            plot([(iterations * burnIn)+1 iterations],[posteriors.degradationRate_protein.mean posteriors.degradationRate_protein.mean],'r','linewidth',2);
            plot([(iterations * burnIn)+1 iterations],[(posteriors.degradationRate_protein.mean - confidence * posteriors.degradationRate_protein.standardDeviation) (posteriors.degradationRate_protein.mean - confidence * posteriors.degradationRate_protein.standardDeviation)],'--r','linewidth',2);
            plot([(iterations * burnIn)+1 iterations],[(posteriors.degradationRate_protein.mean + confidence * posteriors.degradationRate_protein.standardDeviation) (posteriors.degradationRate_protein.mean + confidence * posteriors.degradationRate_protein.standardDeviation)],'--r','linewidth',2);
            hold off;
            title('protein degradation rate');
            ylim([minD*0.9 maxD*1.1]);
            
            %model size frequencies
            subplot(rows,cols,cols*3 + mainPlotCols+1);
            firstNonZero = find(posteriors.switchFrequency.frequencies>0,1,'first');
            lastNonZero = find(posteriors.switchFrequency.frequencies>0,1,'last');
            r = lastNonZero - firstNonZero;
            bar(posteriors.switchFrequency.frequencies(firstNonZero:lastNonZero));
            if r < 7
                set(gca,'xtick',1:(r+1),'xticklabel',(firstNonZero-1):(lastNonZero-1));
            elseif r < 14
                set(gca,'xtick',1:2:(r+1),'xticklabel',(firstNonZero-1):2:(lastNonZero-1));
            else
                set(gca,'xtick',1:3:(r+1),'xticklabel',(firstNonZero-1):3:(lastNonZero-1));
            end
            xlim([0 r+2]);
            title('Model size frequency');
            xlabel('Number of switches');
        else
            %thin the samples
            samplesToPlot = round(linspace(1,iterations,params.maxSamplesToPlot));
            startIterationIdx = find(startIteration < samplesToPlot,1);
            
            %swiches samples
            mainPlotArea = (cols*3+1):(cols*3+mainPlotCols);
            subplot(rows,cols,mainPlotArea);
            hold on;
            swData = cell(maxSwitches,1);
            for x = 1:maxSwitches
                swData{x} = zeros(maxSamplesToPlot,1) * NaN;
            end
            
            count = 1;
            for x = samplesToPlot
                %plot( ones(length(swTimesHistory{x}),1) * x,swTimesHistory{x},'o','MarkerSize',1);
                currSwData = swTimesHistory{x};
                for y = 1:length(currSwData)
                    swData{y}(count) = currSwData(y);
                end
                count = count + 1;
            end
            for x = 1:maxSwitches
                plot(swData{x},'ob','MarkerSize',1);
                %    %plot(swTimesHistory(:,x),'o','MarkerSize',1);
            end
            
            plot([startIterationIdx startIterationIdx],[timescale(1) timescale(end)],'--k','linewidth',2');
            if iterations > 1
                xlim([1 maxSamplesToPlot]);
            end
            
            %plot X% for each switch (+/- confidence SDs)
            for x = 1:length(swTimes)
                if taus(x+1) > taus(x)
                    swCol = 'r';
                else
                    swCol = 'g';
                end
                plot([startIterationIdx maxSamplesToPlot], [swTimes(x) swTimes(x)], swCol,'linewidth',2);
                plot([startIterationIdx maxSamplesToPlot], [swTimes(x) - confidence*swSigmas(x) swTimes(x) - confidence*swSigmas(x)], ['--' swCol],'linewidth',2);
                plot([startIterationIdx maxSamplesToPlot], [swTimes(x) + confidence*swSigmas(x) swTimes(x) + confidence*swSigmas(x)], ['--' swCol],'linewidth',2);
            end
            view(-90,90);
            set(gca,'ydir','reverse');
            ylim([timescale(1) timescale(end)]);
            xlabel('Switch samples (thinned)');
            ylabel('Timescale (hours)');
            %set(gca,'OuterPosition',[(mainPlotWidth + 2 / cols) 0 (2 / cols) 0.5]);
            %     pos = get(gca,'OuterPosition');
            %     pos(3) = 0.21;%65;
            %     set(gca,'OuterPosition',pos);
            hold off;
            
            %             %alpha history
            %             subplot(rows,cols,mainPlotCols + 1);
            %             plot(alphaHistory(samplesToPlot));
            %             hold on;
            %             maxA = max(alphaHistory(startIteration:iterations));
            %             minA = min(alphaHistory(startIteration:iterations));
            %             plot([startIterationIdx startIterationIdx],[minA maxA],'--k','linewidth',2);
            %             %median + percentiles
            %             plot([startIterationIdx maxSamplesToPlot],[posteriors.alpha.median posteriors.alpha.median],'k','linewidth',2);
            %             plot([startIterationIdx maxSamplesToPlot],[posteriors.alpha.lowerQuartile posteriors.alpha.lowerQuartile],'--k','linewidth',2);
            %             plot([startIterationIdx maxSamplesToPlot],[posteriors.alpha.upperQuartile posteriors.alpha.upperQuartile],'--k','linewidth',2);
            %             %mean + standard deviation
            %             plot([startIterationIdx maxSamplesToPlot],[posteriors.alpha.mean posteriors.alpha.mean],'r','linewidth',2);
            %             plot([startIterationIdx maxSamplesToPlot],[(posteriors.alpha.mean - confidence * posteriors.alpha.standardDeviation) (posteriors.alpha.mean - confidence * posteriors.alpha.standardDeviation)],'--r','linewidth',2);
            %             plot([startIterationIdx maxSamplesToPlot],[(posteriors.alpha.mean + confidence * posteriors.alpha.standardDeviation) (posteriors.alpha.mean + confidence * posteriors.alpha.standardDeviation)],'--r','linewidth',2);
            %             hold off;
            %             ylim([minA*0.9 maxA*1.1]);
            %             title('alpha (thinned)');
                        
            %sigma history
            subplot(rows,cols,mainPlotCols + 1);
            if strcmpi(params.regressionMethod,'leastsquares')
                plot(sigma2History(samplesToPlot));
                hold on;
                maxS = max(sigma2History(startIteration:iterations));
                minS = min(sigma2History(startIteration:iterations));
                plot([startIterationIdx startIterationIdx],[minS maxS],'--k','linewidth',2);
                %median + percentiles
                plot([startIterationIdx maxSamplesToPlot],[posteriors.sigma2.median posteriors.sigma2.median],'k','linewidth',2);
                plot([startIterationIdx maxSamplesToPlot],[posteriors.sigma2.lowerQuartile posteriors.sigma2.lowerQuartile],'--k','linewidth',2);
                plot([startIterationIdx maxSamplesToPlot],[posteriors.sigma2.upperQuartile posteriors.sigma2.upperQuartile],'--k','linewidth',2);
                %mean + standard deviation
                plot([startIterationIdx maxSamplesToPlot],[posteriors.sigma2.mean posteriors.sigma2.mean],'r','linewidth',2);
                plot([startIterationIdx maxSamplesToPlot],[(posteriors.sigma2.mean - confidence * posteriors.sigma2.standardDeviation) (posteriors.sigma2.mean - confidence * posteriors.sigma2.standardDeviation)],'--r','linewidth',2);
                plot([startIterationIdx maxSamplesToPlot],[(posteriors.sigma2.mean + confidence * posteriors.sigma2.standardDeviation) (posteriors.sigma2.mean + confidence * posteriors.sigma2.standardDeviation)],'--r','linewidth',2);
                hold off;
                ylim([minS*0.9 maxS*1.1]);
                
                title('sigma^2 (thinned)');
            elseif strcmpi(params.regressionMethod,'weightedleastsquares')
                [hAx,hLine1,hLine2] = plotyy(1:length(samplesToPlot),sigma2History(samplesToPlot),1:length(samplesToPlot),samples.regressionWeightsPhiHistory(samplesToPlot));
                hLine1.LineStyle = '-';
                hLine2.LineStyle = '-';
                title('sigma^2 and phi (thinned)');
            end
            
            %degradation rates
            subplot(rows,cols,cols + mainPlotCols + 1);
            plot(degRate_mRNAHistory(samplesToPlot));
            hold on;
            maxD = max(degRate_mRNAHistory((iterations * burnIn+1):iterations));
            minD = min(degRate_mRNAHistory((iterations * burnIn+1):iterations));
            plot([startIterationIdx startIterationIdx],[ minD maxD], '--k','linewidth',2);
            %median and percentiles
            plot([startIterationIdx maxSamplesToPlot],[ posteriors.degradationRate_mRNA.median posteriors.degradationRate_mRNA.median], 'k','linewidth',2);
            plot([startIterationIdx maxSamplesToPlot],[ posteriors.degradationRate_mRNA.lowerQuartile posteriors.degradationRate_mRNA.lowerQuartile], '--k','linewidth',2);
            plot([startIterationIdx maxSamplesToPlot],[ posteriors.degradationRate_mRNA.upperQuartile posteriors.degradationRate_mRNA.upperQuartile], '--k','linewidth',2);
            %mean and standard deviation
            plot([startIterationIdx maxSamplesToPlot],[posteriors.degradationRate_mRNA.mean posteriors.degradationRate_mRNA.mean],'r','linewidth',2);
            plot([startIterationIdx maxSamplesToPlot],[(posteriors.degradationRate_mRNA.mean - confidence * posteriors.degradationRate_mRNA.standardDeviation) (posteriors.degradationRate_mRNA.mean - confidence * posteriors.degradationRate_mRNA.standardDeviation)],'--r','linewidth',2);
            plot([startIterationIdx maxSamplesToPlot],[(posteriors.degradationRate_mRNA.mean + confidence * posteriors.degradationRate_mRNA.standardDeviation) (posteriors.degradationRate_mRNA.mean + confidence * posteriors.degradationRate_mRNA.standardDeviation)],'--r','linewidth',2);
            hold off;
            title('mRNA degradation rate (thinned)');
            ylim([minD*0.9 maxD*1.1]);
            %set(gca,'OuterPosition',[mainPlotWidth 0 (1 / cols) 0.5]);
            
            subplot(rows,cols,cols*2 + mainPlotCols+1);
            plot(degRate_proteinHistory(samplesToPlot));
            hold on;
            maxD = max(degRate_proteinHistory((iterations * burnIn+1):iterations));
            minD = min(degRate_proteinHistory((iterations * burnIn+1):iterations));
            plot([startIterationIdx startIterationIdx],[ minD maxD], '--k','linewidth',2);
            %median and percentiles
            plot([startIterationIdx maxSamplesToPlot],[ posteriors.degradationRate_protein.median posteriors.degradationRate_protein.median], 'k','linewidth',2);
            plot([startIterationIdx maxSamplesToPlot],[ posteriors.degradationRate_protein.lowerQuartile posteriors.degradationRate_protein.lowerQuartile], '--k','linewidth',2);
            plot([startIterationIdx maxSamplesToPlot],[ posteriors.degradationRate_protein.upperQuartile posteriors.degradationRate_protein.upperQuartile], '--k','linewidth',2);
            %mean and standard deviation
            plot([startIterationIdx maxSamplesToPlot],[posteriors.degradationRate_protein.mean posteriors.degradationRate_protein.mean],'r','linewidth',2);
            plot([startIterationIdx maxSamplesToPlot],[(posteriors.degradationRate_protein.mean - confidence * posteriors.degradationRate_protein.standardDeviation) (posteriors.degradationRate_protein.mean - confidence * posteriors.degradationRate_protein.standardDeviation)],'--r','linewidth',2);
            plot([startIterationIdx maxSamplesToPlot],[(posteriors.degradationRate_protein.mean + confidence * posteriors.degradationRate_protein.standardDeviation) (posteriors.degradationRate_protein.mean + confidence * posteriors.degradationRate_protein.standardDeviation)],'--r','linewidth',2);
            hold off;
            title('protein degradation rate (thinned)');
            ylim([minD*0.9 maxD*1.1]);
            
            %model size frequencies
            subplot(rows,cols,cols*3 + mainPlotCols+1);
            firstNonZero = find(posteriors.switchFrequency.frequencies>0,1,'first');
            lastNonZero = find(posteriors.switchFrequency.frequencies>0,1,'last');
            r = lastNonZero - firstNonZero;
            bar(posteriors.switchFrequency.frequencies(firstNonZero:lastNonZero));
            if r < 7
                set(gca,'xtick',1:(r+1),'xticklabel',(firstNonZero-1):(lastNonZero-1));
            elseif r < 14
                set(gca,'xtick',1:2:(r+1),'xticklabel',(firstNonZero-1):2:(lastNonZero-1));
            else
                set(gca,'xtick',1:3:(r+1),'xticklabel',(firstNonZero-1):3:(lastNonZero-1));
            end
            xlim([0 r+2]);
            title('Model size frequency');
            xlabel('Number of switches');
        end
    else
        xlabel('Timescale (hours)');
    end
else
    disp('Cannot plot: Unknown model type');
end

end