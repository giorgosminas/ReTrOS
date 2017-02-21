function [f x h baseline mus sigmas heights] = getSwitchFit(switchSamples,timescale,estimationPoints,switchDelay,baselineFactor,isGMM,isPlot)

%first get density estimate
points = estimationPoints;
%incSize = (timescale(end)-timescale(1)) / points;

%samples = [switchSamples{ startIteration:end }];
samples = [switchSamples{:}];

if isempty(samples)
    disp('No switch samples');
    f = zeros(points,1);
    x = zeros(points,1);
    h = [];
    baseline = zeros(points,1);
    mus = [];
    sigmas = [];
    heights = [];
    return;
end


minSwTimes = timescale(1) + switchDelay;
maxSwTimes = timescale(end) - switchDelay;
% minSwTimes = 3*(timescale(2) - timescale(1))/2;
% maxSwTimes = timescale(end) - (3*(timescale(2) - timescale(1))/2);
% minSwTimes = 3/2;
% maxSwTimes = timescale(end) - (3/2);



deltaTime = zeros(1,length(timescale)-1);
for x=1:length(deltaTime)
    deltaTime(x) = timescale(x+1) - timescale(x);
end
timePerTimepoint = min(deltaTime);

% % Should the density be weighted by the size of transcriptional jump
% % Possible not if the rates are somewhat specified by the RJ gemoetrical
% % mean criteria
% bData = cell(maxSwitches+1,1);
% for x = 1:(maxSwitches + 1)
%     bData{x} = zeros(usableIterations,1) * NaN;
% end
% for x = 1:usableIterations
%     currBData = rates{startIteration + x - 1};
%     for y = 1:length(currBData)
%         bData{y}(x) = currBData(y);
%     end
% end
% bmatrix = nan(usableIterations, maxSwitches+1);
% for x = 1:(maxSwitches + 1)
%     bmatrix(:,x) = bData{x};
% end
%     
%     
%  swData = cell(maxSwitches,1);
% for x = 1:maxSwitches
%     swData{x} = zeros(usableIterations,1) * NaN;
% end
% 
% for x = 1:usableIterations
%     currSwData = switchSamples{startIteration + x - 1};
%     for y = 1:length(currSwData)
%         swData{y}(x) = currSwData(y);
%     end
% end   
% swmatrix = nan(usableIterations, maxSwitches);
% for x = 1:(maxSwitches)
%     swmatrix(:,x) = swData{x};
% end
% 
% % bdiff = (abs(bmatrix(:,2:end) - bmatrix(:,1:(end-1)))).^2;
% % bdiff_vec = reshape(bdiff, numel(bdiff), 1);
% % sw_vec = reshape(swmatrix, 1, numel(swmatrix));
% % bdiff_vec = bdiff_vec(~isnan(bdiff_vec));
% % sw_vec = sw_vec(~isnan(sw_vec));
% % bweights = (bdiff_vec)/sum((bdiff_vec));
% % n = length(sw_vec);
% % [f x h] = ksdensity(sw_vec, linspace(timescale(1),timescale(end),points),'weights', bweights);

% or if not
%[f x h] = ksdensity(samples,linspace(timescale(1),timescale(end),points));%, 'support', [timescale(1) timescale(end)]);
%h = (points / (timescale(end)-timescale(1))) ^ -0.4;

%oversmoothing bandwidth for Gaussian kernel as given in Scott (1992), equation 6.66 and section 6.5.1.2
h_oversmooth = 1.14 * sqrt(var(samples)) * (length(samples) ^ (-1/5));

h = testBandwidth(h_oversmooth,samples,timescale);
[f x] = ksdensity(samples,linspace(timescale(1),timescale(end),points),'kernel', 'normal','width',h);%, 'support', [timescale(1) timescale(end)]);
n = length(samples);


funif = nan(1, length(x));
xunif = linspace(minSwTimes, maxSwTimes, n);
for i = 1:length(x)
    funif(i) = sum(normpdf((x(i) - xunif)/h))/(n*h);
end

% have a baseline percentage, of which the density must be above
baseline = max(funif .* baselineFactor,0);
% baseline = 1/(maxSwTimes - minSwTimes);
% [p idx] = findpeaks(f,'sortstr','descend');

origF = f;

%30.07.2014
%Actually remove baseline from density, rather than setting to 0 if below
%f = max(f .* (baseline <= f),0);
f = max(f - baseline,0);

auc_f = trapz(x,f);
% are we fitting a mixture model to more than 5% of the overall density
if auc_f < 0.05/2 
    disp('Mixture model cannot be fitted to density less than 2.5% of overall switch density')
    disp('Not fitting mixture model')
    f = zeros(points,1);
    x = zeros(points,1);
    baseline = zeros(points,1);
    mus = [];
    sigmas = [];
    heights = [];
    return;
end
[p idx] = getPeaks(f);
% disp(p)
%remove all peaks which are < X% of largest peak
%once we have a good estimate of the noise/baseline we WONT want to do this
%no need now we have the uniform + %age function for baseline
% bigP = find(p( p / p(1) >= 0) );
% % bigP = find(p( p >= baselineMax) );
% p = p(bigP);
% idx = idx(bigP);

bigIdx = idx;

% For flat distributions the baseine will be very high and associated peaks
% very small
% Therefore discard any peak that is less than 0.2*baseline, may want
% increasing - is there a general rule one could use
% baseP = find(p(p >= 0*baseline));

% p = p(baseP);
% idx = idx(baseP);
% bigIdx = idx;

    
mus = x(idx);
%we now need to check if any switches are too close
m = 1;

while m < length(p)
    extremes = getFullWidthHalfMaximum(idx(m),p(m),f);
    if extremes(1) == 0 || extremes(2) == 0
        overlaps = 1:p;
    else
        %%% this *can* cause index out of bounds
        % Attempted to access x(1001); index out of bounds because numel(x)=1000.
        withinEx = ( mus >= x(extremes(1)) ) .* ( mus <= x(extremes(2)));
        %%%
        withinTime = ( mus >= (mus(m) - switchDelay) ) .* ( mus <= (mus(m) + switchDelay));
        within = withinEx + withinTime;
        overlaps = find( within >= 1);
    end
    if length(overlaps) > 1
        %we want to remove those indicies which overlap with switch m (can only be AFTER index m)
        keep = setdiff(1:length(p), setdiff(overlaps,1:m) );
        p = p(keep);
        idx = idx(keep);
        mus = mus(keep);
    end
    m = m+1;
end

numPeaks = length(p);
idx = sort(idx,'ascend');

%just set the troughs to be halfway between subsequent switches
tIdx = zeros(numPeaks+1,1);
tIdx(1) = 1;
tIdx(end) = length(x);

for i = 2:numPeaks
    %we want to take the minimum value between the two switches, NOT the mid point index
%     offset = idx(i-1);
%     [minValue minIdx] = min( f(idx(i-1):idx(i)) );
%     tIdx(i) = minIdx + offset -1;

    %we will actually take the mid point, as if the pdf is non-0 then it is now part of a peak at this point
    tIdx(i) = round((idx(i) - idx(i-1)) / 2 + idx(i-1));
end


mus = zeros(numPeaks,1);
sigmas = zeros(numPeaks,1);
heights = zeros(numPeaks,1);

% %fit a gaussian to each mode individually
for z = 1:numPeaks
    try
        %how about we shift the baseline, as it currently tries to fit to 0 baseline
        swF = f(tIdx(z):tIdx(z+1));
        swX = x(tIdx(z):tIdx(z+1));
        %sFit = fit(x(tIdx(z):tIdx(z+1))',swF','gauss1');
        [m idx] = max(swF);
        sFit = fit(swX',swF',createGaussEquation(swX(idx)));
        mus(z) = sFit.b;
        sigmas(z) = sFit.c;
        heights(z) = sFit.a;
        
        if isPlot
            figure('windowstyle','docked');
            plot(swX,swF,'b');
            hold on;
            plot(swX,normFunc(swX,mus(z),sigmas(z),heights(z)),'-r');
            title(['Start point: ' num2str(swX(idx)) ' Mean: ' num2str(mus(z))]);
        end
    catch
        disp('Cannot create fit - individual gaussian');

        mus = [];
        sigmas = [];
        heights = [];
        return;
    end
end

%check for overlaps between gaussians
%minSw = mus - 1.96 * sigmas;
%maxSw = mus + 1.96 * sigmas;
minSw = mus - 1 * sigmas;
maxSw = mus + 1 * sigmas;
overlaps = zeros(numPeaks-1,1);

for c = 1:(numPeaks-1)
    overlaps(c) = maxSw(c) > minSw(c+1);
end

if sum(overlaps) > 0
    numPeaks = numPeaks - sum(overlaps);
    newT = zeros(numPeaks+1,1);

    newT(1) = 1;
    newT(end) = points;

    currT = 2;
    for c = 1:length(overlaps)
        if ~overlaps(c)
            newT(currT) = tIdx(c+1);
            currT = currT + 1;
        end
    end
    tIdx = newT;

    mus = zeros(numPeaks,1);
    sigmas = zeros(numPeaks,1);
    heights = zeros(numPeaks,1);
    for z = 1:numPeaks
        try
            %how about we shift the baseline, as it currently tries to fit to 0 baseline
            swF = f(tIdx(z):tIdx(z+1));
            swX = x(tIdx(z):tIdx(z+1));
            %sFit = fit(x(tIdx(z):tIdx(z+1))',swF','gauss1');
            [m idx] = max(swF);
            sFit = fit(swX',swF',createGaussEquation(swX(idx)));
            mus(z) = sFit.b;
            sigmas(z) = sFit.c;
            heights(z) = sFit.a;
        catch
            disp('Cannot create fit - individual gaussian');

            mus = [];
            sigmas = [];
            heights = [];
            return;
        end
    end
end


if isGMM
    disp('fitting a mixture model');
    %try fitting a mixture
    %however, we wont let the means be selected - we will use the existing means from the individual fits
    try
        indSigmas = sigmas;
        indHeights = heights;
        sFit = fit(x', f', createGaussFixedMeanEquation(mus,indSigmas,indHeights));
        for m = 1:numPeaks
            sigmas(m) = getfield(sFit,['c' int2str(m)]);
            heights(m) = getfield(sFit,['a' int2str(m)]);
        end
    catch
        disp('Cannot create fit - mixture of gaussians');
        mus = [];
        sigmas = [];
        heights = [];
        return;
    end
end

%sort the means
[sMus idx] = sort(mus);
sSigmas = zeros(numPeaks,1);
sHeights = zeros(numPeaks,1);

for m = 1:numPeaks
    sSigmas(m) = sigmas(idx(m));
    sHeights(m) = heights(idx(m));
end
mus = sMus;
sigmas = sSigmas;
heights = sHeights;


if isPlot
    figure('Name','Switch fitting','NumberTitle','off','WindowStyle','docked');
    subplot(3,1,1)
    plot(samples, 1:length(samples), 'b.', 'markersize', 1)
    ylim([0 length(samples)])
    xlim([timescale(1) timescale(end)]);

    subplot(3,1,2);
    plot(x,origF,'b','linewidth',2);
    hold on;
    plot(x,baseline,'g','linewidth',2);
    hold off;
    subplot(3,1,3);
%     plot(x,f);
    hold on;
    if isGMM
        plot(x,normFunc(x,mus,sigmas,heights),'r','linewidth',2);
    else
        for z = 1:numPeaks
            sF = heights(z) * exp( -((x-mus(z)).^2/(2 *sigmas(z).^2)));
            plot(x,sF,'r');
        end
    end
    plot(x(bigIdx),f(bigIdx),'og','linewidth',1.5);
    plot(x(idx),f(idx),'or','linewidth',1.5);
    plot(x(tIdx),f(tIdx),'ok','linewidth',1.5);
    hold off;
end

% disp(numPeaks)

end

function [ftype] = createGaussEquation(suggestMu)

%gauss1 model
eqn = 'a*exp(-((x-b)^2/(2*c^2)))';

foptions = fitoptions('Method','NonlinearLeastSquares','Lower',zeros(3,1),'StartPoint',[1 suggestMu 1]);
ftype = fittype(eqn,'options',foptions);

end

function [ftype] = createGaussFixedMeanEquation(mus,suggestSigma,suggestHeight)

numCoeffs = length(mus) * 2;

eqn = '';

coeffStr = cell(numCoeffs,1);
for x = 1:(length(mus)-1)
    eqn = [eqn 'a' int2str(x) '*exp(-((x-' num2str(mus(x)) ')^2/(2*c' int2str(x) '^2)))+'];
    coeffStr{x*2-1} = ['a' int2str(x)];
    coeffStr{x*2} = ['c' int2str(x)];
end

%last one
eqn = [eqn 'a' int2str(numCoeffs / 2) '*exp(-((x-' num2str(mus(end)) ')^2/(2*c' int2str(numCoeffs / 2) '^2)))'];
coeffStr{numCoeffs-1} = ['a' int2str(numCoeffs / 2)];
coeffStr{numCoeffs} = ['c' int2str(numCoeffs / 2)];

suggestStarts = nan(numCoeffs,1);
for x = 1:length(mus)
    suggestStarts(x*2-1) = suggestHeight(x);
    suggestStarts(x*2) = suggestSigma(x);
end
% suggestStarts = rand(numCoeffs,1);
foptions = fitoptions('Method','NonlinearLeastSquares','Lower',zeros(numCoeffs,1),'Startpoint',suggestStarts);
ftype = fittype(eqn,'options',foptions,'coefficients',coeffStr);
end

function [peaks loc] = getPeaks(f)

%we'll have a neighbourhood of 7, rather than 3
neighbourhood = 9;
nHalf = floor(neighbourhood / 2);
pks = zeros(length(f),1);
for x = (nHalf+1):(length(f)-nHalf)
    if max(f( (x-nHalf):(x+nHalf) ) ) == f(x) && f(x) > 0
        %this is a local peak
        pks(x) = 1;
    end
end

pks = find(pks);
peaks = f(pks);
loc = pks;

%sort by largest peak
[peaks idx] = sort(peaks,'descend');
nLoc = loc;
for x = 1:length(loc)
    nLoc(x) = loc(idx(x));
end
loc = nLoc;

end

function [extremes] = getFullWidthHalfMaximum(loc,peak, f)

halfPeak = peak / 2;
%we now need to find where the mid points lie
extremes = zeros(2,1);
%find mid-point to left
leftF = f(1:loc);
leftF = fliplr(leftF);
leftHalfPeak = find(leftF < halfPeak,1);
if ~isempty(leftHalfPeak)
    extremes(1) = loc - leftHalfPeak +1;
end

%find mid-point to right
rightF = f(loc:end);
rightHalfPeak = find(rightF < halfPeak,1);
if ~isempty(rightHalfPeak)
    extremes(2) = rightHalfPeak + loc-1;
end

end

function [hBest] = testBandwidth(h_oversmooth,samples,timescale)
%this function uses an implementation of unbiased cross-validation, as given in Scott (1992) equation 6.67

n = length(samples);
range = (timescale(end) - timescale(1)) * 1.01;
deltaTime = min(timescale(2:end)-timescale(1:(end-1)));
nbins = ceil(range*10/deltaTime); %1000;%ceil(range * 5 / deltaTime)
d = range / nbins;

bins = timescale(1):d:(d*(nbins-1)+timescale(1));
bins(1) = -inf;
bins(end) = inf;
counts = histc(samples,bins);
cnt = zeros(nbins,1);

%convert counts to distance
for i = 1:(nbins-1)
    for j = 0:(i-1)
        cnt(i-j+1) = cnt(i-j+1) + counts(i+1) * counts(j+1);
    end
    cnt(1) = cnt(1) + (counts(i+1) * (counts(i+1)-1))/2;
end

hTest = d:(d/10):h_oversmooth;

u = nan(length(hTest),1);

for i = 1:length(hTest)
    u(i) = bw_ucv(n,d,cnt,hTest(i));
end

[minU idx] = min(u);
hBest = hTest(idx);

end

function u = bw_ucv(n,d,counts,h)
%unbiased cross-validation function for bandwidth selection from: Multivariate Density Estimation: Theory, Practice, and Visualization (Scott 1992)) eqn 6.67
DELMAX = 1000;

nbin = length(counts);
delta = ((0:(nbin-1)) .* d ./ h) .^ 2;
term = exp(-delta / 4) - sqrt(8) * exp(-delta / 2);
delta = delta < DELMAX;
sumTerms = sum(term.*counts'.*delta);

u = 1 / (2 * n * h * sqrt(pi)) + sumTerms / (n*n*h * sqrt(pi));
end
