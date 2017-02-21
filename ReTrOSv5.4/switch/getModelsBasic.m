function [models] = getModelsBasic(swFreq, samples, timescale,estimationPoints, switchDelay,switchBaselineFactor)

[sFreq sIdx] = sort(swFreq,'descend');
modelSize = zeros(length(samples),1);
useSwSamples = samples;
for x = 1:length(useSwSamples)
    modelSize(x) = length(useSwSamples{x});
end

minCoverage = 0.9;
minFreq = 0.1;
coverage = 0;
currModel = 1;

models = [];

while coverage < minCoverage && currModel <= length(swFreq) && swFreq(sIdx(currModel)) >= minFreq
    disp(['Dimension: ' int2str(sIdx(currModel)-1) ' - ' num2str(swFreq(sIdx(currModel)))]);
    
    swSamples = useSwSamples(modelSize == (sIdx(currModel)-1));
    
    if sIdx(currModel) > 1
        [f x bandwidth baseline mus sigmas heights] = getSwitchFit(swSamples,timescale,estimationPoints,switchDelay,switchBaselineFactor,true,false);
        
        if length(mus) > (sIdx(currModel)-1)
            %need to split models
            %length(mus)
            disp('Need to split models');
        end
    else
        mus = [];
        sigmas = [];
        bandwidth = [];
        f = [];
    end
    
    coverage = coverage + swFreq(sIdx(currModel));
    currModel = currModel + 1;
    models = [models ;  {sIdx(currModel)-1, mus, sigmas, heights, length(swSamples), bandwidth, swFreq(sIdx(currModel)), f}];
end
% clf;
% plotModels(models,parameters,profileInfo);

end
