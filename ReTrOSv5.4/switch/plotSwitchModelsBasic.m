function [] = plotSwitchModelsBasic(posteriors,samples,parameters,profileInfo,figHandle)

if isempty(figHandle)
    figHandle = figure;
end

figure(figHandle);
clf(figHandle);
models = posteriors.switchModelsBasic.models;

swEstimationPoints = linspace(profileInfo.timescale(1),profileInfo.timescale(end),parameters.estimationPoints);
for x = 1:size(models,1)
    subplot(size(models,1),1,x);
    plot(swEstimationPoints,models{x,8},'-k','linewidth',2);
    hold on;
    proportion = models{x,5} / (parameters.iterations - parameters.startIteration + 1);
    plot(swEstimationPoints,normFunc(swEstimationPoints,models{x,2},models{x,3},models{x,4}),'-b','linewidth',2);
    title([profileInfo.profileName ': Switch model size: ' int2str(models{x,1}) ' - ' int2str(models{x,5}) ' samples (' num2str(proportion) ')']);
    if x == size(models,1)
        xlabel('Timescale (hours');
    end
    hold off;
end

end