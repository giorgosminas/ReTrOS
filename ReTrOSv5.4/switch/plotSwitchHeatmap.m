function [swData] = plotSwitchHeatmap(varargin)

timescale = [];
if isempty(varargin)
    %display file browser for user to select file
    
    %load switch file
    [FileName,PathName] = uigetfile('*.txt');
    
    if FileName == 0 %user cancelled
        disp('An input file must be selected');
        return;
    end
    
    contents = importdata(fullfile(PathName,FileName),'\t');
    fileName = fullfile(PathName,FileName);
else
    %use the input argument as file
    try
        contents = importdata(varargin{1},'\t');
        fileName = varargin{1};
        if length(varargin) > 1
            %may have timescale
            if isnumeric(varargin{2}) && length(varargin{2}) >= 2
                timescale = varargin{2};
            end
        end
    catch exception
        swData = [];
        disp(['Cannot load file "' varargin{1} '"']);
        return;
    end
end
if ~isstruct(contents)
    swData = [];
    disp('Invalid switch time file');
    return;
end
info = contents.textdata;
data = contents.data;

%we should have 5 columns: gene, switch type, switch mean time, switch sigma, switch strength
if size(info,2) ~= 5
    disp('incorrect number of columns in file - 5 columns required');
    return;
end

allGenes = info(2:end,1);
swType = info(2:end,2);
swTimes = data(:,1);
swSigmas = data(:,2);
swStrength = data(:,3);

%genes = unique(allGenes);
genes = uniqueRetainOrder(allGenes);
NoG = length(genes);

swData = cell(NoG,1);

for x = 1:NoG
    idx = strcmpi(genes{x},allGenes);
    currSw = struct;
    currSw.gene = genes{x};
    currSw.switchTimes = swTimes(idx);
    currSw.switchSigmas = swSigmas(idx);
    currSw.switchStrengths = swStrength(idx);
    currSw.switchType = swType(idx);
    swData{x} = currSw;
end
clf;
% set(gcf,'name','ReTrOS-switch heatmap','numbertitle','off');

ylabels = cell(NoG,1);
for x = 1:NoG
    currSw = swData{x};
    for y = 1:length(currSw.switchTimes)
        strength = currSw.switchStrengths(y);
        if strcmpi(currSw.switchType{y},'increase')
            edgeCol = [1 0 0];
            col = [1 1-strength 1-strength];
        else
            edgeCol = [0 1 0];
            col = [1-strength 1 1-strength];
        end
        mu = currSw.switchTimes(y);
        sigma = currSw.switchSigmas(y);
        p = fill([mu - 1.96*sigma mu + 1.96 * sigma mu + 1.96 * sigma mu - 1.96 * sigma],[x-0.45 x-0.45 x+0.45 x+0.45],col);
        set(p,'EdgeColor',edgeCol);
        hold on;
        plot([mu mu],[x-0.45 x+0.45],'color',edgeCol);
    end
    ylabels{x} = currSw.gene;
end
hold off;
[path name] = fileparts(fileName);
title(name,'interpreter','none');
xlabel('Timescale (hours)');
ylabel('Genes');
set(gca,'ytick',1:NoG,'yticklabel',ylabels,'ydir','reverse');
ylim([0.5 NoG+0.5]);
if ~isempty(timescale)
    xlim([timescale(1) timescale(end)]);
end

end

function u = uniqueRetainOrder(i)

u = {i{1}};

for j = 2:length(i)
    if sum(strcmpi(u,i{j})) == 0
        u = [u ; i{j}];
    end
end

end