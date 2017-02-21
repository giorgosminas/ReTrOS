function [] = plotSmooth(fit,data,figHandle)

%want the same style of plotting as from ReTrOS-switch
%protein/luc fit (optional)
%mRNA backcalc
%transcription backcalc
%native mRNA backcalc (optional)

if isempty(figHandle)
    figHandle = figure('name','RetrOS-smooth','numbertitle','off');
end

%do we have mRNA or protein data as input?
isProteinInput = true;
%if fit.algorithmParameters.reporterDegRates(3)==0 && fit.algorithmParameters.reporterDegRates(4) == 0
if strfind(fit.algorithmParameters.expressionType,'mRNA');
    %mRNA
    isProteinInput = false;
end

%do we have native mRNA degradation rates?
hasNatDeg = false;
if ~isempty(fit.algorithmParameters.nativeDegradationRate_mRNA)
    hasNatDeg = true;
end
timescale = fit.algorithmParameters.time;
timescaleUnique = unique(timescale);
maxReps = 0;
for x = 1:length(timescaleUnique)
    if sum(timescaleUnique(x) == timescale) > maxReps
        maxReps = sum(timescaleUnique(x) == timescale);
    end
end

dataReps = nan(maxReps,length(timescaleUnique));
for x = 1:length(timescaleUnique)
    d = data(timescaleUnique(x) == timescale);
    for y = 1:length(d)
        dataReps(y,x) = d(y);
    end
end

time = timescale(1):fit.algorithmParameters.timeResolution:timescale(end);
%now what plots do we need to make?

%we have at least 2/3 time-series plots, and possibly an extra
if isProteinInput
    rows = 3;
else
    rows = 2;
end
cols = 1;
if hasNatDeg
    rows = rows+1;
end

figure(figHandle);
%protein expression
subplot(rows,cols,1);
dataCol = 'b';
%initial data input
for x = 1:maxReps
    plot(timescaleUnique,dataReps(x,:),dataCol,'linewidth',2);
    hold on;
end
plot(timescale,data,[dataCol 'o'],'MarkerSize',6);

%smoothed protein expression
plot(timescaleUnique,fit.smoothFit.smoothedReporter.fit(1:length(timescaleUnique)),'-g','linewidth',2);

%protein expression from backcalculated mRNA
plot(time,fit.smoothFit.reporter(:,2),'-k','linewidth',2);
plot(time,fit.smoothFit.reporter(:,1),'--k','linewidth',2);
plot(time,fit.smoothFit.reporter(:,3),'--k','linewidth',2);
hold off;
xlim([timescale(1) timescale(end)]);
title(fit.profileName,'Interpreter','none');
if isProteinInput
    ylabel('Protein expression');
else
    ylabel('mRNA expression');
end
currY = ylim;
if currY(1) < 0
    currY(1) = 0;
    ylim(currY);
end

currPlot = 2;
if isProteinInput
    %backcalculated mRNA expression
    subplot(rows,cols,currPlot);
    plot(time(1:(end-1)),fit.smoothFit.reporter_mRNA(:,2),'-r','linewidth',2);
    hold on;
    plot(time(1:(end-1)),fit.smoothFit.reporter_mRNA(:,1),'--r','linewidth',2);
    plot(time(1:(end-1)),fit.smoothFit.reporter_mRNA(:,3),'--r','linewidth',2);
    hold off;
    xlim([timescale(1) timescale(end)]);
    ylabel('mRNA expression');
    currY = ylim;
    if currY(1) < 0
        currY(1) = 0;
        ylim(currY);
    end
    currPlot = currPlot+1;
end

%backcalculated transcription rates
subplot(rows,cols,currPlot);
plot(time(1:(end-2)),fit.smoothFit.tau(:,2),'-b','linewidth',2);
hold on;
plot(time(1:(end-2)),fit.smoothFit.tau(:,1),'--b','linewidth',2);
plot(time(1:(end-2)),fit.smoothFit.tau(:,3),'--b','linewidth',2);
hold off;
xlim([timescale(1) timescale(end)]);
ylabel('Transcription rate');
currY = ylim;
if currY(1) < 0
    currY(1) = 0;
    ylim(currY);
end

currPlot = currPlot + 1;

if hasNatDeg
    subplot(rows,cols,currPlot);
    plot(timescaleUnique,fit.smoothFit.native_mRNA(:,2),'-m','linewidth',2);
    hold on;
    plot(timescaleUnique,fit.smoothFit.native_mRNA(:,1),'--m','linewidth',2);
    plot(timescaleUnique,fit.smoothFit.native_mRNA(:,3),'--m','linewidth',2);
    hold off;
    xlim([timescale(1) timescale(end)]);
    xlabel('Timescale (hours)');
    ylabel('Native mRNA expression');
    currY = ylim;
    if currY(1) < 0
        currY(1) = 0;
        ylim(currY);
    end
else
    xlabel('Timescale (hours)');
end

end

% subplot(4,2,[1 3]); set(gca,'FontSize',12,'LineWidth',2);
% plot(time,luc,time,r.f); legend('data','model'); title(label);
%
% ylabel(['Luc ' detrend_str]);
% xlabel('Time (hours)');
% subplot(4,2,[5 7]); set(gca,'FontSize',12,'LineWidth',2);
% plot(h,mse,h(jhmin),mse(jhmin),'rs');
% legend('sum sq. error','bottom'); xlabel('band width (hr)');
% % plot(time,err,'k'); legend('error'); xlabel('Time (hours)');
%
% if ~isempty(deg_nat)
%
%     subplot(4,2,2); set(gca,'FontSize',12,'LineWidth',2);
%     plot(time,mrnanat_qt(:,2),'-m',time,mrnanat_qt(:,1),'--m',time,mrnanat_qt(:,3),'--m');
%     ylabel('nat. mRNA');
%
%     subplot(4,2,4); set(gca,'FontSize',12,'LineWidth',2);
%     plot(time_int(1:end-2),trans_qt(:,2),'b',time_int(1:end-2),trans_qt(:,1),'--b',time_int(1:end-2),trans_qt(:,3),'--b');
%     ylabel('Transcription');
%
%     subplot(4,2,6); set(gca,'FontSize',12,'LineWidth',2);
%     plot(time_int(1:end-1),mrnaluc_qt(:,2),'r',time_int(1:end-1),mrnaluc_qt(:,1),'--r',time_int(1:end-1),mrnaluc_qt(:,3),'--r');
%     ylabel('LUC mRNA');
%
%     subplot(4,2,8); set(gca,'FontSize',12,'LineWidth',2);
%     plot(time_int,luc_qt(:,2),'k',time_int,luc_qt(:,1),'--k',time_int,luc_qt(:,3),'--k');
%     hold on;
%     plot(time,luc,'-ok','MarkerSize',5);
%     hold off;
%     ylabel('LUC'); xlabel('Time (hours)');
%
% else
%
%     subplot(4,2,[2 4]); set(gca,'FontSize',12,'LineWidth',2);
%     plot(time_int(1:end-2),trans_qt(:,2),'b',time_int(1:end-2),trans_qt(:,1),'--b',time_int(1:end-2),trans_qt(:,3),'--b');
%     ylabel('Transcription');
%
%     subplot(4,2,6); set(gca,'FontSize',12,'LineWidth',2);
%     plot(time_int(1:end-1),mrnaluc_qt(:,2),'r',time_int(1:end-1),mrnaluc_qt(:,1),'--r',time_int(1:end-1),mrnaluc_qt(:,3),'--r');
%     ylabel('LUC mRNA');
%
%     subplot(4,2,8); set(gca,'FontSize',12,'LineWidth',2);
% %     plot(time_int,luc_qt(:,2),'k',time_int,luc_qt(:,1),'--k',time_int,luc_qt(:,3),'--k');
%     plot(time,luc_orig,'-b');%,'MarkerSize',5);
%     hold on;
%     if ~isempty(trend)
%         plot(time,trend,'--k');
%     end
%     hold off;
%     ylabel('LUC'); xlabel('Time (hours)');
%
% end