% Function to compute confidence envelopes for a time series

function structure = ReTrOSsmooth(params,luc,time,deg_nat,var_nat,label,figHandle)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT
% luc:      luminescence traces
% time:     observed time points (in hours)
% deg_nat   (mean) deg. rate of native mRNA
% var_nat   variance of deg. rate of native mRNA
%
% OUTPUT
% structure: structure cell array with reconstructed profiles for LUC,
%               LUC mRNA, transcription and possibly also native mRNA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%boot_samples  = 99; % = 999;
boot_samples = params.bootstrapSamples;
% n             = length(luc);
delta         = params.timeResolution;
time_orig = time;
time          = time - time(1); % = (0:sampling:sampling*n-1)';
n = length(unique(time));
time_int      = (0:delta:time(end))';
n_int         = length(time_int);

luc_orig = luc;

detrend_str = '';
switch params.detrend
    case 'log-linear'
        maxReps = 0;
        t = unique(time);
        for x = 1:length(t)
            if sum(t(x) == time) > maxReps
                maxReps = sum(t(x) == time);
            end
        end
        detrendLuc = nan(1,length(luc));
        trend = nan(1,length(luc));
        for x = 1:maxReps
            [l tr] = detrendLinear(log(luc( ((x-1)*length(t)+1):(x*length(t)) )),t,1);
            detrendLuc( ((x-1)*length(t)+1):(x*length(t)) ) = l;
            trend( ((x-1)*length(t)+1):(x*length(t)) ) = tr;
        end
        luc = exp(detrendLuc)';
        trend = exp(trend)';
        detrend_str = '(detrend: log-linear)';
    case 'linear'
        maxReps = 0;
        t = unique(time);
        for x = 1:length(t)
            if sum(t(x) == time) > maxReps
                maxReps = sum(t(x) == time);
            end
        end
        detrendLuc = nan(1,length(luc));
        trend = nan(1,length(luc));
        for x = 1:maxReps
            [l tr] = detrendLinear(luc( ((x-1)*length(t)+1):(x*length(t)) ),t,1);
            detrendLuc( ((x-1)*length(t)+1):(x*length(t)) ) = l;
            trend( ((x-1)*length(t)+1):(x*length(t)) ) = tr;
        end
        luc = detrendLuc';
        trend = trend';
        detrend_str = '(detrend: linear)';
    otherwise
        trend = [];
        
end

% STORAGE MATRICES INITIALISATION
luc_boot      = zeros(n_int,   boot_samples);
mrnaluc_boot  = zeros(n_int-1, boot_samples);
trans_boot    = zeros(n_int-2, boot_samples);
mrna_nat_boot = zeros(n, boot_samples);

disp('Optimising bandwidth');
% FIND OPTIMAL BAND WIDTH FOR KERNEL REGRESSION

%MODIFICATION: 13.10.2014 - remove NaN entries from luc and time
luc_detrend = luc;
idx = luc >= 0;
luc = luc(idx);
time = time(idx);

[h,mse,jhmin] = optimalBandWidth(params,time,luc);
r = ksrlin(params,time,luc,h(jhmin),time);
% err = luc - r.f; sigma = abs(err);
bandwidthStruct = struct('fit',r.f,'bandwidthRange',h,'SSE',mse,'minBandwidth',h(jhmin));

disp('Running bootstrap');
for i = 1:boot_samples	% SIMS WITH BOOT-SAMPLED DEG. RATES

    [degM_boot,degP_boot] = drawDegRates(params);
    
    % find bootstrap mRNA and transcription traces
    [mrna_nat,mrna,transcription,luc_int] =...
        backcalculationKR(params,delta,time,time_int,luc,h(jhmin),...
        degM_boot,degP_boot,deg_nat,var_nat);
    
    mrnaluc_boot(:,i)  = mrna;
    trans_boot(:,i)    = transcription;
    luc_boot(:,i)      = luc_int;
    mrna_nat_boot(:,i) = mrna_nat;
    
    if rem(i,10)==9, fprintf('%4d: done\n',i); end
           
end

% CONFIDENCE ENVELOPES
luc_qt     = quantile(luc_boot,      [0.025 0.50 0.975], 2);
mrnaluc_qt = quantile(mrnaluc_boot,  [0.025 0.50 0.975], 2);
mrnanat_qt = quantile(mrna_nat_boot, [0.025 0.50 0.975], 2); 
trans_qt   = quantile(trans_boot,    [0.025 0.50 0.975], 2);

fitStructure = struct('reporter',luc_qt,'reporter_mRNA',mrnaluc_qt,'native_mRNA',mrnanat_qt,'tau',trans_qt,'smoothedReporter',bandwidthStruct);

params.nativeDegradationRate_mRNA = [deg_nat,var_nat];
params.time = time_orig;
r = RandStream.getGlobalStream;
r = RandStream.getGlobalStream;
params.randomNumberGenerator = r.Type;
params.randomNumberSeed = r.Seed;
params.randomStreamIndex = r.StreamIndex;
structure = struct('profileName',label,'smoothFit',fitStructure,'algorithmParameters',params);

plotSmooth(structure,luc_detrend,figHandle);

% % PLOTS
% %figure(1);
% clf;
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


end


function [y ax] = detrendMovingAvg(x)

sW13 = [1/24;repmat(1/12,11,1);1/24];
Ys = conv(Y,sW13,'same');
Ys(1:6) = Ys(7); Ys(N-5:N) = Ys(N-6);

xt = Y-Ys;

end

function [y ax] = detrendLinear(x,timescale,replicates)

%adapted from MATLAB 'detrend' function

t = timescale ./ timescale(end);

%no breakpoints, so bp would be [1 (length(obs)-1)], lp is therefore 1
%bp = unique([0;double(bp(:));N-1]);	% Include both endpoints
%lb = length(bp)-1;
% Build regressor with possibly non-linear pieces + DC
a  = [repmat(t,1,replicates) ones(length(t) * replicates,1)];
y = x - a*(a\x);		% Remove best fit
ax = a*(a\x);

end

% function [y ax] = detrendLinear(x)
% 
% %adapted from MATLAB 'detrend' function
% 
% N = length(x);
% 
% %no breakpoints, so bp would be [1 (length(obs)-1)], lp is therefore 1
% %bp = unique([0;double(bp(:));N-1]);	% Include both endpoints 
% %lb = length(bp)-1;
% 
% % Build regressor with linear pieces + DC
% a  = [(1:N)'/N ones(N,1)];
% y = x - a*(a\x);		% Remove best fit
% ax = a*(a\x);
% 
% end