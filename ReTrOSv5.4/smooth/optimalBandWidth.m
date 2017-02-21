function [h,mse,jhmin] = optimalBandWidth(params,time,yval)
%global gBWVECTOR gD2LCHECK

%%%%%%%  BAND WIDTH OPTIMISATION

NoC = length(time);

h = params.bandwidthRange(1) : params.bandwidthRange(3) : params.bandwidthRange(2);                            	% BW vector


mse = zeros(length(h),1);

for j = 1:length(h)                             % leave-one-out CV
    err = zeros(NoC,1);
    for i = 2:NoC-1 % = 4:NoC-3
        sime = [time(1:i-1); time(i+1:NoC)];        % time v. w/o time(i)
        xval = [yval(1:i-1); yval(i+1:NoC)];        % yval v. w/o yval(i)
        
        r = ksrlin(params,sime,xval,h(j),time(i));         % model at time(i)
        err(i) = yval(i) - r.f;                     % err estimate
    end
    mse(j) = (err' * err);
%     mse(j) = (err' * err) * (1+1*min(h(j)-1.5,0).^2);
%     fprintf('%f\t%f\n',h(j),mse(j));
end
%clear r err


% DEPRECATED FUNCTIONALITY
% if params.smoothnessCheck
%     smoothEnough = zeros(length(h),1);
%     for j = 1:length(h)
%         r = ksrlin(params,time,yval,h(j),time);
%         smoothEnough(j) =...
%             min( diff(diff(r.f)) + (dPdM*r.f(2:end-1)) ) >= 0;
%     end
%     h = h(smoothEnough > 0);
%     mse = mse(smoothEnough > 0);
% end

[dummy,jhmin] = min(mse);
BW = h(jhmin);              fprintf('\nopt. band width = %f\n',BW);

%%%%%%%

end