% Function to simulate mRNA and transcription profiles

function [mrna_nat,mrna,transcription,luc_int] =...
    backcalculationKR(params,delta,time,time_int,luc,bandWidth,...
    degM,degP,deg_nat,var_nat)
%global gLUCDEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT
% delta         time step for Euler method
% time          vector of observed time points (including 0)
% time_int      time points for fitted LUC spline
% luc
% bandWidth
% degM          degradation rate for LUC mRNA
% degP          degradation rate for LUC
% deg_nat       ...
% var_nat       ...
%
% OUTPUT
% mrna_nat      approximate native mRNA
% mrna:         approximate LUC mRNA
% transcription approximate transcription function
% luc_int       fitted LUC spline evaluated at 0.1h intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% model a smooth curve
r = ksrlin(params,time,luc,bandWidth,time);
err = luc - r.f; er2 = err .^ 2;
w = ksrlin(params,time,er2,bandWidth,time);
sigma = sqrt(w.f); % plot(time,er2,time,sigma.^2); pause;

%%%%%%%  SIMULATION w/ SIGN-CONSERVING NORMAL DISTRIB [1]
noise = sigma .* randn(size(time,1),1);
noise = noise .* sign(noise) .* sign(err);              % noise def.
luc_sim = r.f + noise;                                  % simulation
r = ksrlin(params,time,luc_sim,bandWidth,time_int);            % re-modelling
luc_int = r.f;
%%%%%%%

%if params.reporterDegRates(3)==0 && params.reporterDegRates(4) == 0
if strfind(params.expressionType,'mRNA')
    mrna = luc_int(1:end-1);                                % mRNA data
else
    mrna = diff(luc_int)/delta + degP * luc_int(1:end-1);   % protein data
end

transcription = diff(mrna)/delta + degM * mrna(1:end-1);

time = unique(time);

switch isempty(deg_nat) + isempty(var_nat)
    case 2                  % no native mRNA deg. rate provided
        mrna_nat = zeros(size(time));
        
    case 1                  % deg. rate of native mRNA at 'temp' known
        mrna_nat = trapeziumHM(transcription,deg_nat,delta,time,mrna(1));
        
    case 0                  % mean deg. rate & assoc. std error provided
        % sample a 'new' native mRNA deg. rate from a Gamma distrib.
        a_gam = (deg_nat^2)/var_nat;
        b_gam = var_nat/deg_nat;
        deg_nat_new = gamrnd(a_gam,b_gam);
        
        % approx. native mRNA
        mrna_nat = trapeziumHM(transcription,deg_nat_new,delta,time,mrna(1));
end

% [1] For SIMULATION w/ GAMMA DISTRIB see ReTrOSv1_20120917

end