function [dM,dP] = drawDegRates(params)
%global gLUCDEG

dMa = params.reporterDegRates(1);
dMv = params.reporterDegRates(2) * params.reporterDegRates(2);
dPa = params.reporterDegRates(3);
dPv = params.reporterDegRates(4) * params.reporterDegRates(4);

% mRNA degradation rate (dM)
a = dMa^2 / dMv;
b = dMv / dMa;
dM = gamrnd(a,b);
%clear a b;

% protein degradation rate (dP)
a = dPa^2 / dPv;
b = dPv / dPa;
dP = gamrnd(a,b);

end