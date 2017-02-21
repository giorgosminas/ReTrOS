function [initialExp birthRates] = calculateModel_mRNA(timescale,swTimes,degRate,Y,timepoints,replicates,useSamples,Xmatrix,regression)
X = nStateSwitchModel_mRNA(timescale,swTimes,degRate);
%X = repmat(X,replicates,1);
X1 = Xmatrix{length(swTimes)+1};
for z = 1:replicates
    X1( (timepoints * (z-1) + 1):(timepoints * z),:) = X;
end
X = X1(useSamples,:);
Y = Y(useSamples);
if strcmpi(regression.type,'weightedleastsquares')
    Qinvroot = regression.Q^(-1/2) ;
    X = Qinvroot*X; 
    Y = Qinvroot*Y;
end

%least squares regression to obtain model initial exp and birth rates
%alpha = inv(X' * X) * X' * Y;
beta = (X' * X) \ X' * Y;

if sum(beta(1:2)<0) > 0
    %negative M(0) or tau_0
    %reformulate model setting the coeffs to 0
    numCoeffs = size(X,2);
    idx = [find(beta(1:2)>=0)' 3:numCoeffs];
    
    X1 = X(:,idx);
    beta1 = (X1' * X1) \ X1' * Y;
    beta = zeros(numCoeffs,1);
    beta(idx) = beta1;
end

initialExp = beta(1);
birthRates = beta(2:length(beta)) * degRate;
for m = 2:length(birthRates)
    birthRates(m) = birthRates(m) + birthRates(m-1);
end

end