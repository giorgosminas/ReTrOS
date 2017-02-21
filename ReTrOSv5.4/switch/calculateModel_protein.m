function [P_0,M_0,taus] = calculateModel_protein(timescale,sTimes,delta_m,delta_p,alpha,Y,timepoints,replicates,useSamples,Xmatrix,regression)

X = nStateSwitchModel_protein(timescale,sTimes,delta_m,delta_p,alpha);
X1 = Xmatrix{length(sTimes)+1};
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

%unconstrained least squares
beta = (X' * X) \ X' * Y;

if sum(beta(1:3)<0) > 0
    %negative P(0), M(0) or tau_0
    %reformulate model setting the coeffs to 0
    numCoeffs = size(X,2);
    idx = [find(beta(1:3)>=0)' 4:numCoeffs];
    
    X1 = X(:,idx);
    beta1 = (X1' * X1) \ X1' * Y;
    beta = zeros(numCoeffs,1);
    beta(idx) = beta1;
end

%constrained (non-negative) least squares
% alphas = lsqnonneg(X,Y);

P_0 = beta(1);
M_0 = beta(2);
taus = beta(3:end) .* delta_m;
for m = 2:length(taus)
    taus(m) = taus(m) + taus(m-1);
end
%taus = max(0,taus);

end
