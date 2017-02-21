function llik = logLikelihood(e_hat,prec,regression)
%function llik = logLikelihood(data,profile,e_hat,prec)
% [nReps timepoints] = size(data);
% e_hat = repmat(profile,nReps,1) - data;
% e_hatT = e_hat';
% e_hat = e_hatT(:);

switch regression.type
    case 'leastsquares'
        llik = logLikelihoodLS(e_hat,prec);
    case 'weightedleastsquares',
        llik = logLikelihoodWLS(e_hat,prec,regression.Q);
end
end

function llik = logLikelihoodLS(e_hat,prec)
sig2 = 1 / prec;
llikvec = log(sqrt(2*pi*sig2)) + ((e_hat).^2)./(2*sig2);
llik = -1 * sum(llikvec);

end


function llik = logLikelihoodWLS(e_hat,prec,Q)
sig2 = 1 / prec ;
Sig = sig2*Q ;
n = length(e_hat) ;
llik = - n*log(2*pi)/2 - n*log(sig2)/2 - sum(log(diag(Q)))/2 - (e_hat'*(Sig\e_hat))/2 ;
end