function [sF] = normFunc(x,mus,sigmas,weights)
sF = zeros(length(x),1);

numPeaks = length(mus);
for m = 1:numPeaks
    sF = sF + (weights(m) * exp( -((x-mus(m)).^2/(2 * sigmas(m).^2))))';
end

end