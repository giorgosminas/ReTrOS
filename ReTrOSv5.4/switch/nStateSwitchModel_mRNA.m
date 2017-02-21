function X = nStateSwitchModel_mRNA(time, sTimes, degRate)
%we dont need the birthrates/thetas, as these will be estimated via regression

nSwitches = length(sTimes);

X = zeros(length(time), nSwitches+2);

X(:,1) = exp(-degRate) .^ time;  %proportion of initial state remaining
X(:,2) = 1 - X(:,1);                %proportion of initial state degraded

for n = 1:nSwitches
    ind = time > sTimes(n);
    X(:,n+2) = (1 - exp(-degRate*(time - sTimes(n)))) .* ind;
end

end