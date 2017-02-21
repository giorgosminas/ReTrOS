function [X] = nStateSwitchModel_protein(time,sTimes,delta_m,delta_p,alpha)

nSwitches = length(sTimes);

delta = delta_p - delta_m + eps;
delta_pInv = 1 / delta_p;
delta_Inv = 1 / delta;

X = zeros(length(time), nSwitches+3);

X(:,1) = exp(-delta_p .* time);
X(:,2) = (alpha / delta) .* (exp(-delta_m .* time) - exp(-delta_p .* time));
X(:,3) = alpha .* ( delta_pInv .* (1 - exp(-delta_p .* time)) - delta_Inv .* (exp(-delta_m .* time) - exp(-delta_p .* time)) );

for n = 1:nSwitches
    ind = time > sTimes(n);
    t = time - sTimes(n);
    X(:,n+3) = (alpha .* ( delta_pInv .* (1 - exp(-delta_p .* t) - delta_Inv .* (exp(-delta_m .* t) - exp(-delta_p .* t))) ) ) .* ind;
end

end