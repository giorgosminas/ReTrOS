function profile = nStateSwitchODE_mRNA(time,initialState,birthRates,degRate,sTimes)

nSwitches = length(sTimes);

theta = zeros(1,2 + nSwitches);
theta(1) = initialState;            %initial state
theta(2) = birthRates(1) / degRate; %birth rate 1

for x = 2:(nSwitches + 1)
    theta(x+1) = (birthRates(x) - birthRates(x-1)) / degRate;   %change in birth rates
end

x_0 = exp(-degRate).^time;  %proportion of initial state remaining
x_1 = 1-x_0;                %proportion of initial state degraded

profile = theta(1) * x_0 + theta(2) * x_1;

for n = 1:nSwitches
    ind = time > sTimes(n);
    profile = profile + (theta(n+2) .* ind) .* (1-exp(-degRate*(time - sTimes(n))));
end

end