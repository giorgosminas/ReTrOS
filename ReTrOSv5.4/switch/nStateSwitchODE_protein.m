function [mRNA protein] = nStateSwitchODE_protein(time,M_0,taus,delta_m,sTimes,P_0,alpha,delta_p)

%compute mRNA and protein ODE profiles using the linear model: alpha_0 * x_0 + ... + alpha_n * x_n
%NOTE: we have renamed alpha_0 to alpha_n as theta_0 to theta_n

nSwitches = length(sTimes);

%model alpha_0 .. x_n (for both ODEs)
theta = zeros(1,3 + nSwitches);
theta(1) = P_0;             %initial protein
theta(2) = M_0;             %initial mRNA
theta(3) = taus(1) / delta_m; %birth rate 1

for x = 2:(nSwitches + 1)
    theta(x+2) = (taus(x) - taus(x-1)) / delta_m;   %change in birth rates
end

%mRNA x_0 .. x_n
xm_0 = exp(-delta_m.*time);  %proportion of initial state remaining
xm_1 = 1-xm_0;                %proportion of initial state degraded

%protein x_0 ... x_n
delta_pm_inv = 1 / ((delta_p - delta_m)+eps);
delta_p_inv = 1 / delta_p;

xp_0 = exp(-delta_p .* time);
xp_1 = alpha / ((delta_p - delta_m)+eps) .* (exp(-delta_m .* time) - exp(-delta_p .* time));
xp_2 = alpha * (delta_p_inv .* (1-exp(-delta_p .* time)) - delta_pm_inv .* (exp(-delta_m .* time) - exp(-delta_p .* time)));

%compute profiles
mRNA = theta(2) * xm_0 + theta(3) * xm_1;
protein = theta(1) * xp_0 + theta(2) * xp_1 + theta(3) * xp_2;

for n = 1:nSwitches
    ind = time > sTimes(n);
    t = time - sTimes(n);
    mRNA = mRNA + (theta(n+3) .* ind) .* (1-exp(-delta_m*t));
    protein = protein + (theta(n+3) .* ind) .* alpha .* (delta_p_inv .* (1-exp(-delta_p.*t)) - delta_pm_inv .* (exp(-delta_m .* t) -exp(-delta_p .* t)));
end

end


% function [mRNA protein] = nStateSwitchProteinODE(time,M_0,taus,delta_m,sTimes,P_0,alpha,delta_p,obsIdx)
% 
% nSwitches = length(sTimes);
% 
% theta = zeros(1,2 + nSwitches);
% theta(1) = M_0;            %initial state
% theta(2) = taus(1) / delta_m; %birth rate 1
% 
% for x = 2:(nSwitches + 1)
%     theta(x+1) = (taus(x) - taus(x-1)) / delta_m;   %change in birth rates
% end
% 
% x_0 = exp(-delta_m).^time;  %proportion of initial state remaining
% x_1 = 1-x_0;                %proportion of initial state degraded
% 
% mRNA = theta(1) * x_0 + theta(2) * x_1;
% 
% for n = 1:nSwitches
%     ind = time > sTimes(n);
%     mRNA = mRNA + (theta(n+2) .* ind) .* (1-exp(-delta_m*(time - sTimes(n))));
% end
% 
% protein = P_0 .* exp(-delta_p .* time);
% 
% for y = 2:length(time)
%     protein(y) = protein(y) + alpha * exp(-delta_p * time(y)) * fastTrapz( time(1:y), exp(delta_p .* time(1:y)) .* mRNA(1:y) );
%     %protein(y) = protein(y) + alpha * exp(-delta_p * time(y)) * trapz(time(1:y),exp(delta_p .* time(1:y)) .* mRNA(1:y));
%     %protein(y) = protein(y) + alpha * exp(-delta_p * time(y)) * quad(@(x)mRNAODE(x,M_0,taus,delta_m,sTimes,delta_p),time(1),time(y));
% end
% 
% mRNA = mRNA(obsIdx);
% protein = protein(obsIdx);
% 
% end
% 
% function z = fastTrapz(x,y)
% 
% z = sum(diff(x) .* (y(1:(end-1)) + y(2:end))/2);
% 
% end
% 
% function [x] = mRNAODE(time,M_0,taus,delta_m,sTimes,delta_p)
% 
% nSwitches = length(sTimes);
% 
% theta = zeros(1,2 + nSwitches);
% theta(1) = M_0;            %initial state
% theta(2) = taus(1) / delta_m; %birth rate 1
% 
% for x = 2:(nSwitches + 1)
%     theta(x+1) = (taus(x) - taus(x-1)) / delta_m;   %change in birth rates
% end
% 
% x_0 = exp(-delta_m).^time;  %proportion of initial state remaining
% x_1 = 1-x_0;                %proportion of initial state degraded
% 
% mRNA = theta(1) * x_0 + theta(2) * x_1;
% 
% for n = 1:nSwitches
%     ind = time > sTimes(n);
%     mRNA = mRNA + (theta(n+2) .* ind) .* (1-exp(-delta_m*(time - sTimes(n))));
% end
% 
% x = exp(delta_p .* time) .* mRNA;
% 
% end