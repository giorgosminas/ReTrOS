% Function to apply the trapezoidal rule

function [rna_nat] = trapeziumHM(tau,deg,delta,time,m0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT %
% tau: reconstrcuted transcription function
% deg: native mRNA degradation rate
% delta: time step for Euler/trapezium rule
% time: time points for numerical integration
% m0: native mRNA at time point 0h (taken to be equal to LUC mRNA)
%
% OUTPUT %
% rna_nat: reconstructed native mRNA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time interval for numerical integration
n_time = length(time);

% preallocation for final rna output
rna_nat = zeros(size(time));
rna_nat(1) = exp(-deg*time(1))*m0;

for t = 2:(n_time-1)
    
    s = (time(1):delta:time(t))';
    %s = (0:delta:time(t))';
%     size(s)
%     size(tau)
%     size(tau(1:length(s)))
%     size(exp(-deg*(time(t)-s)))
    Y = tau(1:length(s)) .* exp(-deg*(time(t)-s));
    
    trapezoid = trapz(s,Y);
    rna_nat(t) = trapezoid;
    
end

% Last integration step
s = (time(1):delta:time(end)-2*delta)';
Y = tau(1:length(s)) .* exp(-deg*(s(end)-s));
trapezoid = trapz(s,Y);
rna_nat(n_time) = trapezoid;

rna_nat(2:end) = rna_nat(2:end) + rna_nat(1);

end