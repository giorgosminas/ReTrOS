function r=ksrlin(params,x,y,h,xp)
%global gKERNEL gPOLDEG

%% WITH P-TH POLYNOMIAL (DESIGN MATRIX)

p = params.polynomialDegree;

if isempty(params.kernel) || strcmpi(params.kernel,'Ga')
    kerf = @(z)exp(-z.*z/2)/sqrt(2*pi);             % Gaussian (default)
elseif  strcmpi(params.kernel,'Ep')
    kerf = @(z)(3/4)*(1-z.^2).*(abs(z)<1);          % Epanechinikov
elseif  strcmpi(params.kernel,'Tr')
    kerf = @(z)(35/32)*((1-z.^2).^3).*(abs(z)<1);   % Triweight
end

N = length(xp);

r.x = xp;
r.f = zeros(N,1);
for k = 1:N
    d = r.x(k) - x;
    z = kerf(d/h);
    if p == 1
        X = [ones(size(d)) d];              % p = 1
    elseif p == 2
        X = [ones(size(d)) d d.^2];         % p = 2
    else
        X = [ones(size(d)) d d.^2 d.^3];    % p = 3
    end
    W = diag(z);                            % weight (each time point)
    beta = (X' * W * X) \ (X' * W * y);     % regression (each time point)
    r.f(k) = beta(1);                       % estimate = intercept
end

%% WITH LIEAR MODEL (FORMULA)

% N = length(xp);
% 
% kerf = @(z)exp(-z.*z/2)/sqrt(2*pi);
% 
% r.x=xp;
% r.f=zeros(N,1);
% for k=1:N
%     d=r.x(k)-x;
%     z=kerf(d/h);
%     s1=d.*z;
%     s2=sum(d.*s1);
%     s1=sum(s1);
%     r.f(k)=sum((s2-s1*d).*z.*y)/(s2*sum(z)-s1*s1);
% end

end