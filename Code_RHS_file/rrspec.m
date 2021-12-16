function [rrtfr coef dcoef tcoef] = rrspec(x,fs,k,f);
% rrtfr = rrspec(x,fs,k,f) computes the recursively reassigned STFT spectrogram with 
% the constants defined
% x    - input signal
% fs   - sampling frequency
% k    - order
% f    - frequency vector
%
%
% Equations refered to are found in the journal paper "Recursive 
% Time-Frequency Reassignment" available in IEEE's Transactions on Signal 
% Processing. This script is found online at
% http://www.ii.uib.no/~geirkn/rrspec/
%
% Author: Geir Kjetil Nilsen (geir.kjetil.nilsen@gmail.com) 2009
% Updated December, 2012.
%

if k < 2 || k > 4
    sprintf('k must be in the range 2 <= k <= 4')
    rrtfr = 0;
    return;
end;

T = 1/fs;
nf = length(f);
N = length(x);
w = zeros(nf, N);
w_d = zeros(nf, N);
w_t = zeros(nf, N);
rrtfr = zeros(nf,N);

Omega = f(2)-f(1);

% Balance time-frequency resolution in number of samples, eq. 16
sigma_p = (sqrt(Omega)*factorial(k-1))/(sqrt(2*pi*T)*(k-1)^(k-1)*exp(-(k-1))); 

coef = cell(2,1);
dcoef = cell(2,1);
tcoef = cell(2,1);

for j = 1:nf;
    omega_p = 2*pi*f(j);
    p = -sigma_p + i*omega_p; % eq. 6
    a = exp(p*T);
     
    % Coefficients from eq. 14, k = 2...5

    % Numerator coefficients in coef{i}, i = 2n+1, n = 1,2,3,4
    % Demoninator coefficients in coef{i}, i = 2n, n = 2,3,4,5
    
    coef{3} = sigma_p^2*T^2*[0 a];
    coef{4} = [1 -2*a a.^2];
    coef{5} = sigma_p^3*T^3*[0 1/2*a 1/2*a.^2];
    coef{6} = [1 -3*a 3*a.^2 -a.^3];
    coef{7} = sigma_p^4*T^4*[0 1/6*a 2/3*a.^2 1/6*a.^3];
    coef{8} = [1 -4*a 6*a.^2 -4*a.^3 a.^4];
    coef{9} = sigma_p^5*T^5*[0 1/24*a 11/24*a.^2 11/24*a.^3 1/24*a.^4];
    coef{10}= [1 -5*a 10*a.^2 -10*a.^3 5*a.^4 -a.^5];
    
    % Coefficients from eq. 23, k = 2...4
    dcoef{3} = sigma_p.^2*T*[1 -a*(1 - p*T)];
    dcoef{4} = coef{4}; 
    dcoef{5} = sigma_p.^3*T.^2*[0 a/2*(2 + p*T) -a.^2/2*(2 - p*T)];
    dcoef{6} = coef{6}; 
    dcoef{7} = sigma_p.^4*T.^3*[0 a/6*(3 + p*T) a.^2*2/3*p*T -a.^3/6*(3 - p*T)];
    dcoef{8} = coef{8};
    
    % Coefficients from eq. 24, k = 2...4
    tcoef{3} = 2/sigma_p*coef{5}; 
    tcoef{4} = coef{6}; 
    tcoef{5} = 3/sigma_p*coef{7}; 
    tcoef{6} = coef{8}; 
    tcoef{7} = 4/sigma_p*coef{9}; 
    tcoef{8} = coef{10}; 
    
    
    w(j,:)   = filter( coef{2*k-1},  coef{2*k},x); % Eq. 16
    w_d(j,:) = filter(dcoef{2*k-1}, dcoef{2*k},x); % Freq. reass. coeffs with eq. 16
    w_t(j,:) = filter(tcoef{2*k-1}, tcoef{2*k},x); % Time reass. coeffs with eq. 16
end;


rtfr = abs(w);
step = f(2)-f(1); 
f0 = f(1);

% eq. 27
for n = 1:N;
    for f = 1:nf;
        that = n - round(1/T*real(w_t(f,n)./w(f,n))); % eq. 26
        that = min(max(that, 1),N); % edge
                      
        fhat = round(1/step*imag(w_d(f,n)./w(f,n))/2.0/pi - 1/step*(f0)); % eq. 25
        fhat = mod(fhat, nf) + 1; % edge
        
        rrtfr(nf - fhat + 1,that) = rrtfr(nf - fhat + 1,that) + rtfr(f,n); 
     end;
 end;
  