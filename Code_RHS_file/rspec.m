function rtfr = rspec(x,fs,k,f)
% rtfr = rtfr(x,fs,k,f) computes the recursive STFT spectrogram with 
% the constants defined
% x    - input signal
% fs   - sampling frequency
% k    - order
% f    - frequency vector
%
% Equations refered to are found in the journal paper "Recursive 
% Time-Frequency Reassignment" available in IEEE's Transactions on Signal 
% Processing. This script is found online at
% http://www.ii.uib.no/~geirkn/rrspec/
%
% Author: Geir Kjetil Nilsen (geir.kjetil.nilsen@gmail.com) 2009
% Updated December, 2012.
%

if k < 2 || k > 5
    sprintf('k must be in the range 2 <= k <= 5')
    rtfr = 0;
    return;
end;

T = 1/fs;
nf = length(f);
N = length(x);
rtfr = zeros(nf, N);

Omega = f(2)-f(1); % frequency spacing

% Balance time-frequency resolution in number of samples, eq. 16
sigma_p = (sqrt(Omega)*factorial(k-1))/(sqrt(2*pi*T)*(k-1)^(k-1)*exp(-(k-1)));
coef = cell(10,1);

for j = 1:nf;
    omega_p = 2*pi*f(j);
    p = -sigma_p + i*omega_p; % eq. 6
 
    a = exp(p*T);
    
    % Coefficients from eq. 14, k = 2...5
    % Numerator coefficients in coef{i}, i = 2n+1, n = 1,2,3,4
    % Demoninator coefficients in coef{i}, i = 2n, n = 1,2,3,4,5

    coef{3} = sigma_p^2*T^2*[0 a];
    coef{4} = [1 -2*a a.^2];
    coef{5} = sigma_p^3*T^3*[0 1/2*a 1/2*a.^2];
    coef{6} = [1 -3*a 3*a.^2 -a.^3];
    coef{7} = sigma_p^4*T^4*[0 1/6*a 2/3*a.^2 1/6*a.^3];
    coef{8} = [1 -4*a 6*a.^2 -4*a.^3 a.^4];
    coef{9} = sigma_p^5*T^5*[0 1/24*a 11/24*a.^2 11/24*a.^3 1/24*a.^4];
    coef{10}= [1 -5*a 10*a.^2 -10*a.^3 5*a.^4 -a.^5];
       
    rtfr(nf-j+1,:) = filter(coef{2*k-1},coef{2*k},x);  % Eq. 16
end;

rtfr = abs(rtfr);