function [sqf] = sqfFFT(im, N, opt)
% generate spherical quadrature filters
%	INPUT 
%  N    - SQF order list, eg 0:Nmax
%  m    - kernel half-size, the kernel size will be 2m+1 x 2m+1 1
%
%   OUTPUT
%  sqf  - spherical quadrature filters
% author: Vojtech Cvrcek, email: cvrcekv@gmail.com, date: 17-Oct-2019

[m, n] = size(im);

[u1, u2] = ndgrid(0:m-1, 0:n-1);
mhalf = floor(m/2);
nhalf = floor(n/2);
u1(u1 > mhalf) = u1(u1 > mhalf) - m;
u2(u2 > nhalf) = u2(u2 > nhalf) - n;
r2 = u1.^2 + u2.^2;
r2(1,1) = 1; % safe division


switch opt.type
    case 'pois'
         s = opt.scl;
         [v1,v2] = ndgrid(-mhalf:mhalf,-mhalf:mhalf);
         pois = s/(2*pi)*(v1.^2 + v2.^2 + s^2).^(-3/2);
         pois = pois/sum(pois(:)); % unit integral
         ibf_fft = fft2(ifftshift(pois));
    case 'cauchy'
        ibf_fft = r2.^(opt.f0*opt.sigma/2).*exp(-opt.sigma*sqrt(r2));        
    case 'gaussderiv'
        ibf_fft = r2.^((opt.f0*opt.sigma)^2/2).*exp(-opt.sigma^2*r2);
    case 'loggabor'
        ibf_fft = exp(-log(sqrt(r2)/opt.f0).^2/(2*log(opt.sigma/opt.f0)^2));
    otherwise
        error('Filter type must be cauchy, gaussderiv, or loggabor.')
end

im_fft = fft2(im);
sqf = NaN(m, n, length(N));
for ni = 1:length(N) % Hilbert operator order
  % Hilbert transform operator
  ght_fft = ( 1i*(u1 + 1i*u2)./sqrt(r2) ).^N(ni);
  
  H = ght_fft.*ibf_fft; % SQF
  %h = ifft2(H);
  %hnorm = sqrt(h(:)'*h(:)); % the normalization should not be needed, all norms are the same
  sqf(:,:,ni) = ifft2(im_fft.*H);
  %sqf(:,:,ni) = ifftshift(ifft2(im_fft.*H));
end