function [h] = sqfGeneralized(N, m, opt)
% Convolution kernels for spherical quadrature filters
%	INPUT 
%  N    - SQF order list, eg 0:Nmax
%  m    - kernel half-size, the kernel size will be 2m+1 x 2m+1 1
%	OUTPUT 
% author: Vojtech Cvrcek, email: cvrcekv@gmail.com, date: 15-Oct-2019


 M = 3*(2*m+1);
 [u1, u2] = ndgrid(0:M-1, 0:M-1);
 mhalf = floor(M/2);
 u1(u1 > mhalf) = u1(u1 > mhalf) - M;
 u2(u2 > mhalf) = u2(u2 > mhalf) - M;
 r2 = u1.^2 + u2.^2;
 r2(1,1) = 1; % safe division
 
 
 switch opt.type
    case 'pois'
         s = opt.scl;
         [v1,v2] = ndgrid(-mhalf:mhalf,-mhalf:mhalf);
         pois = s/(2*pi)*(v1.^2 + v2.^2 + s^2).^(-3/2);
         pois = pois/sum(pois(:)); % unit integral
         ibf_fft = fft2(pois);
    case 'cauchy'
        ibf_fft = r2.^(opt.f0*opt.sigma/2).*exp(-opt.sigma*sqrt(r2));        
    case 'gaussderiv'
        ibf_fft = r2.^((opt.f0*opt.sigma)^2/2).*exp(-opt.sigma^2*r2);        
    case 'loggabor'
        ibf_fft = exp(-log(sqrt(r2)/opt.f0).^2/(2*log(opt.sigma/opt.f0)^2));
        %f = zeros(M);
        %f((M+1)/2, (M+1)/2) = 1;
        %ibf_fft = ibf_fft.*fft2(f); % centering ft        
    otherwise
        error('Filter type must be pois, cauchy, gaussderiv, or loggabor.')
 end

if strcmp(opt.type, 'cauchy') || strcmp(opt.type, 'gaussderiv') || strcmp(opt.type, 'loggabor')
    f = zeros(M);
    f((M+1)/2, (M+1)/2) = 1;
    ibf_fft = ibf_fft.*fft2(f); % spatial centering ft        
end
 
 
if 0
    f = zeros(M);
    f((M+1)/2, (M+1)/2) = 1;
    ibf_fft_alt = ibf_fft.*fft2(f);     
   % ibf_fft = fftshift(ibf_fft); % centering ft
    ibf = ifft2(ibf_fft);
    ibf_alt = ifft2(ibf_fft_alt);
    figure(1);
    imagesc(abs(ibf));
    figure(2);
    imagesc(abs(ibf_alt));
    
    keyboard;
end

 
 h = cell(1,length(N));
 for ni = 1:length(N) % Hilbert operator order
  n = N(ni);
  % Hilbert transform operator
  ght_fft = ( 1i*(u1 + 1i*u2)./sqrt(r2) ).^n;
  
  H = ght_fft.*ibf_fft; % SQF - riesz transformation .* isometric bandwidth filter 
  h{ni} = spar(ifft2(H),[3 3],'sum');
 end
  
%  A = NaN(length(h));
%  for k1 = 1:length(h)
%   for k2 = 1:length(h)
%    A(k1,k2) = sum(sum((conj(h{k1}).*h{k2} + conj(h{k2}).*h{k1})/2));
%   end
%  end
% 
%  figure
%  imager(A)
%  title('orhogonality matrix')
end
