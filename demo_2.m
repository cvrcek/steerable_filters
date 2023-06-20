%% load image
data = load('im.mat');
im = data.im;
phi = data.phi;
n = data.n;

    
figure(1);
imagesc(im);
%% use spatial filter (poisson)
N = 0:2:8;
opt = struct;
opt.type = 'pois';
opt.maskHalfSize = 25;
opt.scl = 1;
opt.doitbyFFT = false;

sqf = steerableMapGeneralized(im, N, opt);
r = reshape(mmsteermap('line', sqf, N, phi), n, n);
imagesc(r);
%% spatial, log gabor
N = 0:2:8;
opt = struct;
opt.type = 'loggabor';
opt.doitbyFFT = false;

opt.maskHalfSize = 25;
opt.sigma = 0.52877;
opt.f0 = 1/18;

sqf = steerableMapGeneralized(im, N, opt);
r = reshape(mmsteermap('line', sqf, N, phi), n, n);
imagesc(r);
%% fourier
N = 0:2:8;
opt = struct;
opt.type = 'loggabor';

opt.maskHalfSize = 25;
opt.sigma = 5.52877;
opt.f0 = 1/18;
opt.doitbyFFT = true;

sqf = steerableMapGeneralized(im, N, opt);
r = reshape(mmsteermap('line', sqf, N, phi), n, n);
imagesc(r);