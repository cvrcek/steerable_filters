%% load image
n = 501;
if 0
    rs0 = [10, 10, 10];
    rs = {rs0, rs0, rs0};
    lens0 = [100, 50, 25];
    lens = {lens0, lens0, lens0};
    [ims, streaks, streaks_mid_gt, streaks_idx_gt, phis_gt, star_locs] = ...
                                         artificial_image(...
                                            'num_of_streaks', 3,...
                                            'num_of_stars', 0,...
                                            'streak_amplitudes', [1.5, 2.3, 2],...
                                            'stars_amplitudes', [],...
                                            'rs', rs,...
                                            'lens', lens,...
                                            'mn', [n, n]...
                                            );
    im = ims{2};
    save('im.mat', 'im');
elseif 1
    [ims, streaks, streaks_mid_gt, streaks_idx_gt, phis_gt, star_locs] = ...
                                     artificial_image(...
                                     'mn', [n, n],...
                                     'streak_amplitudes', [1.5]);
    im = ims{2};
    phi = phis_gt{1}{1} + 90;
    save('im.mat', 'im', 'phi', 'n');
else
    data = load('im.mat');
    im = data.im;
    phi = data.phi;
    n = data.n;
end
    
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