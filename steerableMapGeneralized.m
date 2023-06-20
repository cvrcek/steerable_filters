function sqf = steerableMapGeneralized(ims, N, opt)
 %
 % sqf = mmsteerablemap(ims,N,scl,opt)
 %
 % Steerable edge correlation map
 %
 %    N - filter orders, eg 0:2:Nmax
 %
 %  sqf - Steerable features, response of image to spherical quadrature
 %        filter basis set.
 %        Reshaped for more efficient steered correlation computation
 %        when we reshape it back, we get:
 %        SQF = reshape(sqf,size(ims,1),[]);
 %        SQF(i,j,:) is the set of responses to all filters at pixel (i,j)
 %
 % See also: mmsteermap
% author: Vojtech Cvrcek, email: cvrcekv@gmail.com, date: 15-Oct-2019
 lopt.doitbyFFT = true;
 lopt.maskHalfSize = 20;
 lopt.normalizeSQF = true;
 
 opt = optionmerge(lopt,'opt');
  
 if opt.doitbyFFT
  sqf = sqfFFT(ims,N, opt);
 else
  sqf = steerablemapConv(ims, N, opt);  % by convolution
 end
 
 if opt.normalizeSQF
  sqfnrm = sqrt(sum(conj(sqf).*sqf, 3) - sqf(:,:,1).^2);
  sqf = sqf./repmat(sqfnrm,[1,1,size(sqf,3)]);
 end
 
 % reshape for faster computations
 sqf = reshape(shiftdim(sqf,2),size(sqf,3),[]);
end

function sqf = steerablemapConv(ims, N, opt)

 h = sqfGeneralized(N, opt.maskHalfSize, opt);

 img = repeat(ims,size(h{1}));
 sqf = NaN(size(ims,1),size(ims,2),length(N));
 for ni=1:length(N)
  sqf(:,:,ni) = cutoff(conv2(img,h{ni},'same'), size(h{1}));
 end
end