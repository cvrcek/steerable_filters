function c = mmsteermap(type,sqf,N,angle)
 %
 % c = mmsteermap(type,sqf,N,angle)
 %
 %  type - feature type
 %   sqf - SQF response map from mmsteerablemap
 %  Nmax - maximum order (must be the same as for mmsteerablemap)
 % angle - orientation angle [degrees], in -90:90
 %         the positive angle is measured from the first axis (i, vertical,
 %         from top down) towards the second axis (j, horizontal from left
 %         to right)
 %
 %     c - correlation map steered to the angle
 %
 % See also: mmsteerablemap
 
 % (c) Radim Sara (sara@cmp.felk.cvut.cz) FEE CTU Prague, 03 Jun 13
 
 
 tL = mmtemplate(type,N);
 
 r = exp(1i*N*angle/180*pi); % rotation vector
 t = r.*tL; % rotated template
 %tnrm = t*t' - t(1)^2; % this should be the same as tL*tL'-tL(1)^2
 % TODO: make the mmtemplate already normalized
 
 
 %sqfr = reshape(shiftdim(sqf,2),size(sqf,3),[]);
 if iscell(sqf)
     c = zeros(size(sqf{1}));
     for i = 2:5         
         c = c + real(t(i))*real(sqf{i}) + imag(t(i))*imag(sqf{i});
     end
 else
    c = real(t)*real(sqf) + imag(t)*imag(sqf) - t(1)*sqf(1,:);
 end
 
%- x = repmat(shiftdim(t,-1),[size(sqf,1), size(sqf,2), 1]);
 %sqfnrm = sqf(:)'*sqf(:) - sum(sum(sqf(:,:,1).^2));
 %n = size(sqf,3);
 %sqfnrm = sum(conj(sqf).*sqf, 3) - sqf(:,:,1).^2; % this could be done at the beginning (37% of time)
 
 %covIMt = (sum((conj(sqf).*x+conj(x).*sqf)/2, 3) - t(1)*sqf(:,:,1));

%- covIMt = sum(real(sqf).*real(x) + imag(sqf).*imag(x),3) - t(1)*sqf(:,:,1);
 
 %c = covIMt./sqrt(sqfnrm)/sqrt(tnrm); % correlation
 
end

% x = rand(3,4) + 1i*rand(3,4);
% y = rand(3,4) + 1i*rand(3,4);
% 
% sum(sum((conj(x).*y + conj(y).*x)/2))
% (x(:)'*y(:) + y(:)'*x(:))/2



