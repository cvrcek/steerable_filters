%       Image processing - neighborhood operation  with reduction
%                                
%       Author:         Vladimir Smutny 
%
%       Language:       Matlab 4.2, (c) MathWorks
%                                                                
%       15.10.94  Computer Vision Laboratory CTU Prague
%
%       function OUT=spar(IMAGE,[ROWMASKSIZE,COLMASKSIZE],FUN)
%             to the input image is applied the operation on the neighborhood
%             size of the output image is reduced by the factor ROWMASKSIZE and
%             COLMASKSIZE (i.e. mask is not overlapped when applied)
%       IMAGE   = input image
%       ROWMASKSIZE = size of the mask, optional, implicit value is 2 
%       COLMASKSIZE = size of the mask, optional, implicit value is ROWMASKSIZE 
%       FUN     = method from the list: 
%                    'min','max'     ... minimum, maximum
%                    'first'         ... first value (upper left pixel)
%		     'central'       ... central value (central (odd size of mask) 
%					 or slightly upper left (even size of mask) pixel)
%                    'aver', 'gaver' ... average, geometric average
%                    'sum', 'prod'   ... sum, product
%                    'median'        ... median 
%                    'mednan'        ... median ignoring NaNs (slow)
%                    'disp'          ... dispersion (sum(y.^2)/n-(sum(y)/n)^2)
%                 optional parameter, implicit is 'first',
%       OUT     = output image
%       
%       see also INTERP4(IM,STEP), MASKFILT

% updated for rgb images by R. Sara Jan 19, 2000
% median ignoring NaNs by R. Sara, Sep 5, 2001

function out=spar(InputImage,MaskSize,fun)


[m,n]=size(InputImage);



if (nargin == 1) 
        RowMaskSize=2;
        ColMaskSize=2;
        fun = 'first';
elseif  (nargin == 2) 
  if isstr(MaskSize)
        fun = MaskSize;
        RowMaskSize=2;
        ColMaskSize=2;
  else        
        RowMaskSize=MaskSize(1);
        if length(MaskSize)==1
          ColMaskSize=MaskSize(1);
        else
          ColMaskSize=MaskSize(2);
        end
        fun = 'first';
  end
elseif  (nargin == 3)
  RowMaskSize=MaskSize(1);
  if length(MaskSize)==1
    ColMaskSize=MaskSize(1);
  else
    ColMaskSize=MaskSize(2);
  end
end


switch ndims(InputImage)
 case {1,2}
  % just continue
 case 3
  out(:,:,1) = spar(InputImage(:,:,1),[RowMaskSize,ColMaskSize],fun);
  out(:,:,2) = spar(InputImage(:,:,2),[RowMaskSize,ColMaskSize],fun);
  out(:,:,3) = spar(InputImage(:,:,3),[RowMaskSize,ColMaskSize],fun);
  return
 otherwise
  error 'Matrices of dimension larger than 3 are not supported'
end

% followInputImageg  cuts the size of the process image only to valid size of the result
m = RowMaskSize * floor(m/RowMaskSize);
n = ColMaskSize * floor(n/ColMaskSize);

InputImage = InputImage(1:m,1:n);

if fun(1) == 'c'
  out=InputImage(ceil(RowMaskSize/2):RowMaskSize:m,ceil(ColMaskSize/2):ColMaskSize:n);
  return
end

if fun(1) == 'f'
  out=InputImage(1:RowMaskSize:m,1:ColMaskSize:n);
  return
end

out=zeros(m/RowMaskSize,n/ColMaskSize);

if fun(1) =='a'
  for i=1:RowMaskSize
    for j=1:ColMaskSize
      out=out+InputImage(i:RowMaskSize:m,j:ColMaskSize:n);
    end
  end
  out=out/(RowMaskSize*ColMaskSize);
  return
end  

if fun(1) =='d'
  sum =zeros(m/RowMaskSize,n/ColMaskSize);
  qsum=zeros(m/RowMaskSize,n/ColMaskSize);
  for i=1:RowMaskSize
    for j=1:ColMaskSize
      sum = sum+InputImage(i:RowMaskSize:m,j:ColMaskSize:n);
      qsum=qsum+InputImage(i:RowMaskSize:m,j:ColMaskSize:n).^2;
    end
  end
  out=qsum/(RowMaskSize*ColMaskSize)-(sum/(RowMaskSize*ColMaskSize)).^2;
  return
end  

if fun(1) =='s'
  for i=1:RowMaskSize
    for j=1:ColMaskSize
      out=out+InputImage(i:RowMaskSize:m,j:ColMaskSize:n);
    end
  end
  return
end  

if fun(1) =='g'
  out=ones(m/RowMaskSize,n/ColMaskSize);
  for i=1:RowMaskSize
    for j=1:ColMaskSize
      out=out.*InputImage(i:RowMaskSize:m,j:ColMaskSize:n);
    end
  end
  out=out^(1/(RowMaskSize*ColMaskSize));
  return
end  

if fun(1) =='p'
  out=ones(m/RowMaskSize,n/ColMaskSize);
  for i=1:RowMaskSize
    for j=1:ColMaskSize
      out=out.*InputImage(i:RowMaskSize:m,j:ColMaskSize:n);
    end
  end
  return
end  

transpose = 0;
if m < n
  transpose = 1;
  temp = n;
  n    = m;
  m    = temp;
  InputImage   = InputImage';
  out  = out';
  temp = RowMaskSize;
  RowMaskSize  = ColMaskSize;
  ColMaskSize    = temp;
end  
 
if strcmpi(fun, 'mednan')
 InputImage = nan2x(InputImage,Inf);
end

if fun(1:3) == 'med'

  ind = [1:RowMaskSize]' * ones(1,ColMaskSize) + ones(RowMaskSize,1)* [0:m:(ColMaskSize-1)*m];

  ind = ind(:)*ones(1,n/ColMaskSize) + ones(RowMaskSize*ColMaskSize,1)*[0:m*ColMaskSize:m*(n-ColMaskSize)];
  ind = ind(:);

  for i=1:m/RowMaskSize
    c = reshape(InputImage((i-1)*RowMaskSize+ind),RowMaskSize* ...
		ColMaskSize,n/ColMaskSize);
    out(i,1:n/ColMaskSize) = median(c);
  end
  
  if transpose out = out';end

  if lower(fun) == 'mednan'
   out = x2nan(out,Inf);
  end
  return
end


% MIN and MAX case  
  
  
 
  
ind = [1:ColMaskSize]' * ones(1,n/ColMaskSize) + ones(ColMaskSize,1) * [0:ColMaskSize:n-ColMaskSize];

for i=1:m/RowMaskSize
  eval(['c = ', fun, '(InputImage((i-1)*RowMaskSize+1:i*RowMaskSize,:));']);
  eval(['out(i,1:n/ColMaskSize) = ', fun,'(reshape(c(ind),ColMaskSize,n/ColMaskSize));']);
end

if transpose out = out';end

