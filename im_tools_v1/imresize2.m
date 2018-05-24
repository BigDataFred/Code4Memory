function [rout,g,b] = imresize(arg1,arg2,arg3,arg4,arg5,arg6)
%IMRESIZE Resize image.
%	B = IMRESIZE(A,M,'method') returns an image matrix that is 
%	M times larger (or smaller) than the image A.  The image B
%	is computed by interpolating using the method in the string
%	'method'.  Possible methods are 'nearest', 'bilinear', 
%	'bicubic', or 'spline'. B = IMRESIZE(A,M) uses 'nearest' when
%	for indexed images and 'bilinear' for intensity images.
%
%	B = IMRESIZE(A,[MROWS NCOLS],'method') returns a matrix of 
%	size MROWS-by-NCOLS.
%
%	[R1,G1,B1] = IMRESIZE(R,G,B,M,'method') or 
%	[R1,G1,B1] = IMRESIZE(R,G,B,[MROWS NCOLS],'method') resizes
%	the RGB image in the matrices R,G,B.  'bilinear' is the
%	default interpolation method.
%
%	When the image size is being reduced, IMRESIZE lowpass filters
%	the image before interpolating to avoid aliasing. By default,
%	this filter is designed using FIR1, but can be specified using 
%	IMRESIZE(...,'method',H). The default filter is 11-by-11.
%	IMRESIZE(...,'method',N) uses an N-by-N filter.
%	IMRESIZE(...,'method',0) turns off the filtering.
%	Unless a filter H is specified, IMRESIZE will not filter
%	when 'nearest' is used.
%	
%	See also IMZOOM, FIR1, INTERP2.

%	Clay M. Thompson 7-7-92
%	Copyright (c) 1992 by The MathWorks, Inc.
%	$Revision: 1.11 $  $Date: 1993/09/01 22:25:27 $

error(nargchk(2,6,nargin));

% Trap imresize(r,b,g,...) calls.
if nargin==4,
  if ~isstr(arg3), % imresize(r,g,b,m)
    r = imresize(arg1,arg4,'bil');
    g = imresize(arg2,arg4,'bil');
    b = imresize(arg3,arg4,'bil');
    if nargout==0, imshow(r,g,b), else, rout = r; end
    return
  end
elseif nargin==5, % imresize(r,g,b,m,'method')
  r = imresize(arg1,arg4,arg5);
  g = imresize(arg2,arg4,arg5);
  b = imresize(arg3,arg4,arg5);
  if nargout==0, imshow(r,g,b), else, rout = r; end
  return
elseif nargin==6, % imresize(r,g,b,m,'method',h)
  r = imresize(arg1,arg4,arg5,arg6);
  g = imresize(arg2,arg4,arg5,arg6);
  b = imresize(arg3,arg4,arg5,arg6);
  if nargout==0, imshow(r,g,b), else, rout = r; end
  return
end

% Determine default interpolation method
if nargin<3, 
  if isgray(arg1), case_method = 'bil'; else case_method = 'nea'; end,
else
  method = [lower(arg3),'   ']; % Protect against short method
  case_method = method(1:3);
end

if prod(size(arg2))==1,
  bsize = floor(arg2*size(arg1));
else
  bsize = arg2;
end

if any(size(bsize)~=[1 2]),
  error('M must be either a scalar multiplier or a 1-by-2 size vector.');
end

if nargin<4,
  nn = 11; h = []; % Default filter size
else
  if length(arg4)==1, nn = h; h = []; else nn = 0; h = arg4; end
end

[m,n] = size(arg1);

if nn>0 & case_method(1)=='b',  % Design anti-aliasing filter if necessary
  if bsize(1)<m, h1 = fir1(nn-1,bsize(1)/m); else h1 = 1; end
  if bsize(2)<n, h2 = fir1(nn-1,bsize(2)/n); else h2 = 1; end
  if length(h1)>1 | length(h2)>1, h = h1'*h2; else h = []; end
end

if ~isempty(h), % Anti-alias filter A before interpolation
  a = filter2(h,arg1);
else
  a = arg1;
end

if case_method(1)=='n', % Nearest neighbor interpolation
  dx = n/bsize(2); dy = m/bsize(1); 
  uu = (dx/2+.5):dx:n+.5; vv = (dy/2+.5):dy:m+.5;
elseif   all(case_method == 'bil') ...
       | all(case_method == 'bic') ...
       | all(case_method == 'spl'),
  uu = 1:(n-1)/(bsize(2)-1):n; vv = 1:(m-1)/(bsize(1)-1):m;
else
  error(['Unknown interpolation method: ',method]);
end

%
% Interpolate in blocks
%
nu = length(uu); nv = length(vv);
blk = bestblk([nv nu]);
nblks = floor([nv nu]./blk); nrem = [nv nu] - nblks.*blk;
mblocks = nblks(1); nblocks = nblks(2);
mb = blk(1); nb = blk(2);

rows = 1:blk(1); b = zeros(nv,nu);
for i=0:mblocks,
  if i==mblocks, rows = (1:nrem(1)); end
  for j=0:nblocks,
    if j==0, cols = 1:blk(2); elseif j==nblocks, cols=(1:nrem(2)); end
    if ~isempty(rows) & ~isempty(cols)
      [u,v] = meshgrid(uu(j*nb+cols),vv(i*mb+rows));
      % Interpolate points
      if case_method(1) == 'n', % Nearest neighbor interpolation
        b(i*mb+rows,j*nb+cols) = interp2(a,u,v,'nearest');
      elseif all(case_method == 'bil'), % Bilinear interpolation
        b(i*mb+rows,j*nb+cols) = interp2(a,u,v,'linear');
      elseif all(case_method == 'bic'), % Bicubic interpolation
        b(i*mb+rows,j*nb+cols) = interp2(a,u,v,'bicubic');
      elseif all(case_method == 'spl'), % spline interpolation
        b(i*mb+rows,j*nb+cols) = interp2(a,u,v,'spline');
      end
    end
  end
end

if nargout==0,
  if isgray(b), imshow(b,size(colormap,1)), else imshow(b), end
  return
end

if isgray(arg1), rout = max(0,min(b,1)); else rout = b; end

