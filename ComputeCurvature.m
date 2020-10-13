function [k,nInfl,mci] = ComputeCurvature(xy, varargin)
%COMPUTECURVATURE  - compute signed curvature of UTI contour
%
%	usage:  [k, NINFL, MCI] = ComputeCurvature(xy, ...)
%
% given contour XY [nPts x X,Y] this procedure computes its signed curvature K
% using a gradient method: (dx .* ddy - dy .* ddx) ./ (dx.^2 + dy.^2).^1.5
%
% specify nonzero CIRCLE parameter to compute curvature instead as the reciprocal 
% of circle radii fit through successive point triplets
%
% truncates (and generates a warning for) any contours that "curl over"; i.e., for which x
% is non-monotonic and y is higher on the open end
%
% plots data and associated curvature if no output arguments requested
%
% optionally returns the Preston et al. (2019) NINFL measure (number of inflection points) from 
% trimmed signed curvature (values whose associated radius is less than TRIM * path integral); 
% an inflection point is counted if the trimmed curvature changes sign, and there is always at 
% least one inflection point for non-collinear points
%
% optionally returns the Dawson et al. (2016) modified curvature index (MCI), the integral of
% the filtered (5-tap Butterworth) unsigned curvature
%
% optional supported 'NAME',VALUE parameter pairs (defaults in {}):
%   CIRCLE - if nonzero compute curvature as the reciprocal of circle radi {0}
%   NPTS   - XY resampled to NPTS equally spaced points along its arc (0 disables; {100})
%   FCLP   - MCI low-pass filter cutoff (0 disables; {.25})
%   TRIM   - NINFL curvature trim factor (0 disables; {.3})
%
% examples
%   compute curvature without resampling:
% k = ComputeCurvature(xy, 'NPTS',0);
%
%   modify TRIM factor in computing NINFL
% [~,NINFL] = ComputeCurvature(xy, 'TRIM',.25)
%
%   compute MCI using circle method
% [~,~,MCI] = ComputeCurvature(xy, 'CIRCLE',1);
%
% Dawson K, Tiede M & Whalen D. (2016). Methods for quantifying tongue shape and complexity using 
% ultrasound imaging. Clinical Linguistics & Phonetics, 30(3-5), 328-344.
%
% Preston J, McCabe P, Tiede M & Whalen D. (2019). Tongue shapes for rhotics in school-age children
% with and without residual speech errors. Clinical Linguistics & Phonetics, 33(4), 334-348.

% specify nonzero NOISY for diagnostics

% mkt 06/15
% mkt 07/17 add mci
% mkt 05/18 add circle k, reorganize
% mkt 10/20 rationalize output plotting

% parse args
if nargin < 1, eval('help ComputeCurvature'); return; end;
[m,n] = size(xy);
if min([m,n]) > 2, error('expecting [nPts x X,Y] array for XY'); end;
if n > m, xy = xy'; end;

circle = 0;
FcLP = .25;
nPts = 100;
trim = .3;
noisy = 0;
for ai = 2 : 2 : length(varargin),
	switch upper(varargin{ai-1}),
		case 'CIRCLE', circle = varargin{ai};
		case 'FCLP', FcLP = varargin{ai};
		case 'NPTS', nPts = varargin{ai};
		case 'TRIM', trim = varargin{ai};
		case 'NOISY', noisy = varargin{ai};
		otherwise, error('unrecognized parameter (%s)', varargin{ai-1});
	end;
end;

% resample to nPts equally spaced points along arc unless NPTS==0
dist = [0 ; cumsum(sqrt(sum(diff(xy).^2,2)))];
if nPts > 0, xy = interp1(dist,xy,linspace(0,dist(end),nPts)','spline'); end;

if noisy, xy0 = xy; end;

% test for "curl over" non-monotonicity
q = find(sign(diff(xy(:,1))) < 1);
if length(q)+1 == size(xy,1),		% all negative:  flip curve to be left -> right
	xy = flipud(xy);
	q = find(sign(diff(xy(:,1))) < 1);
	if noisy, fprintf('all negative (flipped)\n'); end;
end;
if ~isempty(q),
	if length(q) > size(xy,1)/2, 	% more than half negative:  flip curve to be left -> right
		xy = flipud(xy); 
		q = find(sign(diff(xy(:,1))) < 1); 
		if noisy, fprintf('more than half negative (flipped)\n'); end;
	elseif noisy,
		fprintf('%d non-monotonic points detected\n', length(q));
	end;
	if q(1) == 1,
		n = find(diff(q)>1);
		if isempty(n), n = length(q); end;
		if xy(1,2) < xy(q(end),2), 
			xy(q(1:n),:) = [];
			fprintf('ComputeCurvature:  %d leading non-monotonic overcurl points deleted\n', n); 
			q = find(sign(diff(xy(:,1))) < 1);
		else,
			q(1:n) = [];
			if noisy, fprintf('%d leading non-monotonic points retained\n', n); end;
		end;
	end;
	if ~isempty(q),
		if xy(end,2) < xy(q(1),2) && xy(end,1) < xy(max(q),1),
			xy(q(1):end,:) = []; 
			fprintf('ComputeCurvature:  %d trailing non-monotonic overcurl points deleted\n', length(q));
		elseif noisy,
			fprintf('%d trailing non-monotonic points retained\n', length(q));
		end;
	end;
end;

% find signed curvature using circle radius method
if circle,
	mag = @(v) (sqrt(sum(v.^2,2)));
	V1 = [xy(1:end-2,:)-xy(3:end,:),zeros(size(xy,1)-2,1)];
	V2 = [xy(2:end-1,:)-xy(3:end,:),zeros(size(xy,1)-2,1)];
	V12 = cross(V1,V2);
	k = 2 * V12(:,3) ./ (mag(V1) .* mag(V2) .* mag(V1-V2));	

% find signed curvature using central differencing
else,
	dx = gradient(xy(:,1)); dy = gradient(xy(:,2));
	ddx = gradient(dx); ddy = gradient(dy);
	k = (dx .* ddy - dy .* ddx) ./ (dx.^2 + dy.^2).^1.5;
end;

% trim curvature to values whose associated radius is less than TRIM * path integral from first to last point
fk = k; 
if trim > 0,
	q = sum(sqrt(sum(diff(xy).^2,2))) * trim;
	fk(abs(1./k) > q) = 0;
end;

% count inflections (nonzero sign changes) plus one
sfk = sign(fk);
xfk = sfk(sfk ~= 0);
if isempty(xfk),
	c = corrcoef(xy(:,1),xy(:,2));
	c = 1 - c(2,1);
	if c < .01,				% correlation > .99
		nInfl = 0;			% collinear points
	else,
		nInfl = 1;			% curvature below threshold
	end;
else,
	nInfl = sum(diff(xfk)~=0) + 1;
end;

% compute MCI (average of Simpson's Rule applied to intervals 1:N-1 plus trapezoid for final interval
% and 2:N plus trapezoid for first interval) on filtered curvature (kk)
n = length(k);
rk = flipud(k);
kk = [rk;k;rk];		% pad to avoid edge effects
if FcLP > 0,
	[b,a] = butter(5,FcLP);
	kk = filtfilt(b,a,kk);
end; 
kk = kk(n+1:n*2);	% filtered curvature
x = [0 ; cumsum(sqrt(sum(diff(xy).^2,2)))];
y = abs(kk);
if mod(n,2),		% even number of intervals
	mci = sr(x,y);
else,				% odd number of intervals
	z1 = sr(x(1:n-1),y(1:n-1)) + .5*diff(x(n-1:n))*sum(y(n-1:n));
	z2 = sr(x(2:n),y(2:n)) + .5*diff(x(1:2))*sum(y(1:2));
	mci = mean([z1,z2]);
end;

if nargout > 0, return; end;

% map inflection points to curvature for plotting
if nInfl < 2,
	N = [];
else,
	z = FindExtents(find(sfk == 0)); 
	if sfk(z(1,1)) == 0,
		sfk(z(1,1):z(1,2)) = sfk(z(1,2)+1);
		z = FindExtents(find(sfk == 0)); 
	end;
	for zi = 1 : size(z,1), 
		sfk(z(zi,1):z(zi,2)) = sfk(z(zi,1)-1); 
	end;
	N = find(diff(sfk))';
end;

% plot
figure; 
subplot(211);
plot(xy(:,1),xy(:,2),'b-');
if ~isempty(N), 
	hold on; 
	plot(xy(N,1),xy(N,2),'ro'); hold on; plot(xy(N,1),xy(N,2),'g*');
end;
axis equal; 
yl = get(gca,'ylim');
r = .05*diff(yl);
set(gca,'ydir','reverse','ylim',r*[-1 1]+yl);
title(inputname(1),'interpreter','none')
subplot(212);
h = plot([k,fk,abs(kk)]);
h(1).Color = 'b'; h(2).Color = 'r'; h(3).Color = [0 .7 0];
set(gca,'xlim',[1 length(k)]);
if ~isempty(N), 
	line([N;N], get(gca,'ylim'), 'color','g', 'linewidth',2); 
end;
line([1 length(k)],[0 0],'color',[.7 .7 .7],'linestyle',':');
title(sprintf('# inflections:  %d       MCI = %.2f', nInfl, mci));
legend(h,'Curvature (K)','Trimmed K','Filtered abs(K)');
clear k

%===================================================================================================
% FINDEXTENTS  - find continguous extents of an indexed signal

function idx = FindExtents(v)

idx = find(diff([-1;v]) > 1);
len = diff([idx;length(v)+1]);
idx = [v(idx) , v(idx)+len-1];			% [nExtents x head,tail]


%===================================================================================================
% SR  - integrate using quadratic interpolation (Simpson's Rule)

function z = sr(x,y)

n = length(y);
dx = diff(x);
dx1 = dx(1:end-1);
dx2 = dx(2:end);

alpha = (dx1+dx2)./dx1/6;
a0 = alpha.*(2*dx1-dx2);
a1 = alpha.*(dx1+dx2).^2./dx2;
a2 = alpha.*dx1./dx2.*(2*dx2-dx1);

z = sum(a0(1:2:end).*y(1:2:n-2) + a1(1:2:end).*y(2:2:n-1) + a2(1:2:end).*y(3:2:n),1);
