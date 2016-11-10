function [fh, ih] = implot(img, mag, map, varargin)
%IMPLOT  - scaled image plotting
%
%	 usage:  [fh, ih] = implot(img, mag, map, range)
%
% Use this function to plot an image scaled correctly
% to its source dimensions
%
% [rows x cols x 1] IMG arrays are plotted using 
% 	specified colorMAP (default = gray)
%
% [rows x cols x 3] images are plotted as truecolor
%
% optional MAG argument performs bicubic interpolation
% 	to specified magnification (default == 1)
%
% optional RANGE argument specifies imagesc intensity scaling
%
% optionally returns created figure, image handle

% mkt 07/97

%	parse args

if nargin < 1,
	eval('help implot');
	return;
end;
if nargin < 2 | isempty(mag), mag = 1; end;
if any(size(mag)~=[1 1] | mag<=0),
	error('magnification argument must be scalar value > 0');
end;
if nargin < 3, map = []; end;


%	convert source to double if necessary

if ~isa(img, 'double') & length(size(img))<3,
	img = double(img);
end;


%	make new figure, size it

[rows, cols, colors] = size(img);
magRows = floor(rows * mag);
magCols = floor(cols * mag);
pos = get(0, 'defaultfigureposition');
fh = figure('units', 'pixels', ...
		'resize', 'off', ...
		'paperpositionmode', 'auto', ...
		'position', [pos(1) pos(2)+pos(4)-magRows magCols magRows]);
if get(0, 'ScreenDepth') < 16, nGrays = 64; else, nGrays = 256; end;
if isempty(map), map = gray(nGrays); end;		% put here since gray generates a figure if none exists
set(fh, 'colormap', map);


%	magnify if necessary

if mag ~= 1,
	if mag < 1,						% build anti-aliasing filter
		h = fir1(10, magRows/rows)' * fir1(10, magCols/cols);
	end;
	if colors > 1,					% handle each plane of truecolor images
		urImg = img;
		img = repmat(uint8(0), [magRows magCols 3]);
		for rgb = 1 : 3,
			i = double(urImg(:,:,rgb));
			if mag < 1,				% filter if shrinking
				i = filter2(h, i);
			end;					% interpolate
			img(:,:,rgb) = uint8(round(resize(i,mag,'cubic')));
		end;
	else,							% single plane image
		if mag < 1,					% filter if shrinking
			img = filter2(h, img);
		end;
		img = resize(img, mag, 'cubic');
	end;
end;


%	plot

ih = imagesc(img, varargin{:});			% plot image; range passed if specified
set(gca, 'position', [0 0 1 1], 'xtick',[],'ytick',[]);
% zoom on

if nargout < 1,
	clear fh ih
end;
