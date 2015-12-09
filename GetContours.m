function varargout = GetContours(varargin)
%GETCONTOURS  - extract contours from US movie frames
%
%	usage:  GetContours(fName, ...)     % to initialize
%
% displays current FRAME from movie FNAME (default ".avi" ext)
% click in order along contour to place anchor points
% click and drag on a point to reposition it
% control (right) click on a point deletes it
% shift (mid) click on a point reports anchor point info
% double-click clears all points
% undo reverts to previous anchor positions
%
% use menu entries to adjust anchor points and image display
%
% supported optional 'NAME',VALUE parameter pairs:
%   CONFIG   - a struct specifying display properties to modify; valid fields are
%              DOTSIZE, DOTCOLOR, LINEWIDTH, LINECOLOR
%   IMGMOD   - procedure applied to every image before loading (e.g. {@Rescale,1.5})
%   KEYFRAMES- key frame list (default none)
%   NPOINTS  - number of output contour points (default 100)
%   SAVEIMG  - save frame images with matched contours ('T' or 'F'; default 'F')
%   TEXTGRID - Praat TextGrid and tier to parse for key frames (see examples)
%               note that only labeled intervals are used but all point labels used
%   TRACKER  - automatic contour-fitting procedure (e.g. {@Snake})
%   VNAME    - output variable name (defaults to FNAME)
%
% VNAME is an array-of-structs for each frame with labeled contours, with fields
%   FRAME    - frame number
%   XY       - ULC-based image coordinates describing contour [nPoints x X,Y] 
%   ANCHORS  - associated contour spline anchor points [nAnchors x X,Y] 
%   NOTE     - optional annotation
%
% if VNAME exists as a variable in the workspace its values are loaded and it is 
% renamed to VNAME_old (VNAME is updated with values from the new session)
%
% the current position of the XY contour and ANCHORS for that FRAME are updated
% to VNAME when the current frame is changed or when closing the window
%
% when the window is closed, VNAME.mat is created within the current working directory
% containing variable VNAME with fields as above with additional
%   TIME     - frame offset in secs from start of movie
%   IMAGE    - grayscale image at current frame (if SAVEIMG was specified as 'T')
%
% contours, associated frame numbers and images may be exported to the workspace 
% before closing the GetContours window using
%   [contours,frames,images] = GetContours('EMIT');   % [nRows x nCols x nFrames] images
%   contour = GetContours('EMIT',FRAME);   % emit only specified frame(s)
%
% to extract contours from VNAME use 
%   contours = reshape(cell2mat({VNAME.XY}),[NPOINTS 2 length(VNAME)]);
%
% Example:  get contours and associated frames from VNAME 'foo' with NPOINTS == 50
%   GetContours('movie.avi', 'NPOINTS',50, 'VNAME','foo', 'SAVEIMG','T');
%   [contours,frames] = GetContours('EMIT');  % [nPoints x X,Y x nFrames]
% or equivalently
%   contours = reshape(cell2mat({foo.XY}),[50 2 length(foo)]);
%   frames = cell2mat({foo.FRAMES})
%
% Example:  specify key frames from point or interval tier "frame" in "foo.TextGrid"
%   GetContours('movie.avi', 'TEXTGRID',{'foo','frame'})
% if no tier name specified loads from the first tier found
%   GetContours('movie.avi', 'TEXTGRID','foo')
% specify key frames directly
%   GetContours('movie.avi', 'KEYFRAMES',[23 47 123 247 319])
% 
% Example:  scale image by 75%, use Snake active contour model for tracking
%   GetContours('movie.avi', 'IMGMOD',{@Rescale,.75}, 'TRACKER',@Snake);

% Copyright (C) 2015 mark tiede <tiede@haskins.yale.edu>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3 or any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% mark tiede fecit
% mkt 12/13 v0.4
% mkt 07/14 v0.5 handle corrupt movie frames, display interval labels
% mkt 08/14 v0.6 support annotation, fix keyframe issues
% mkt 10/15 v0.7 support for r2014b+
% mkt 12/15 v0.8

GCver = 'v0.8';	% current version

if nargin < 1,
	eval('help GetContours');
	return;
end;

% strip src,evt on callbacks
if ~ischar(varargin{1}), varargin(1:2) = []; end;

% branch by action
switch upper(varargin{1}),

%-----------------------------------------------------------------------------
% ABOUT:  display version and brief help

	case 'ABOUT',
		state = get(gcbf,'userData');
		blurb = sprintf('%s\n\n%s\n%s\n%s\n%s\n%s\n\n%s', ...
			'Extract contours from US movie frames:', ...
			'click in order along contour to place anchor points', ...
			'click and drag on a point to reposition it', ...
			'control-click on a point deletes it', ...
			'shift-click on a point reports its position', ...
			'double-click clears all points that frame', ...
			sprintf('Output variable name: %s', state.VNAME));
			
		width = 400;
		height = 300;
		fPos = get(gcf, 'position');
		pos = [fPos(1)+(fPos(3)-width)/2 , fPos(2)+(fPos(4)-height)/2 , width , height];

		cfg = dialog('name', 'About GetContours', ...
			'menubar', 'none', ...
			'position', pos);

		uicontrol(cfg, ...
			'Position', [width/2-150,height-50,300,30], ...
			'Style', 'text', ...
			'fontName', 'Arial', ...
			'fontSize', 24, ...
			'String', sprintf('GetContours  %s',GCver));
			
		uicontrol(cfg, ...
			'Position', [width/2-180,height-90,360,25], ...
			'Style', 'edit', ...
			'fontName', 'Arial', ...
			'fontSize', 12, ...
			'String', 'https://github.com/mktiede/GetContours');

		uicontrol(cfg, ...
			'Style', 'frame', ...
			'Position', [20 20 width-39 170]);
		if ismac, fs = 14; else, fs = 11; end;
		uicontrol(cfg, ...
			'Style', 'text', ...
			'HorizontalAlignment', 'left', ...
			'String', blurb, ...
			'fontName', 'Arial', ...
			'fontSize', fs, ...
			'Position', [25 25 width-46 160]);


%-----------------------------------------------------------------------------
% ANNOTATE:  add label to current frame

	case 'ANNOTATE',
		state = get(gcbf,'userData');
		v = evalin('base',state.VNAME);
		frames = cell2mat({v.FRAME});
		k = (state.CURFRAME == frames);
		ts = GetText(v(k).NOTE, frames(k));
		set(state.LH,'string',ts);
		v(k).NOTE = ts;
		assignin('base',state.VNAME,v);

		
%-----------------------------------------------------------------------------
% AVERAGING:  specify frame averaging

	case 'AVERAGING',
		state = get(gcbf,'userData');
		[enabled,win] = GetAvg(state);
		if ~isempty(enabled),
			state.USEAVG = enabled;
			state.AVG = win;
			set(state.IH,'cdata',GetImage(state));
			set(gcbf,'userData',state);
		end;
		

%-----------------------------------------------------------------------------
% CLOSE:  shutdown handler

	case 'CLOSE',
		state = get(gcbf,'userData');
		if ishandle(state.RH), delete(state.RH); end;		% close contrast adjustment (if open)
		delete(gcbf);										% close window

		v = evalin('base',state.VNAME);						% update output variable state
		frames = cell2mat({v.FRAME});
		k = (state.CURFRAME == frames);
		v(k).XY = state.XY;
		v(k).ANCHORS = state.ANCHORS;
		[n,k] = sort(frames);
		v = v(k);											% impose sequential order
		v(cellfun(@isempty,{v.XY})) = [];					% delete empty frames

% add additional annotation fields		
		sr = state.MH.FrameRate;
		for vi = 1 : length(v),
			v(vi).TIME = (v(vi).FRAME-1)/sr;							% frame offset in secs from start of movie
			if state.SAVEIMG,
				try,
					img = uint8(mean(GetMovieFrame(state.MH,v(vi).FRAME,state.FRATE),3));	% grayscale image at frame
					if ~isempty(state.IMGMODP), img = state.IMGMODP(img,state.IMGMODA); end;
					v(vi).IMAGE = img;							% grayscale image at frame
				catch,
					v(vi).IMAGE = zeros(get(state.MH,'Height'),get(state.MH,'Width'),'uint8');
				end
			else,
				v(vi).IMAGE = [];
			end;
		end;
		
		assignin('base',state.VNAME,v);						% update in workspace
		SaveVar(v, state.VNAME);							% write MAT file		


%-----------------------------------------------------------------------------
% DOWN:  mouseDown handler
%
% click in image creates new anchor point
% click on an anchor repositions it
% control-click on a point deletes it
% double-click clears all points

	case 'DOWN',
		state = get(gcbf,'userData');
		state.PREVAP = state.ANCHORS;
		gotPoint = (length(varargin) > 1);	% nonzero for click on existing point
		mod = get(gcbf,'selectionType');

		switch mod,		

			case 'normal',			% unmodified

% move existing point
				if gotPoint,		
			
					set(gcbf, 'windowButtonMotionFcn',{@GetContours,'MOVE',gcbo}, ...
								'windowButtonUpFcn',{@GetContours,'UP',gcbo}, ...
								'pointer','crosshair', ...
								'userData',UpdateContour(state));
% add new point								
				else,				
					cp = get(gca, 'currentPoint');
					cp = cp(1,1:2);
					lh = MakePoint(cp(1),cp(2),state.CONFIG);

% if new point is within existing points (less than half distance between nearest two points)
% then add it between those points, else append it to nearest end
					if isempty(state.ANCHORS),
						k = 1;
					else,
						d = sqrt(sum((ones(size(state.ANCHORS,1),1)*cp - state.ANCHORS).^2,2));
						[~,k] = min(d);
					end;
					if k == 1,				% add to beginning
						state.ALH = [lh , state.ALH];
						state.ANCHORS = [cp ; state.ANCHORS];
					elseif k == length(d),	% add to end
						state.ALH(end+1) = lh;
						state.ANCHORS(end+1,:) = cp;
					else,
						if d(k-1) < d(k+1), k2 = k; k = k-1; else, k2 = k+1; end;
						state.ALH = [state.ALH(1:k) , lh , state.ALH(k2:end)];
						state.ANCHORS = [state.ANCHORS(1:k,:) ; cp ; state.ANCHORS(k2:end,:)];
					end;
					set(gcbf, 'windowButtonMotionFcn',{@GetContours,'MOVE',lh}, ...
								'windowButtonUpFcn',{@GetContours,'UP',lh}, ...
								'pointer','crosshair', ...
								'userData',UpdateContour(state));
				end;

% delete all points
			case 'open',			% double-click
				delete(findobj(gca,'tag','CONTOUR'));
				state.CLH = []; state.ALH = []; state.ANCHORS = [];
				set(gcf,'userData',UpdateContour(state));

% delete existing point (ctl) or echo position (shift)
			otherwise,
				if gotPoint,
					k = find(gcbo == state.ALH);	% clicked point index
					if strcmp(mod,'extend'),	% shift (mid)
						fprintf('%d anchor points, current point located [ %d , %d ]\n', ...
								length(state.ALH), round(get(state.ALH(k),'Xdata')), round(get(state.ALH(k),'Ydata')));
					else,						% ctl (right)
						state.ANCHORS(k,:) = [];
						state.ALH(k) = [];
						delete(gcbo);
						if isempty(state.ANCHORS), delete(state.CLH); state.CLH = []; end;
						set(gcf,'userData',UpdateContour(state));
					end;
				end;
		end;
		
	
%-----------------------------------------------------------------------------
% EMIT:  export contours to workspace

	case 'EMIT',
		state = get(gcf,'userData');
		v = evalin('base',state.VNAME);						% update output variable state
		frames = cell2mat({v.FRAME});
		k = (state.CURFRAME == frames);
		v(k).XY = state.XY;
		v(k).ANCHORS = state.ANCHORS;
		[n,k] = sort(frames);
		v = v(k);	
		v(cellfun(@isempty,{v.XY})) = [];					% delete empty frames
		contours = reshape(cell2mat({v.XY}),[state.NPOINTS 2 length(v)]);
		frames = cell2mat({v.FRAME});
		if length(varargin) > 1,
			idx = varargin{2};
			[v,k] = intersect(frames,idx);
			contours = contours(:,:,k);
			frames = frames(k);
		end;
		varargout{1} = contours;
		if nargout > 1,
			varargout{2} = frames;
			if nargout > 2,		% include images
				q = get(state.IH,'cdata');
				q = zeros(size(q,1),size(q,2),length(frames),'uint8');
				for f = 1 : length(frames),
					img = uint8(mean(GetMovieFrame(state.MH,frames(f),state.FRATE),3));
					if ~isempty(state.IMGMODP), img = state.IMGMODP(img,state.IMGMODA); end;
					q(:,:,f) = img;
				end;
				varargout{3} = q;
			end;
		end;
		return;
		
		
%-----------------------------------------------------------------------------
% EXPORT:  save contours to tab delimited output file

	case 'EXPORT',
		state = get(gcbf,'userData');
		[fName,pName] = uiputfile('*.txt','Export Contours to tab-delimited text file',[state.VNAME,'.txt']);
		if fName == 0, return; end;		% cancel
		
		v = evalin('base',state.VNAME);			
		frames = cell2mat({v.FRAME});			% frames with data
		k = (state.CURFRAME == frames);
		v(k).XY = state.XY;
		v(k).ANCHORS = state.ANCHORS;
		v(k).NOTE = get(state.LH,'string');
		[n,k] = sort(frames);
		v = v(k);	
		v(cellfun(@isempty,{v.XY})) = [];					% delete empty frames
		contours = reshape(cell2mat({v.XY}),[state.NPOINTS 2 length(v)]);
		frames = cell2mat({v.FRAME});
		notes = {v.NOTE};
		sr = state.MH.FrameRate;
		fid = fopen(fName,'wt');
		if fid < 0, fprintf('Error attempting to open %s\n',fName); return; end;
		fprintf(fid,'FRAME\tTIME\tPOINT\tNOTE\tX\tY\n');
		for fi = 1 : length(frames),
			for ci = 1 : size(contours,1),
				fprintf(fid,'%d\t%f\t%s\t%d\t%f\t%f\n', frames(fi), (frames(fi)-1)/sr, notes{fi}, ci, contours(ci,1,fi), contours(ci,2,fi));
			end;
		end;
		fclose(fid);
		fprintf('wrote %s\n',fName);
		
		
%-----------------------------------------------------------------------------
% FILTER:  show image forces

	case 'FILTER',
		state = get(gcbf,'userData');
		img = ComputeImageForces(im2double(get(state.IH,'cdata')));
		img = im2uint8((img-min(img(:)))./range(img(:)));
		set(state.IH,'cdata',img);
		
		
%-----------------------------------------------------------------------------
% FLIP:  invert image

	case 'FLIP',
		if strcmp(varargin{2},'HORIZONTAL'),
			if strcmp(get(gca,'xdir'),'normal'), set(gca,'xdir','reverse'); else, set(gca,'xdir','normal'); end;
		else,
			if strcmp(get(gca,'ydir'),'normal'), set(gca,'ydir','reverse'); else, set(gca,'ydir','normal'); end;
		end;
		
		
%-----------------------------------------------------------------------------
% FRAME:  set image frame

	case 'FRAME',
		state = get(gcbf,'userData');
		v = evalin('base',state.VNAME);
		frames = cell2mat({v.FRAME});			% frames with data
		f = state.CURFRAME;						% current frame before change
		switch varargin{2},
			case 'PREV', if f > 1, f = f - 1; end;
			case 'NEXT', if f < state.NFRAMES, f = f + 1; end;
			case 'SPECIFY',
				f = GetFrame(f, state.NFRAMES);
				if isempty(f), return; end;
				if f<1, f = 1; elseif f>state.NFRAMES, f = state.NFRAMES; end;
			case {'KPREV','KNEXT'},
				[~,k] = min(abs(state.KEYFRAMES(:,1)-f));
				if strcmp(varargin{2},'KPREV'),
					if k > 1, k = k - 1; end;
				else,
					if k < size(state.KEYFRAMES,1), k = k + 1; end;
				end;
				f = state.KEYFRAMES(k,1);
			case {'DPREV','DNEXT'},
				fr = sort(frames);
				if strcmp(varargin{2},'DPREV'),
					fr = fliplr(fr(fr<f));
				else
					fr = fr(fr>f);
				end
				if isempty(fr), return; end
				f = fr(1);
		end;
		if state.CURFRAME == f, return; end;	% nothing to do
		
% update output variable state
		k = (state.CURFRAME == frames);
		v(k).XY = state.XY;
		v(k).ANCHORS = state.ANCHORS;
		v(k).NOTE = get(state.LH,'string');
		assignin('base',state.VNAME,v);

% change frame
		state.CURFRAME = f;

% get current annotation if any
		ts = '';
		k = find(f == frames);
		if ~isempty(k), ts = v(k).NOTE; end

% use Praat labels if no annotation
		if isempty(ts) && ~isempty(state.KEYFRAMES),
			if size(state.KEYFRAMES,2) == 1,
				k = find(f == state.KEYFRAMES);
				if ~isempty(k), ts = state.LABELS{k}; end
			else,
				k = find(f>=state.KEYFRAMES(:,1) & f<=state.KEYFRAMES(:,2));
				if ~isempty(k), k = k(end); ts = state.LABELS{k(1)}; end
			end
		end

% update
		set(state.IH,'cdata',GetImage(state));
		set(state.TH,'string',sprintf('%04d',f));
		set(state.LH,'string',ts);
		delete(findobj(gca,'tag','CONTOUR'));
		state.ALH = []; state.CLH = [];
		k = find(f == frames);
		if isempty(k),				% virgin frame
			if strcmpi('off',get(state.NH,'checked')),
				state.ANCHORS = [];		% don't inherit
				state.XY = [];
			end;
			vv = struct('XY',state.XY,'ANCHORS',state.ANCHORS,'FRAME',f,'NOTE',ts);
			v(end+1) = vv;			
		else,						% update from existing anchors
			state.PREVAP = state.ANCHORS;
			state.ANCHORS = v(k).ANCHORS;
			state.XY = v(k).XY;
			v(k).NOTE = ts;
		end;
		assignin('base',state.VNAME,v);

		for k = 1 : size(state.ANCHORS,1),
			state.ALH(k) = MakePoint(state.ANCHORS(k,1),state.ANCHORS(k,2),state.CONFIG);
		end;
		if ~isempty(state.XY),
			state.CLH = line(state.XY(:,1),state.XY(:,2),'color',state.CONFIG.LINECOLOR,'linewidth',state.CONFIG.LINEWIDTH,'tag','CONTOUR','hitTest','off');
			uistack(state.CLH,'bottom'); uistack(state.CLH,'up');
		end;
		
		set(gcbf,'name',sprintf('%s  [%d of %d]',state.FNAME,f,state.NFRAMES),'userData',state);
		
	
%-----------------------------------------------------------------------------
% MAP:  set colormap

	case 'MAP',
		set(get(get(gcbo,'parent'),'children'),'checked','off');
		set(gcbo,'checked','on');
		mapName = varargin{2};
		if strcmp(mapName,'Inv Gray'),
			map = 1-gray;
		else,
			map = eval(lower(get(gcbo,'label')));
		end;
		set(gcbf,'colormap',map);
		

%-----------------------------------------------------------------------------
% MOVE:  mouseMvt handler

	case 'MOVE',
		cp = get(gca, 'currentPoint');
		cp = cp(1,1:2);
		lh = varargin{2};
		if ~ishandle(lh), return; end;
		x = get(lh,'xdata');
		y = get(lh,'ydata');
		if x ~= cp(1) || y ~= cp(2), 
			set(lh, 'xdata',cp(1), 'ydata',cp(2));
			state = get(gcbf,'userData');
			k = find(lh == state.ALH);
			state.ANCHORS(k,:) = [x,y];
			set(gcbf,'userData',UpdateContour(state));
		end;


%-----------------------------------------------------------------------------
% RANGE:  adjust range

	case 'RANGE',
		if strcmp(varargin{2},'FULL'),
			set(gca,'clim',[0 255]);
			return;
		end;
		state = get(gcbf,'userData');
		if ishandle(state.RH),
			figure(state.RH);
			return;
		end;
		state.RH = imcontrast(state.IH);
		set(gcbf,'userData',state);

	
%-----------------------------------------------------------------------------
% REDISTRIBUTE:  redistribute anchor points along current contour

	case 'REDISTRIBUTE',
		state = get(gcbf,'userData');
		nAnchors = size(state.ANCHORS,1);
		if nAnchors < 3, return; end;	% need at least three anchors
		state.PREVAP = state.ANCHORS;
		k = [0 ; cumsum(sqrt(sum(diff(state.XY).^2,2)))];
		state.ANCHORS = interp1(k,state.XY,linspace(0,k(end),nAnchors),'linear');
		for k = 1 : nAnchors,
			set(state.ALH(k),'xdata',state.ANCHORS(k,1),'ydata',state.ANCHORS(k,2));
		end;
		set(gcbf,'userData',UpdateContour(state));
		

%-----------------------------------------------------------------------------
% RESET:  reset to original image

	case 'RESET',
		state = get(gcbf,'userData');
		state.USEAVG = 0;
		set(state.IH,'cdata',GetImage(state));
		set(gcbf,'userData',state);
		set(gca,'xdir','normal','ydir','reverse','clim',state.CLIM);
		
		
%-----------------------------------------------------------------------------
% TRACK:  apply automatic tracker

	case 'TRACK',
		state = get(gcbf,'userData');
		nAnchors = size(state.ANCHORS,1);
		if nAnchors < 4, return; end;
		state.PREVAP = state.ANCHORS;
		set(gcbf,'pointer','watch'); drawnow;
		img = im2double(get(state.IH,'cdata'));
		if isempty(state.TRACKERA),
			xy = state.TRACKERP(img, state.ANCHORS, state.NPOINTS);
		else,
			xy = state.TRACKERP(img, state.ANCHORS, state.NPOINTS, state.TRACKERA{:});
		end;
		set(gcbf,'pointer','arrow');
		state.XY = xy;
		set(state.CLH,'xdata',xy(:,1),'ydata',xy(:,2));
		k = [0 ; cumsum(sqrt(sum(diff(xy).^2,2)))];
		state.ANCHORS = interp1(k,xy,linspace(0,k(end),nAnchors),'linear');
		for k = 1 : nAnchors,				% redistribute anchors
			set(state.ALH(k),'xdata',state.ANCHORS(k,1),'ydata',state.ANCHORS(k,2));
		end;
		set(gcbf,'userData',state);


%-----------------------------------------------------------------------------
% UNDO:  unwind last change

	case 'UNDO',
		state = get(gcbf,'userData');
		curA = state.ANCHORS;
		state.ANCHORS = state.PREVAP;
		state.PREVAP = curA;
		delete(findobj(gca,'tag','CONTOUR'));
		state.ALH = []; state.CLH = [];
		set(gcbf,'userData',UpdateContour(state));
		
		
%-----------------------------------------------------------------------------
% UP:  mouseUp handler

	case 'UP',
		set(gcbf, 'windowButtonMotionFcn','', 'windowButtonUpFcn','', 'pointer','arrow');
		cp = get(gca, 'currentPoint');
		cp = cp(1,1:2);
		lh = varargin{2};
		if ~ishandle(lh), return; end;
		x = get(lh,'xdata');
		y = get(lh,'ydata');
		set(lh, 'xdata',cp(1), 'ydata',cp(2));
		state = get(gcbf,'userData');
		k = find(lh == state.ALH);
		state.ANCHORS(k,:) = [x,y];
		set(gcbf,'userData',UpdateContour(state));
		
		
%-----------------------------------------------------------------------------
% INITialize

	otherwise,
		fh = findobj('tag','GETCONTOURS');
		if isempty(fh),
			Initialize(varargin{:});
		else,
			figure(fh);				% allow only one active GETCONTOURS object
		end;
		if isunix, [s,r] = unix('osascript -e ''tell application "MATLAB" to activate'''); end;
end;


%=============================================================================
% COMPUTEIMAGEFORCES  - Canny edge filter based on Gaussian image derivatives
% cf. Image Toolbox edge

function img = ComputeImageForces(img)

% tweaks
sigma = 4;
GaussianDieOff = .0001;

% design the filters - a gaussian and its derivative
pw = 1:30; 			% supported widths
ssq = sigma^2;
width = find(exp(-(pw.*pw)/(2*ssq))>GaussianDieOff,1,'last');
if isempty(width), width = 1; end

t = (-width:width);
gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq);     % the Gaussian 1D filter

% find the directional derivative of 2D Gaussian
[x,y] = meshgrid(-width:width,-width:width);
dgau2D = -x.*exp(-(x.*x+y.*y)/(2*ssq))/(pi*ssq);

% convolve the filters with the image in each direction
img = imfilter(img,gau,'conv','replicate');   		% run the filter across rows
img = imfilter(img,gau','conv','replicate'); 		% and then across columns

% apply directional derivatives
ax = imfilter(img, dgau2D, 'conv','replicate');
ay = imfilter(img, dgau2D', 'conv','replicate');

% normalize
img = sqrt((ax.*ax) + (ay.*ay));
if max(img(:)) > 0, img = img ./ max(img(:)); end
img = 1 - img;


%=============================================================================
% GETAVG  - get averaging options
%
% returns [] on cancel

function [enabled,win] = GetAvg(state)

width = 220;
height = 120;
fPos = get(gcf, 'position');
pos = [fPos(1)+(fPos(3)-width)/2 , fPos(2)+(fPos(4)-height)/2 , width , height];

cfg = dialog('name', 'Frames to average:', ...
	'menubar', 'none', ...
	'position', pos, ...
	'keyPressFcn', 'set(gcbf,''UserData'',1);uiresume', ...
	'UserData', 0);

% enabled button
h = height - 30;
ebH = uicontrol(cfg, ...
	'Position', [40 h 150 25], ...
	'Style', 'checkbox', ...
	'Value', state.USEAVG, ...
	'String', 'Averaging enabled');
	
% win field
h = h - 30;
winH = uicontrol(cfg, ...
	'Position', [40 h 40 25], ...
	'Style', 'edit', ...
	'String', num2str(state.AVG));
uicontrol(cfg, ...
	'Position', [82 h-4 100 25], ...
	'Style', 'text', ...
	'String', '# frames pre/post');

% OK, cancel buttons
uicontrol(cfg, ...		% buttons
	'Position',[width/2-70 15 60 25], ...
	'String','OK', ...
	'Callback','set(gcbf,''UserData'',1);uiresume');
uicontrol(cfg, ...
	'Position',[width/2+10 15 60 25], ...
	'String','Cancel', ...
	'Callback','set(gcbf,''UserData'',0);uiresume');

% wait for input
uiwait(cfg);
if ishandle(cfg) && get(cfg, 'UserData'),
	enabled = get(ebH,'value');
	win = str2num(get(winH,'string'));
	if isnan(win) || win<0, win = 0; end;
else,
	enabled = []; win = [];
end;
if ishandle(cfg), delete(cfg); end;


%=============================================================================
% GETFRAME  - get frame number
%
% returns [] on cancel

function frame = GetFrame(curFrame, nFrames)

width = 200;
height = 100;
fPos = get(gcf, 'position');
pos = [fPos(1)+(fPos(3)-width)/2 , fPos(2)+(fPos(4)-height)/2 , width , height];

cfg = dialog('name', 'Go to frame...', ...
	'menubar', 'none', ...
	'position', pos, ...
	'keyPressFcn', 'set(gcbf,''UserData'',1);uiresume', ...
	'UserData', 0);

% entry field
eh = uicontrol(cfg, ...
	'position', [20 60 80 25], ...
	'style', 'edit', ...
	'horizontalAlignment', 'left', ...
	'string', num2str(curFrame));
uicontrol(cfg, ...
	'position', [100 60 60 20], ...
	'style', 'text', ...
	'string', sprintf('of %d', nFrames));

% OK, cancel buttons
uicontrol(cfg, ...		% buttons
	'Position',[width/2-70 15 60 25], ...
	'String','OK', ...
	'Callback','set(gcbf,''UserData'',1);uiresume');
uicontrol(cfg, ...
	'Position',[width/2+10 15 60 25], ...
	'String','Cancel', ...
	'Callback','set(gcbf,''UserData'',0);uiresume');

% wait for input
uiwait(cfg);
if ishandle(cfg) && get(cfg, 'UserData'),
	frame = str2num(get(eh, 'string'));
else,
	frame = [];
end;
if ishandle(cfg), delete(cfg); end;


%=============================================================================
% GETIMAGE  - load (possibly averaged) image

function img = GetImage(state)

fail = 0;
if state.USEAVG,
	f = state.CURFRAME;
	pre = f - state.AVG;
	post = f + state.AVG;
	if pre < 1, pre = 1; end;
	if post > state.NFRAMES, post = state.NFRAMES; end;
	try,
		img1 = uint8(mean(GetMovieFrame(state.MH,pre,state.FRATE),3));
		if ~isempty(state.IMGMODP), img1 = state.IMGMODP(img1,state.IMGMODA{:}); end;
		img = zeros(size(img1,1),size(img1,2),post-pre+1,'uint8');
		img(:,:,1) = img1;
		for k = 1 : post-pre,
			img1 = uint8(mean(GetMovieFrame(state.MH,k+pre,state.FRATE),3));
			if ~isempty(state.IMGMODP), img1 = state.IMGMODP(img1,state.IMGMODA{:}); end;
			img(:,:,k) = img1;
		end;
		img = uint8(mean(img,3));
	catch
		fail = 1;
	end
else,
	try,
		img = uint8(mean(GetMovieFrame(state.MH,state.CURFRAME,state.FRATE),3));
		if ~isempty(state.IMGMODP), img = state.IMGMODP(img,state.IMGMODA{:}); end;
	catch
		fail = 1;
	end
end;
if fail,
	img = zeros(get(state.MH,'Height'),get(state.MH,'Width'),'uint8');
	fprintf('unable to load frame %d from %s\n', state.CURFRAME, get(state.MH,'name'));
end


%=============================================================================
% GETMOVIEFRAME  - get frame from open video object

function img = GetMovieFrame(mh, frame, fRate)

if verLessThan('matlab','8.5.0'),
	try, img = read(mh, frame); catch, error('unable to read frame %d from specified movie handle',frame); end;
else,
	mh.CurrentTime = (frame-1)/fRate;
	try, img = readFrame(mh); catch, error('unable to read frame %d from specified movie handle',frame); end;
end;


%=============================================================================
% GETTEXT  - get frame annotation

function ts = GetText(ts,fn)

width = 250;
height = 120;
fPos = get(gcf, 'position');
pos = [fPos(1)+(fPos(3)-width)/2 , fPos(2)+(fPos(4)-height)/2 , width , height];

cfg = dialog('Name', 'Annotate', ...
	'menubar', 'none', ...
	'position', pos, ...
	'keyPressFcn', 'set(gcbf,''UserData'',1);uiresume', ...
	'UserData', 0);

% text field
uicontrol(cfg, ...
	'position', [20 80 width-40 20], ...
	'style', 'text', ...
	'horizontalAlignment', 'left', ...
	'string', sprintf('Enter annotation for frame %d',fn));
eh = uicontrol(cfg, ...
	'position', [20 60 width-40 25], ...
	'style', 'edit', ...
	'horizontalAlignment', 'left', ...
	'callback', 'set(gcbf,''UserData'',1);uiresume', ...
	'string', ts);

% OK, cancel buttons
uicontrol(cfg, ...		% buttons
	'Position',[width/2-70 15 60 25], ...
	'String','OK', ...
	'Callback','set(gcbf,''UserData'',1);uiresume');
uicontrol(cfg, ...
	'Position',[width/2+10 15 60 25], ...
	'String','Cancel', ...
	'Callback','uiresume');

% wait for input
uiwait(cfg);
if get(cfg, 'UserData'),
	ts = get(eh, 'string');
else,
	ts = [];
end;
delete(cfg);


%===============================================================================
% INITIALIZE

function Initialize(fName, varargin)

% defaults
inherit = 0;
keyFrames = [];
nPoints = 100;
saveImg = 0;
textgrid = [];
ImageMod = [];
Tracker = [];
cfg = [];
defCfg = struct('DOTSIZE',15,'DOTCOLOR','r','LINEWIDTH',1,'LINECOLOR','y');
if ismac, defCfg.DOTSIZE = 25; end;

[~,vName] = fileparts(fName);
for ai = 2 : 2 : length(varargin),
	switch upper(varargin{ai-1}),
		case 'CONFIG', cfg = varargin{ai}; 
		case 'IMGMOD', ImageMod = varargin{ai}; 
		case 'KEYFRAMES', keyFrames = varargin{ai}; keyFrames = keyFrames(:);
		case 'NPOINTS', nPoints = varargin{ai};
		case 'SAVEIMG', saveImg = strcmpi(varargin{ai}(1),'T');
		case 'TEXTGRID',textgrid = varargin{ai};
		case 'TRACKER', Tracker = varargin{ai};
		case 'VNAME',vName = varargin{ai};
		otherwise, error('unsupported option (%s)', varargin{ai-1});
	end;
end;

% override display defaults
if ~isempty(cfg),
	fnc = upper(fieldnames(cfg));
	fnd = fieldnames(defCfg);
	for fi = 1 : length(fnc),
		k = find(strcmp(fnc{fi},fnd));
		if ~isempty(k), defCfg.(fnd{k}) = cfg.(fnc{fi}); end;
	end;
end;
cfg = defCfg;

% reload data from existing mat file
if exist([vName,'.mat']) == 2,
	v = load([vName,'.mat']);
	v = v.(vName);
	if evalin('base',sprintf('exist(''%s'',''var'')',vName)),
		v = evalin('base',vName);
		assignin('base',[vName,'_old'],v);		% any previous variable renamed to <NAME>_old
	end
	if ~isfield(v,'NOTE'), v(1).NOTE = ''; end	% bw compatibility
	assignin('base',vName,v);
end
	
if evalin('base',sprintf('exist(''%s'',''var'')',vName)),
	if isunix, [s,r] = unix('osascript -e ''tell application "MATLAB" to activate'''); end;
	if strcmp(questdlg(sprintf('%s exists; overwrite its values?\n(%s_old backup will be created)',vName,vName), ... 
			'Verify...', 'Yes', 'No', 'Yes'), 'No'), return; end;
end;

% open the movie
fNameExt = fName;
[p,fName,e] = fileparts(fName);
if isempty(e), fNameExt = fullfile(p,[fName,'.avi']); end;
try, 
	mh = VideoReader(fNameExt); 
catch, 
	error('unable to open movie file %s', fNameExt); 
end;
if verLessThan('matlab','8.5.0'),
	nFrames = get(mh,'NumberOfFrames');
	frameRate = get(mh,'FrameRate');
else,
	frameRate = mh.FrameRate;
	nFrames = floor(mh.Duration*frameRate);
end;

% parse the TextGrid
keyFrames = []; labs = [];
if ~isempty(textgrid), [keyFrames,labs] = ParseTextGrid(textgrid,frameRate); end;
if isempty(keyFrames), 
	frame = 1; ts = ''; enableS = 'off';
else, 
	frame = keyFrames(1); ts = labs{1}; enableS = 'on';
end

% load first image
try,
	img = uint8(mean(GetMovieFrame(mh,frame,frameRate),3));
	if ~isempty(ImageMod), 
		p = ImageMod{1};
		a = ImageMod(2:end);
		img = p(img,a{:}); 
	end;
catch,
	error('unable to load frame %d in %s (%d frames available)', frame, fNameExt, nFrames);
end;
[rows,cols] = size(img);

% display image
pos = get(0, 'defaultfigureposition');
fh = figure('units','pixels', 'resize','off', 'tag','GETCONTOURS', 'position',[pos(1) , pos(2)+pos(4)-rows , cols, rows]);
set(fh, 'colormap', gray(256));
ih = imagesc(img);
clim = get(gca,'clim');
set(gca, 'position',[0 0 1 1], 'xtick',[],'ytick',[]);

% display frame number
th = text(50,80,sprintf('%04d',frame),'fontsize',18,'color','w');

% display frame label
lh = text(50,100,ts, 'fontsize',18,'color','w','interpreter','none');

% init menu
cbs = 'if strcmp(''on'',get(gcbo,''checked'')),set(gcbo,''checked'',''off'');else,set(gcbo,''checked'',''on'');end';
menu = uimenu(fh, 'label','GETCONTOURS');
uimenu(menu,'label','About...', ...
			'callback',{@GetContours,'ABOUT'});
uimenu(menu,'label','Export Contours...', ...
			'callback',{@GetContours,'EXPORT'});

uimenu(menu,'label','Previous Frame', ...
			'separator','on', ...
			'accelerator', 'B', ...
			'callback',{@GetContours,'FRAME','PREV'});
uimenu(menu,'label','Next Frame', ...
			'accelerator', 'F', ...
			'callback',{@GetContours,'FRAME','NEXT'});
uimenu(menu,'label','Specify Frame...', ...
			'callback',{@GetContours,'FRAME','SPECIFY'});

uimenu(menu,'label','Previous Key Frame', ...
			'separator','on', ...
			'enable', enableS, ...
			'accelerator', '1', ...
			'callback',{@GetContours,'FRAME','KPREV'});
uimenu(menu,'label','Next Key Frame', ...
			'enable', enableS, ...
			'accelerator', '2', ...
			'callback',{@GetContours,'FRAME','KNEXT'});
uimenu(menu,'label','Prev Frame With Data', ...
			'accelerator', '3', ...
			'callback',{@GetContours,'FRAME','DPREV'});
uimenu(menu,'label','Next Frame With Data', ...
			'accelerator', '4', ...
			'callback',{@GetContours,'FRAME','DNEXT'});
uimenu(menu,'label','Annotate Frame', ...
			'accelerator', 'A', ...
			'callback',{@GetContours,'ANNOTATE'});

nh = uimenu(menu,'label','Inherit Anchors', ...
			'separator','on', ...
			'accelerator', 'I', ...
			'callback',cbs);
uimenu(menu,'label','Redistribute Anchors', ...
			'accelerator','R', ...
			'callback',{@GetContours,'REDISTRIBUTE',1});
if isempty(Tracker), s = 'off'; else, s = 'on'; end;
uimenu(menu,'label','Apply Tracking', ...
			'accelerator','T', ...
			'enable',s, ...
			'callback',{@GetContours,'TRACK'});
uimenu(menu,'label','Undo Last Change', ...
			'accelerator','Z',...
			'callback',{@GetContours,'UNDO',1});
			
uimenu(menu,'label','Filter', ...
			'separator','on', ...
			'accelerator', 'L', ...
			'callback',{@GetContours,'FILTER'});
uimenu(menu,'label','Original', ...
			'accelerator', 'G', ...
			'callback',{@GetContours,'RESET'});
uimenu(menu,'label','Averaging...', ...
			'callback',{@GetContours,'AVERAGING'});

h = uimenu(menu,'label','Contrast','separator','on');
uimenu(h,'label','Full Range','callback',{@GetContours,'RANGE','FULL'});
uimenu(h,'label','Adjust Range','callback',{@GetContours,'RANGE','ADJUST'});

h = uimenu(menu,'label','Colormap');
uimenu(h,'label','Gray','checked','on','callback',{@GetContours,'MAP','Gray'});
uimenu(h,'label','Inv Gray','callback',{@GetContours,'MAP','Inv Gray'});
cMaps = {'Autumn','Bone','Colorcube','Cool','Copper','Flag','Hot','HSV','Jet','Lines','Pink','Prism','Spring','Summer','Winter'};
for k = 1 : length(cMaps),
	if k == 1, sep = 'on'; else, sep = 'off'; end;
	uimenu(h,'label',cMaps{k},'separator',sep,'callback',{@GetContours,'MAP',cMaps{k}});
end;

h = uimenu(menu,'label','Flip');
uimenu(h,'label','Horizontal', ...
			'callback',{@GetContours,'FLIP','HORIZONTAL'});
uimenu(h,'label','Vertical', ...
			'callback',{@GetContours,'FLIP','VERTICAL'});
		
% init internal state
state = struct('IH', ih, ...				% image handle
				'MH', mh, ...				% movie handle
				'TH', th, ...				% frame text handle
				'LH', lh, ...				% annotation text handle
				'NH', nh, ...				% inherit anchors menu handle
				'CLIM', clim, ...			% original contrast limits
				'CONFIG', cfg, ...			% display configuration
				'IMGMODP', [], ...			% image preprocessor 
				'IMGMODA', [], ...			% and its arguments
				'TRACKERP', [], ...			% automatic tracking procedure 
				'TRACKERA', [], ...			% and its arguments
				'CURFRAME', frame, ...		% currently displayed frame
				'KEYFRAMES', keyFrames, ...	% keyframe list
				'NFRAMES', nFrames, ...		% number of available frames
				'FRATE', frameRate, ...		% movie frame rate
				'LABELS', [], ...			% Praat keyframe labels
				'NPOINTS', nPoints, ...		% number of contour points
				'ANCHORS', [], ...			% current anchor points
				'ALH', [], ...				% and their line handles
				'PREVAP', [], ...			% previous anchor points
				'XY', [], ...				% current contour points
				'CLH', [], ...				% and their line handle
				'USEAVG', 0, ...			% frame averaging enabled
				'SAVEIMG', saveImg, ...		% save images
				'AVG', 0, ...				% pre/post # frames to average
				'RH', [], ...				% contrast adjustment handle
				'FNAME', fName, ...			% movie name
				'VNAME', vName);			% emit array name
state.LABELS = labs;
if ~isempty(ImageMod),
	if iscell(ImageMod), state.IMGMODP = ImageMod{1}; else, state.IMGMODP = ImageMod; end;
	if length(ImageMod) > 1, state.IMGMODA = ImageMod(2:end); end;
end;
if ~isempty(Tracker),
	if iscell(Tracker), state.TRACKERP = Tracker{1}; else, state.TRACKERP = Tracker; end;
	if length(Tracker) > 1, state.TRACKERA = Tracker(2:end); end;
end;

set(fh,'name',sprintf('%s  [%d of %d]',fName,frame,nFrames), ...
		'numberTitle','off', ...
		'closeRequestFcn',{@GetContours,'CLOSE'}, ...
		'tag','GETCONTOURS', ...
		'userData',state);
set(ih,'buttonDownFcn', {@GetContours,'DOWN'});		% trap buttonDown

% create or link to output variable
v = []; anchors = [];
if evalin('base',sprintf('exist(''%s'',''var'')',vName)),
	v = evalin('base',vName);
	if isstruct(v) && isfield(v,'ANCHORS'),
		if isfield(v,'TIME'), v = rmfield(v,{'TIME','IMAGE'}); end;
		k = find(frame == cell2mat({v.FRAME}));
		if isempty(k),
			v(end+1) = struct('XY',[],'ANCHORS',[],'FRAME',frame,'NOTE','');
		else,
			anchors = v(k).ANCHORS;
		end;
	end;	
end;

% initialize contour from passed-in anchor points
if ~isempty(anchors),
	state.ANCHORS = anchors;
	set(fh,'userData',UpdateContour(state)); 
end;
if isempty(v),
	v = struct('XY',state.XY,'ANCHORS',anchors,'FRAME',frame,'NOTE','');
end;
assignin('base',vName,v);


%===============================================================================
% MAKEPOINT  - create anchor point

function lh = MakePoint(x,y,cfg)

if ismac, ms = 25; else, ms = 15; end;
lh = line(x,y,'marker','.','markerSize',cfg.DOTSIZE,'color',cfg.DOTCOLOR,'tag','CONTOUR','buttonDownFcn',{@GetContours,'DOWN','POINT'});


%===============================================================================
% PARSETEXTGRID  - load specified tier from Praat TextGrid

function [frames,labs] = ParseTextGrid(fName, sr)

if iscell(fName),
	tName = fName{2};
	fName = fName{1};
else,
	tName = [];										% use first tier found
end
[segs,labs] = ReadPraatTier(fName, tName);
if size(segs,2) > 1,								% interval tier
	k = find(cellfun(@isempty,labs));
	segs(k,:) = [];
	labs(k) = [];
	if isempty(segs),
		[p,f,e] = fileparts(fName);
		if isempty(e), fName = fullfile(p,[f,'.TextGrid']); end
		fprintf('Warning:  no labeled intervals found in %s (so no key frames available)\n', fName);
		frames = [];
		return;
	end
end
frames = floor(segs*sr) + 1;						% convert to 1-based frame numbers
k = find(diff(frames,1,2) < 1);
if ~isempty(k),
	fprintf('warning:  %d intervals were not included as key frames (less than one frame in length)\n', length(k));
	frames(k,:) = [];
	segs(k) = [];
end


%===============================================================================
% SAVEVAR  - save variable as VNAME to specified MAT file

function SaveVar(var, fName)

vName = fName;
vs = vName;
eval([vName '=var;save ' fName ' ' vs]);
fprintf('wrote %s.mat\n', fName);


%===============================================================================
% UPDATECONTOUR  - update displayed contour

function state = UpdateContour(state)

% create anchor points
if isempty(state.ALH) && ~isempty(state.ANCHORS),
	for k = 1 : size(state.ANCHORS,1),
		state.ALH(k) = MakePoint(state.ANCHORS(k,1),state.ANCHORS(k,2),state.CONFIG);
	end;
end;

% create contour line
if isempty(state.CLH) && ~isempty(state.ANCHORS),
	state.CLH = line(state.ANCHORS(1,1),state.ANCHORS(1,2),'color',state.CONFIG.LINECOLOR,'linewidth',state.CONFIG.LINEWIDTH,'tag','CONTOUR','hitTest','off');
	if ~verLessThan('matlab','8.5.0'), state.CLH.PickableParts = 'none'; end;
	uistack(state.CLH,'bottom'); uistack(state.CLH,'up');
end;

% update contour line from anchor points
switch size(state.ANCHORS,1),
	case 0,
		state.XY = [];
	case 1,
		state.XY = ones(state.NPOINTS,1) * state.ANCHORS;
		set(state.CLH,'xdata',state.ANCHORS(1,1),'ydata',state.ANCHORS(1,2));
	case 2,
		k = [0 ; sqrt(sum(diff(state.ANCHORS).^2,2))];
		state.XY = interp1(k,state.ANCHORS,linspace(0,k(end),state.NPOINTS),'linear');
		set(state.CLH,'xdata',state.XY(:,1),'ydata',state.XY(:,2));
	otherwise,
		k = [0 ; cumsum(sqrt(sum(diff(state.ANCHORS).^2,2)))];
		state.XY = interp1(k,state.ANCHORS,linspace(0,k(end),state.NPOINTS),'pchip');
		set(state.CLH,'xdata',state.XY(:,1),'ydata',state.XY(:,2));
end;
drawnow;

