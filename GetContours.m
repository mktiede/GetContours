function varargout = GetContours(varargin)
%GETCONTOURS  - extract contours from US movie frames
%
%	usage:  GetContours(fName, ...)     % to initialize
%
% displays current FRAME from movie or DICOM FNAME
% click along contour to place anchor points
% click and drag on a point to reposition it
% control (right) click on a point deletes it
% shift (mid) click on a point reports anchor point info
% undo reverts to previous anchor positions
%
% use menu entries to adjust anchor points and image display
%
% supported optional 'NAME',VALUE parameter pairs:
%   ANCHORS  - seed initial frame with these anchor points [nAnchors x X,Y]
%   AUDIO    - display available audio track if true (default); click on it to set frame
%   CLEAN    - ignore any previous data if true (logical 1)
%   CONFIG   - a struct specifying display properties to modify; valid fields are
%                DOTSIZE, DOTCOLOR, LINEWIDTH, LINECOLOR
%   CROP     - crop displayed image to cropped rectangle [Xmin Ymin width height]
%   IMGMOD   - procedure applied to every image before loading (e.g. @histeq)
%   FRAME    - initial frame to display (default is first keyframe, first if none)
%   KEYFRAME - key frame list (default none)
%   MPP      - mm/pixel ratio (default none)
%   NPOINTS  - number of points / contour (default 100)
%   ORIGIN   - probe origin [X,Y] in ULC pixel coordinates (default none)
%   RESIZE   - resize displayed image by specified factor (applies after any cropping)
%   TEXTGRID - Praat TextGrid and tier to parse for key frames (see examples)
%               note that only labeled intervals are used but all point labels used
%   TRACKER  - automatic contour-fitting procedure (e.g. @gct_SLURP)
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
%   IMAGE    - grayscale image at current frame (if INCLUDE IMAGES parameter enabled)
%
% additional GetContours windows may be opened, but only the first permits editing
%
% contours, associated frame numbers and images may be exported to the workspace 
% before closing the GetContours window using
%   [contours,frames,images] = GetContours('EMIT');   % [nRows x nCols x nFrames] images
%   contour = GetContours('EMIT',FRAME);   % emit only specified frame(s)
% coordinates are in ULC-origin pixel units [NPOINTS x X,Y]; if origin and mm/pixel values 
% available two additional columns give origin-centered mm values
%
% to extract contours from VNAME use 
%   contours = reshape(cell2mat({VNAME.XY}),[NPOINTS 2 length(VNAME)]);
%
% Example:  get contours and associated frames from VNAME 'foo'
%   GetContours('movie.avi', 'VNAME','foo');
%   [contours,frames] = GetContours('EMIT');  % [nPoints x X,Y x nFrames]
% or equivalently
%   contours = reshape(cell2mat({foo.XY}),[nPoints 2 length(foo)]);
%   frames = cell2mat({foo.FRAMES})
% contours may also be EXPORTed (in long format) to a text file with GetContours, 
% and subsequently converted to wide or EdgeTrak format using the RESHAPECONTOURS procedure
%
% Example:  specify key frames from point or interval tier "frame" in "foo.TextGrid"
%   GetContours('movie.avi', 'TEXTGRID',{'foo','frame'})
% if no tier name specified loads from the first tier found
%   GetContours('movie.avi', 'TEXTGRID','foo')
% specify key frames directly, ignore any previous data
%   GetContours('movie.avi', 'KEYFRAMES',[23:47 123 247:319], 'CLEAN',true)
% 
% Example:  crop image, use thick green lines
%   cfg = struct('LINEWIDTH',2, 'LINECOLOR','g');
%   GetContours('movie.avi', 'CROP',[100 43 670 503], 'CONFIG',cfg);
% 
% Example:  resize image by 75%, use SLURP tracker
%   GetContours('movie.avi', 'RESIZE',.75, 'TRACKER',@gct_SLURP);

% Copyright (C) 2015-20 mark tiede <tiede@haskins.yale.edu>
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
% mkt 12/15 v0.8 UltraFest2015
% mkt 08/16 v0.9 bug fixes, better tracking support
% mkt 10/16 v1.0 release, fix inherit anchors into empty frame
% mkt 11/16 v1.1 fix explicit keyframes issues
% mkt 02/17 v1.2 fix point addition order
% mkt 06/18 v2.0 support sequence processing. tracking plugins, mpp detection
% mkt 07/18 v2.1 fix TextGrid issues
% mkt 10/18 v2.2 fix initialization overwrite bug
% mkt 11/19 v2.3 mods for internal improvements
% mkt 03/20 v2.4 bug fixes, scroller, DICOM, SLURP support
% mkt 08/20 v2.5 minor bug fixes, COEF support, add gct_snake (UltraFest IX release)
% mkt 08/20 v2.6 Fourier coef shape fitting support
% mkt 08/20 v2.7 support preseeded ANCHORS
% mkt 09/20 v2.8 fix DICOM close
% mkt 09/20 v3.0 support audio panel
% mkt 09/20 v3.1 support draw mode, multiple panels
% mkt 09/20 v3.2 support info, frame differencing, anchor deletion issue
% mkt 10/20 v3.3 bug fixes
% mkt 10/20 v3.4 test all files as DICOM first, fix cfg bug

% STATE (gcf userData) defines internal state for currently displayed frame
% VNAME (defined in base ws) defines values for each visited frame
% when switching frames VNAME is updated from STATE and vice versa

persistent PLAYERH

GCver = 'v3.4';	% current version

if nargin < 1,
	eval('help GetContours');
	return;
end;

% strip src,evt on callbacks
if ~ischar(varargin{1}), 
	varargin(1:2) = []; 
	if get(gca,'userdata') && ~strcmp(varargin{1},'CYCLE'), return; end;		% ignore callbacks while cycling
end;

% branch by action
switch upper(varargin{1}),

%-----------------------------------------------------------------------------
% ABOUT:  display version and brief help

	case 'ABOUT',
		state = get(gcbf,'userData');
		blurb = sprintf('%s\n\n%s\n%s\n%s\n%s\n\n%s', ...
			'Extract contours from US movie frames:', ...
			'click along contour to place anchor points', ...
			'click and drag on a point to reposition it', ...
			'control-click on a point deletes it', ...
			'shift-click on a point reports its position', ...
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
			set(gcbf,'userData',state);
		end;
		set(state.IH,'cdata',GetImage(state));
		

%-----------------------------------------------------------------------------
% CLOSE:  shutdown handler

	case 'CLOSE',
		state = get(gcbf,'userData');
		if ishandle(state.RH), delete(state.RH); end;		% close contrast adjustment (if open)
		delete(gcbf);										% close window
		if state.LOCK, return; end;

		v = evalin('base',state.VNAME);						% update output variable state
		frames = cell2mat({v.FRAME});
		k = (state.CURFRAME == frames);
		v(k).XY = state.XY;
		v(k).ANCHORS = state.ANCHORS;
		v(k).TRKRES = state.TRKRES;
		[~,k] = sort(frames);
		v = v(k);											% impose sequential order
		v(cellfun(@isempty,{v.XY})) = [];					% delete empty frames

% add additional annotation fields		
		if ischar(state.MH),
			info = dicominfo(state.MH);						% DICOM
			sr = info.CineRate;
		else,
			sr = state.MH.FrameRate;
		end;
		for vi = 1 : length(v),
			v(vi).TIME = (v(vi).FRAME-1)/sr;				% frame offset in secs from start of movie
			if state.PARAMS.SAVEIMG,
				try,
					v(vi).IMAGE = GetMovieFrame(state.MH,v(vi).FRAME,state.CROP,state.RESIZE,state.IMGMODP,state.IMGMODA{:});	% grayscale image at frame
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
% COEF:  map shape to Fourier coefficients (requires defined contour, origin & MPP)

	case 'COEF',
		state = get(gcbf,'userData');
		nc = 3;		% # coefficients
		nAnchors = size(state.ANCHORS,1);
		np = state.NPOINTS;
		gLen = 5;
		if nAnchors < 3, 
			fprintf('Coefficient estimation requires at least three anchor points\n');
			return;
		end;
		xy = state.XY;
		k = [0 ; cumsum(sqrt(sum(diff(xy).^2,2)))];		
		xy = interp1(k, xy, linspace(0, k(end), np)', 'pchip');
		[h,w] = size(get(state.IH,'cdata'));
		xy(:,1) = xy(:,1) - state.ORIGIN(1);
		xy(:,2) = state.ORIGIN(2) - xy(:,2);
		xy = (xy * state.MPP) / (.8 * 10);	% convert to cm:  distance along Liljencrants tongue with schwa is 14.4 "cm", 14.4/18 = .8
		xy(:,2) = xy(:,2) - max(xy(:,2)) + 3;	% default DC = 3 cm
		[glx,gly] = MakeGrid(np, gLen, .1, 0);
		if xy(1,1) < 0, xy = flipud(xy); end;	% want first point to be anterior
		xy0 = xy;

		if xy(1,2) > 0,				% extrapolate to last polar gridline
			d = hypot(xy(1,1),xy(1,2));
			xy = [[d,0] ; xy];
		end;
		if xy(end,2) > gly(2,end),	% extrapolate to last pharyngeal gridline
			xy = [xy ; [xy(end,1),gly(2,end)]];
		end;
		k = find(xy(:,2) < gly(2,end));
		xy(k,:) = [];
		d = ShapeToDist(xy, glx, gly);			% distance function (lips to larynx)
		[DC,C,S,mag,phi] = DistToCoef(d, nc);	% fitted coefficients
		
		fprintf('\nFourier Coefficients describing current contour (Frame %d):',state.CURFRAME);
		fprintf('\n  DC: %5.2f',DC);
		fprintf('\n Cos:'); fprintf(' %5.2f',C);
		fprintf('\n Sin:'); fprintf(' %5.2f',S);
		fprintf('\n Mag:'); fprintf(' %5.2f',mag); fprintf('  (constriction degree)');
		fprintf('\n Phi:'); fprintf(' %5.1f',phi*180/pi); fprintf('  (constriction location, in degrees)\n');
	
% adjust semi-polar gridlines
		idx = find(diff(gly),1,'last');		% index of first pharyngeal gridline
		th = cart2pol(glx(1,1:idx),gly(1,1:idx));
		[x,y] = pol2cart(th,d(1:idx));
		[glx(2,1:idx),gly(2,1:idx)] = pol2cart(th,.5);

% adjust pharyngeal gridlines
		N = size(glx,2);						% # grid lines
		y(idx:N) = gly(1,idx:N);
		x(idx:N) = d(idx:N);
		xy = [x(:),y(:)];
		flipped = diff(glx(:,end)) > 0;
		if flipped,			% faces right
			xy(idx:end,1) = -xy(idx:end,1); 
			glx(2,idx:end) = -.5;
		end;
		
% captions
		s{1} = sprintf(' DC = %5.2f', DC);
		s{2} = ['Cos =', sprintf('%6.2f',C)];
		s{3} = ['Sin =', sprintf('%6.2f',S)];
		s{4} = ['Mag =', sprintf('%6.2f',DC+sqrt(C.^2 + S.^2))];
		s{5} = ['Phi =', sprintf('%6.2f',phi*180/pi)];

		[xy,d] = CoefToShape(C, S, DC, glx, gly);

		figure('position',[35 571 560 420],'name',sprintf('Frame %04d', state.CURFRAME));
		axes('position',[.1 .75 .82 .2])
		h = plot(d,'linewidth',2);
		hh = line(get(gca,'xlim'),[DC;DC],'color','k');
		cs = {'DC'}; for k = 1 : length(h), cs{end+1} = sprintf('Coef %d',k); end; 
		legend([hh;h],cs{:});
		ylabel('Magnitude (cm)'); grid on;
		x = round([0:8]*N/8);
		for k = 1 : length(x), xs{k} = sprintf('%d',(k-1)*45-180); end;
		set(gca,'ylim',[0 gLen], 'xlim',[1 N], 'xtick',x, 'xticklabel',xs);
		text(1,gLen+.1,'Lips','fontSize',11,'verticalAlignment','bottom');
		text(N,gLen+.1,'Larynx','fontSize',11,'verticalAlignment','bottom','horizontalAlignment','right');
		text(N/2,gLen+.1,'Phase (degrees)','fontSize',12,'verticalAlignment','bottom','horizontalAlignment','center');

		axes('position',[.1 .08 .82 .6]);
		line(glx,gly,'color',[.8 .8 .8]); 
		k = round([0:8]*N/8);
		k(1) = 1;
		line(glx(:,k),gly(:,k),'color','b');
		box on; hold on;
		h = plot(xy0(:,1),xy0(:,2),'r','linewidth',2);
		k = [0 ; cumsum(sqrt(sum(diff(xy).^2,2)))];		
		xy = interp1(k, xy, linspace(0, k(end), state.NPOINTS)', 'pchip');
		h(2) = plot(xy(:,1),xy(:,2),'color',[0 .7 0],'linewidth',2);
		set(gca,'xlim',(gLen+.1)*[-1 1],'ylim',(gLen+.1)*[-1 1]);
		axis equal; ylabel('cm');
		legend(h,'contour','fit','location','northeast');
		for k = 1 : 5, text(2,-(k/2+1),s{k},'fontname','Courier','fontsize',10,'color','k'); end
	
	
%-----------------------------------------------------------------------------
% CONFIG:  configure

	case 'CONFIG',
		state = get(gcbf,'userData');
		par = DoConfig(state.PARAMS, state.DEFPAR);
		if isempty(par), return; end;
		state.PARAMS = par;
		state.MPP = state.PARAMS.MPP;
		state.ORIGIN = state.PARAMS.ORIGIN;
		set(gcbf,'userData',state);
		if (~isempty(state.MPP) && ~isempty(state.ORIGIN)), es = 'on'; else, es = 'off'; end;
		set(state.FCH,'enable',es);
		
	
%-----------------------------------------------------------------------------
% CYCLE:  cycle through movie frames

	case 'CYCLE',
		if get(gca,'userdata'),		% if already cycling stop
			set(gca,'userdata',0); 
			return;
		end;
		state = get(gcbf,'userdata');
		set(state.TMH(3),'checked','off');		% ensure tracking disabled
		v = evalin('base',state.VNAME);
		frames = cell2mat({v.FRAME});			% frames with data
		f = state.CURFRAME;
		switch varargin{2},
			case 'FWD', dir = 1;
			case 'BCK', dir = -1;
			otherwise, set(gca,'userdata',0); return;	% stop cycling
		end;
		
% update output variable state
		if dir && ~isempty(frames),
			k = (state.CURFRAME == frames);
			v(k).XY = state.XY;
			v(k).ANCHORS = state.ANCHORS;
			v(k).NOTE = get(state.LH,'string');
			v(k).TRKRES = state.TRKRES;
			assignin('base',state.VNAME,v);
			set(state.LH,'string','');
			delete(findobj(gca,'tag','CONTOUR'));
			state.ALH = []; state.CLH = [];
			set(gcbf,'userdata',state);
		end;

% movie loop
		state.IMGH.UserData = dir;
		while state.IMGH.UserData,
			if ~state.IMGH.UserData, break; end;
			if dir > 0,
				f = state.CURFRAME + 1;
				if f > state.NFRAMES, f = state.NFRAMES; break; end;
			else,
				f = state.CURFRAME - 1;
				if f < 1, f = 1; break; end;
			end;
			state.CURFRAME = f;
			set(state.IH,'cdata',GetImage(state));
			if state.ISKF(f), c = 'g'; else, c = 'y'; end;
			set(state.TH,'color',c,'string',sprintf('%s',FmtFrame(f,state.FRATE)));
			n = sprintf('%s  [%d of %d]',state.FNAME,f,state.NFRAMES);
			if ~state.LOCK, n = [n , ' (editable)']; end;
			set(gcbf,'name',n, 'userdata',state);
			set(state.FRAMEH, 'string', num2str(f));
			set(state.SCROLLERH, 'value', f);

% update audio cursor
			if ~isempty(state.AUDIO),
				ff = floor(((f-1)/state.FRATE)*state.AUDIO.SRATE) + 1;
				set(state.AUDIO.ACH,'xdata',[ff;ff]);
			end;
			drawnow;
		end;
		state.IMGH.UserData = 0;
		axes(state.IMGH);
% clean up

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
		k = find(f == frames);
		if isempty(k),					% virgin frame
			if strcmpi('off',get(state.NH,'checked')),
				state.ANCHORS = [];		% don't inherit
				state.XY = [];
			end;
			vv = struct('XY',state.XY,'ANCHORS',state.ANCHORS,'FRAME',f,'NOTE',ts,'TRKRES',[]);
			v(end+1) = vv;
		else,							% update from existing anchors
			if isempty(v(k).ANCHORS) && strcmpi('on',get(state.NH,'checked')),
				v(k).XY = state.XY;		% inherit into empty frame
				v(k).ANCHORS = state.ANCHORS;
				state.PREVAP = [];
			else,
				state.PREVAP = state.ANCHORS;
			end;			
			state.ANCHORS = v(k).ANCHORS;
			state.XY = v(k).XY;
			v(k).NOTE = ts;
			state.TRKRES = v(k).TRKRES;
		end;
		assignin('base',state.VNAME,v);

		for k = 1 : size(state.ANCHORS,1),
			state.ALH(k) = MakePoint(state.ANCHORS(k,1),state.ANCHORS(k,2),state.CONFIG);
		end;
		if ~isempty(state.XY),
			state.CLH = line(state.XY(:,1),state.XY(:,2),'color',state.CONFIG.LINECOLOR,'linewidth',state.CONFIG.LINEWIDTH,'tag','CONTOUR','hitTest','off');
			uistack(state.CLH,'bottom'); uistack(state.CLH,'up');
		end;
% update contour line using active Tracker handler
		if ~isempty(state.TRACKER),
			ns = state.TRACKER('PLOT', state);
			if ~isempty(ns), state = ns; end;
		end;

		set(gcbf,'userData',state);
		

%-----------------------------------------------------------------------------
% DELETE:  delete existing anchors

	case 'DELETE',
		state = get(gcbf,'userData');
		state.PREVAP = state.ANCHORS;
		if strcmp(questdlg('Clear all anchors...', 'Verify...', 'Yes', 'No', 'Yes'), 'Yes'),
			delete(findobj(gca,'tag','CONTOUR'));
			state.CLH = []; state.ALH = []; state.ANCHORS = []; state.XY = [];
			set(gcbf,'userData',UpdateContour(state));
		end;
		

%-----------------------------------------------------------------------------
% DIFF:  toggle frame differencing

	case 'DIFF',
		if strcmp('on',get(gcbo,'checked')),
			set(gcbo,'checked','off');
		else,
			set(gcbo,'checked','on');
		end;
		state = get(gcbf,'userData');
		set(state.IH,'cdata',GetImage(state));		
	
	
%-----------------------------------------------------------------------------
% DOWN:  mouseDown handler
%
% click in image creates new anchor point
% click on an anchor repositions it
% control-click on a point deletes it
% click and drag with no previous anchor points creates XY to which anchors are fitted

	case 'DOWN',
		state = get(gcbf,'userData');
		state.PREVAP = state.ANCHORS;
		gotPoint = (length(varargin) > 1);	% nonzero for click on existing point
		mod = get(gcbf,'selectionType');
		cp = get(gca, 'currentPoint');
		cp = cp(1,1:2);

		switch mod,		

			case 'normal',			% unmodified

% move existing point
				if gotPoint,		
			
					set(gcbf, 'windowButtonMotionFcn',{@GetContours,'MOVE',gcbo}, ...
								'windowButtonUpFcn',{@GetContours,'UP',gcbo}, ...
								'pointer','crosshair', ...
								'userData',UpdateContour(state));

% draw mode:  no existing anchors
				elseif length(state.ANCHORS) == 0,
				
					lh = line(cp(1),cp(2),'color','c','linewidth',2,'tag','TEMPLINE');
					set(gcbf, 'windowButtonMotionFcn',{@GetContours,'DRAW','MOVE',lh}, ...
								'windowButtonUpFcn',{@GetContours,'DRAW','UP',lh}, ...
								'pointer','crosshair', ...
								'userData',state);

% add new point								
				else,
					if isempty(state.TRACKER), 	
						trackerAddPt = 0;	% default processing
					else,					% defer to tracker
						trackerAddPt = state.TRACKER('ADDPT',state,cp);
					end;
					if trackerAddPt < 0, return; end;	% -1 flags ignore new point
					lh = MakePoint(cp(1),cp(2),state.CONFIG);

% if new point is within existing points (less than half distance between nearest two points)
% then add it between those points, else append it to nearest end
					n = length(state.ALH);
					
% if tracker returns 1 then anchor is appended to end
					if trackerAddPt>0 && n>1, n = 1; end;	% force append to end
					switch n,
						case 0,			% first point
							state.ALH = lh;
							state.ANCHORS = cp;
						case 1,			% second point
							state.ALH = [state.ALH , lh];
							state.ANCHORS = [state.ANCHORS ; cp];
						otherwise,		
							d = sqrt(sum((ones(size(state.ANCHORS,1),1)*cp - state.ANCHORS).^2,2));
							[~,k] = min(d);			% index of closest point
							if k == 1,				% prefix to first point
								state.ALH = [lh , state.ALH];
								state.ANCHORS = [cp ; state.ANCHORS];
							elseif k == length(d),	% append to last point
								state.ALH = [state.ALH, lh];
								state.ANCHORS = [state.ANCHORS ; cp];
							else,					% insert betweeen existing points
								if d(k-1) < d(k+1), k2 = k; k = k-1; else, k2 = k+1; end;
								state.ALH = [state.ALH(1:k) , lh , state.ALH(k2:end)];
								state.ANCHORS = [state.ANCHORS(1:k,:) ; cp ; state.ANCHORS(k2:end,:)];
							end;
					end;
					set(gcbf, 'userData',UpdateContour(state));
				end;

% ignore double-click
			case 'open',
				;

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
% DRAW:  draw contour
%
%	varargin{2}:  MOVE or UP
%	varargin{3};  line handle
%
% contour drawn if no existing anchor points; on mouseUp anchors are distributed along drawn contour

	case 'DRAW',
		cp = get(gca, 'currentPoint');
		cp = cp(1,1:2);
		lh = varargin{3};
		xy = [lh.XData ; lh.YData]';
		if xy(end,1)~=cp(1) || xy(end,2)~=cp(2), 
			xy(end+1,:) = cp; 
			set(lh,'xdata',xy(:,1),'ydata',xy(:,2));
		end;

% MOVE		
		if strcmp(varargin{2},'MOVE'),
			drawnow;
			return;
		end;

% UP		
		set(gcbf, 'windowButtonMotionFcn','', 'windowButtonUpFcn','', 'pointer','arrow');
		delete(findobj(gca,'tag','TEMPLINE'));
		state = get(gcbf,'userData');
		if size(xy,1) < 2,		% add anchor point
			state.ANCHORS = cp;
			set(gcbf,'userData',UpdateContour(state));
			return;
		end;
		
% resample drawn contour
		k = [0 ; cumsum(sqrt(sum(diff(xy).^2,2)))];
		state.XY = interp1(k, xy, linspace(0, k(end), state.NPOINTS)', 'pchip');	

% find signed curvature using central differencing
		dx = gradient(state.XY(:,1)); dy = gradient(state.XY(:,2));
		ddx = gradient(dx); ddy = gradient(dy);
		k = (dx .* ddy - dy .* ddx) ./ (dx.^2 + dy.^2).^1.5;

% trim curvature to values whose associated radius is less than TRIM * path integral from first to last point
		fk = k; 
		trim = .1;
		if trim > 0,
			q = sum(sqrt(sum(diff(xy).^2,2))) * trim;
			fk(abs(1./k) > q) = 0;
		end;

% count inflections (nonzero sign changes)
		sfk = sign(fk);
		xfk = sfk(sfk ~= 0);
		if isempty(xfk),
			xfk = sign(k); xfk = xfk(xfk~=0);
			if isempty(xfk),
				nInfl = 0;			% collinear points
			else,
				nInfl = 1;			% curvature below threshold
			end;
		else,
			nInfl = sum(diff(xfk)~=0) + 1;
		end;
		if nInfl > 9, nInfl = 9; end;

		nAnchors = nInfl + 2;
		k = [0 ; cumsum(sqrt(sum(diff(state.XY).^2,2)))];
		state.ANCHORS = interp1(k,state.XY,linspace(0,k(end),nAnchors),'linear');
		set(gcbf,'userData',UpdateContour(state));

	
%-----------------------------------------------------------------------------
% EMIT:  export contours to workspace
%
% coordinates include mm values relative to origin if mm/pixel factor and origin available

	case 'EMIT',
		fh = findobj('tag','GETCONTOURS');
		state = get(fh,'userData');
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
		if (~isempty(state.MPP) && ~isempty(state.ORIGIN)),		% emit mm coordinates
			mmc(:,1,:) = contours(:,1,:) - state.ORIGIN(1);
			mmc(:,2,:) = state.ORIGIN(2) - contours(:,2,:);		% flip	
			mmc = mmc * state.MPP;
			varargout{1} = [contours , mmc];
		else,													% emit pixel coordinates
			varargout{1} = contours;
		end;
		if nargout > 1,
			varargout{2} = frames;
			if nargout > 2,		% include images
				q = get(state.IH,'cdata');
				q = zeros(size(q,1),size(q,2),length(frames),'uint8');
				for f = 1 : length(frames),
					q(:,:,f) = GetMovieFrame(state.MH,frames(f),state.CROP,state.RESIZE,state.IMGMODP,state.IMGMODA{:});
				end;
				varargout{3} = q;
			end;
		end;
		return;
		
		
%-----------------------------------------------------------------------------
% EXPORT:  save contours to tab delimited output file
%
% coordinates include mm values relative to origin if mm/pixel factor and origin available

	case 'EXPORT',
		fh = gcbf;
		if isempty(fh), fh = findobj('tag','GETCONTOURS'); end;
		state = get(fh,'userData');
		[fName,pName] = uiputfile('*.tsv','Export Contours to tab-delimited text file',state.PARAMS.FNAME);
		if fName == 0, return; end;		% cancel
		
		v = evalin('base',state.VNAME);			
		frames = cell2mat({v.FRAME});			% frames with data
		k = (state.CURFRAME == frames);
		if isempty(k), fprintf('no data for export\n'); return; end;
		gotOrigin = (~isempty(state.MPP) && ~isempty(state.ORIGIN));
		
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
		
		if gotOrigin,
			mmc(:,1,:) = contours(:,1,:) - state.ORIGIN(1);
			mmc(:,2,:) = state.ORIGIN(2) - contours(:,2,:);		% flip	
			mmc = mmc * state.MPP;
		end;
		
		fid = fopen(fName,'wt');
		if fid < 0, fprintf('Error attempting to open %s\n',fName); return; end;
		fprintf(fid,'FRAME\tTIME\tNOTE\tPOINT\tX\tY');
		if gotOrigin, fprintf(fid,'\tmmX\tmmY'); end;
		fprintf(fid,'\n');
		for fi = 1 : length(frames),
			for ci = 1 : size(contours,1),
				fprintf(fid,'%d\t%f\t%s\t%d\t%.2f\t%.2f', frames(fi), (frames(fi)-1)/sr, notes{fi}, ci, contours(ci,1,fi), contours(ci,2,fi));
				if gotOrigin, fprintf(fid,'\t%.2f\t%.2f',mmc(ci,1,fi), mmc(ci,2,fi)); end;
				fprintf(fid,'\n');
			end;
		end;
		fclose(fid);
		fprintf('wrote %s\n',fName);
		
		
%-----------------------------------------------------------------------------
% FILTER:  show image forces

	case 'FILTER',
		state = get(gcbf,'userData');
		img = ComputeImageForces(im2double(get(state.IH,'cdata')),state.PARAMS.SIGMA);
		img = im2uint8((img-min(img(:)))./range(img(:)));
		set(state.IH,'cdata',img);
		
		
%-----------------------------------------------------------------------------
% FLIP:  invert image

	case 'FLIP',
		if strcmp(varargin{2},'HORIZONTAL'),
			if strcmp(get(gca,'xdir'),'normal'), 
				set(gca,'xdir','reverse'); 
				set(get(gcbo,'userdata'),'HorizontalAlignment','right');
			else, 
				set(gca,'xdir','normal'); 
				set(get(gcbo,'userdata'),'HorizontalAlignment','left');
			end;
		else,
			if strcmp(get(gca,'ydir'),'normal'), 
				set(gca,'ydir','reverse'); 
			else, 
				set(gca,'ydir','normal'); 
			end;
		end;
		
		
%-----------------------------------------------------------------------------
% FRAME:  set image frame

	case 'FRAME',
		fh = gcbf;
		if isempty(fh), fh = findobj('tag','GETCONTOURS'); end;
		state = get(fh,'userData');
		f = state.CURFRAME;						% current frame before change
		external = 0;
		switch varargin{2},
			case 'EXPLICIT',
				external = 1;
				f = varargin{3};	% called externally
			case 'PREV', if f > 1, f = f - 1; end;
			case 'NEXT', if f < state.NFRAMES, f = f + 1; end;
			case 'SPECIFY',
				f = str2num(get(state.FRAMEH,'string'));
				if isempty(f),
					f = state.CURFRAME;
				elseif f < 1,
					f = 1;
				elseif f > state.NFRAMES,
					f = state.NFRAMES;
				end;
			case 'SCROLL',
				set(state.TMH(3),'checked','off');		% ensure tracking disabled
				f = round(get(state.SCROLLERH,'value'));
				if f < 1, f = 1; elseif f > state.NFRAMES, f = state.NFRAMES; end;
			case {'KPREV','KNEXT'},
				[~,k] = min(abs(state.KEYFRAMES(:,1)-f));
				if strcmp(varargin{2},'KPREV'),
					if k > 1, k = k - 1; end;
				else,
					if k < size(state.KEYFRAMES,1), k = k + 1; end;
				end;
				f = state.KEYFRAMES(k,1);
			case {'DPREV','DNEXT'},
				v = evalin('base',state.VNAME);
				fr = sort(cell2mat({v.FRAME}));			% frames with data
				if strcmp(varargin{2},'DPREV'),
					fr = fliplr(fr(fr<f));
				else
					fr = fr(fr>f);
				end
				if isempty(fr), return; end
				f = fr(1);
		end;
		
		if state.CURFRAME == f, return; end;	% nothing to do
		state = NewFrame(state, f, 0, external);
 		set(fh,'userData',state);
		
	
%-----------------------------------------------------------------------------
% GETMOVIEFRAME:  hook to internal GetMovieFrame; returns movie frame(s) based 
%					on current state parameters
%
%	varargin{2}:  frame(s) to return (uint8)

	case 'GETMOVIEFRAME',
		fh = findobj('tag','GETCONTOURS');
		state = get(fh,'userData');
		varargout{1} = GetMovieFrame(state.MH, varargin{2}, state.CROP, state.RESIZE, state.IMGMODP, state.IMGMODA{:});

		
%-----------------------------------------------------------------------------
% INFO:  report current frame information

	case 'INFO',
		state = get(gcbf,'userData');
		fprintf('Frame %s  %s\n', FmtFrame(state.CURFRAME, state.FRATE), get(state.LH,'string'));
		if size(state.XY,1)<state.NPOINTS || size(state.ANCHORS,1)<3, return; end;
		xy = state.XY;
		[~,k] = min(xy(:,2));
		highPt = xy(k,:);									% contour vertical max
		xy2 = [xy(end,:) ; xy(1:end-1,:)];
		areas = xy2(:,1).*xy(:,2) - xy(:,1).*xy2(:,2);
		sa = sum(areas);
		cenPt = sum((xy+xy2).*repmat(areas,1,2))./(3*sa);	% centroid
		cumDist = [0;cumsum(sqrt(sum(diff(xy).^2,2)))];
		len = cumDist(end);									% contour length
		if ~isempty(state.MPP) && ~isempty(state.ORIGIN),
			mmPts = [highPt ; cenPt];
			mmPts(:,1) = mmPts(:,1) - state.ORIGIN(1);
			mmPts(:,2) = state.ORIGIN(2) - mmPts(:,2);
			mmPts = mmPts * state.MPP;
			mmLen = len * state.MPP;
			fprintf('  high point:  [%4.0f,%4.0f] (pixels)   [%5.1f,%5.1f] (mm)\n', highPt, mmPts(1,:));
			fprintf('    centroid:  [%4.0f,%4.0f] (pixels)   [%5.1f,%5.1f] (mm)\n', cenPt, mmPts(2,:));
			fprintf('      length:  %.0f (pixels)   %.1f (mm)\n', len, mmLen);
		else,
			fprintf('  high point:  [%4.0f,%4.0f] (pixels)\n', highPt);
			fprintf('    centroid:  [%4.0f,%4.0f] (pixels)\n', cenPt);
			fprintf('      length:  %.0f (pixels)\n', len);
		end;
		try,
			[~,ninfl,mci] = ComputeCurvature(xy);
			fprintf('       NINFL:  %d    MCI: %.2f\n', ninfl, mci);
		catch,
			;
		end;
		if ~isempty(state.TRKRES) && isnumeric(state.TRKRES),
			fprintf('      TRKRES:  %.1f\n', state.TRKRES);
		end;

		
%-----------------------------------------------------------------------------
% MAP:  set colormap

	case 'MAP',
		if strcmp(get(gcbo,'text'),'Reset Original Image'),
			mh = get(gcbo,'userdata');
		else,
			mh = gcbo;
		end;
		set(get(get(mh,'parent'),'children'),'checked','off');
		set(mh,'checked','on');
		mapName = varargin{2};
		if strcmp(mapName,'Inv Gray'),
			map = 1-gray;
		else,
			map = eval(lower(get(mh,'label')));
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
			state.TRKRES = [];					% invalidate tracker result
			state.ANCHORS(k,:) = [x,y];
			set(gcbf,'userData',UpdateContour(state,1));
		end;


%-----------------------------------------------------------------------------
% MPP:  find mm/pixel ratio

	case 'MPP',
		fh = gcbf;
		if isempty(fh), fh = findobj('tag','GETCONTOURS'); end;
		state = get(fh,'userData');
		state.USEAVG = 0;
		img = get(state.IH,'cdata');
		cm = GetMPP;		% get mpp cm units (1 or .5)
		if isempty(cm), return; end;

% get rectangular selection
		rect = getrect(state.IH.Parent);
		r = round([rect(2),rect(2)+rect(4)-1]);
		c = round([rect(1),rect(1)+rect(3)-1]);
		[h,w] = size(img);
		if r(1) < 1, r(1) = 1; end;
		if r(2) > h, r(2) = h; end;
		if c(1) < 1, c(1) = 1; end;
		if c(2) > w, c(2) = w; end;
		img = img(r(1):r(2),c(1):c(2));
		img = img > 100;		% threshold

% find difference between nonzero stripes along sweep
% across height of rect at horizontal point of greatest intensity
		[v,k] = max(sum(img));
		s = img(:,k)*2;
		idx = find(diff([0;s]) > 1);
		len = diff([idx;length(s)+1]);
		dy = median(len(1:end-1));

% mpp
		mpp = 10*cm/dy;
		state.MPP = mpp;
		state.PARAMS.MPP = mpp;
		set(fh, 'userdata', state);
		if (~isempty(state.MPP) && ~isempty(state.ORIGIN)), es = 'on'; else, es = 'off'; end;
		set(state.FCH,'enable',es);
		if nargout < 1,
			fprintf('mm / pixel ratio: %.3f\n', mpp);
		else,
			varargout{1} = mpp;
		end;
		
		
%-----------------------------------------------------------------------------
% ORIGIN:  find origin of US probe arc

	case 'ORIGIN',
		fh = gcbf;
		if isempty(fh), fh = findobj('tag','GETCONTOURS'); end;
		state = get(fh,'userData');
		state.USEAVG = 0;
		img = get(state.IH,'cdata');
		if strcmp(questdlg('Use the interactive crosshairs to select a point below the probe surface arc', ... 
			'Find Origin...', 'Proceed', 'Cancel', 'Proceed'), 'Cancel'), return; end;

% get seed point
		isBlk = (sum(img(:)>128) < sum(img(:)<128));	% nonzero for predominantly black
		if isBlk, set(state.IH,'CDATA',255-img); end;	% invert for visibility with ginput
		xy = round(ginput(1));
		if isBlk, set(state.IH,'CDATA',img); end;
		if isempty(xy), return; end;
		
% set threshold		
		[c,x] = imhist(img);		
		[~,k] = max(c(5:end-4));	% skip lowest/highest values
		thresh = x(k+15);			% set thresh 10 levels higher than max

% find probe surface
		img = img > thresh;
		y = find(img(1:xy(2),xy(1)),1,'last');		% point on outline
		try,
			img = bwselect(img, xy(1), y);
		catch,
			fprintf('unable to detect probe surface with specified seed point\n');
			return;
		end;

% find arc
		r = find(any(img,2));						% find rows with pixels
		clip = r(end) - 100;						% keep bottom 100 rows
		img(1:clip,:) = 0;
		split = round(mean(find(any(img,1))));		% left/right split
		r = find(any(img(:,1:split),2),1,'last');	% row with lowest pixel left side
		c = find(img(r,1:split));					% leftmost col(s) this row
		c = round(mean(c));							% take center of multiple matches
		minL = [c,r];								% left edge of arc
		r = find(any(img(:,split:end),2),1,'last');	% row with lowest pixel right side
		c = find(img(r,split:end)) + split - 1;		% rightmost col(s) this row
		c = round(mean(c));							% take center of multiple matches
		minR = [c,r];								% right edge of arc
		B = fliplr(bwtraceboundary(img,fliplr(minR),'W'));	% trace outline [nPts x x,y]
		k = find(all(B==ones(size(B,1),1)*minL,2));
		arc = B(1:k,:);								% surface right -> left, [nPts x x,y]
		minLR = [minL ; minR];

% fit circle to arc 
		M = [arc(:,1).^2 + arc(:,2).^2 , arc , ones(size(arc,1),1)];
		[~,~,v] = svd(M);
		origin = [-v(2:3,4) / (2*v(1,4))]';
		
% plot (tagged as CONTOUR)
		line(arc(:,1),arc(:,2),'color','g','tag','CONTOUR');
		line([arc([1 end],1)';origin([1 1])],[arc([1 end],2)';origin([2 2])],'color','g','linestyle',':','tag','CONTOUR');
		line(origin(1),origin(2),'color','g','marker','+','markersize',20,'tag','CONTOUR');
	
		state.ORIGIN = origin;
		state.PARAMS.ORIGIN = origin;
		set(fh, 'userdata', state);
		if (~isempty(state.MPP) && ~isempty(state.ORIGIN)), es = 'on'; else, es = 'off'; end;
		set(state.FCH,'enable',es);
		if nargout < 1,
			fprintf('Origin (center of circle fit to probe surface, pixels):  [%.1f %.1f]\n', origin);
		else,
			varargout{1} = origin;
		end;
		
%-----------------------------------------------------------------------------
% PLAY:  play audio & video +/- PLAYINT centered on current frame

	case 'PLAY',
		state = get(gcbf,'userData');
		ht = (state.CURFRAME-1)/state.FRATE + state.PARAMS.PLAYINT/1000;

% audio
		try,
			info = audioinfo(state.FNAMEEXT);
			hts = floor(ht*info.SampleRate) + 1;	% audio interval (samples)
			if hts(1) < 1, hts(1) = 1; end;
			if hts(2) > info.TotalSamples, hts(2) = info.TotalSamples; end;
			[s,sr] = audioread(state.FNAMEEXT,hts);
			PLAYERH = audioplayer(s,sr);
			play(PLAYERH);
		catch,
			;
		end;

% video
		hts = floor(ht*state.FRATE) + 1;		% movie interval (frames)
		if hts(1) < 1, hts(1) = 1; end;
		if hts(2) > state.NFRAMES, hts(2) = state.NFRAMES; end;
		f0 = state.CURFRAME;
		for f = hts(1) : hts(2),
			state.CURFRAME = f;
			set(state.IH,'cdata',GetImage(state));
			if state.ISKF(f), c = 'g'; else, c = 'y'; end;
			set(state.TH,'color',c,'string',sprintf('%s',FmtFrame(f,state.FRATE)));
			drawnow;
		end;
		state.CURFRAME = f0;
		set(state.IH,'cdata',GetImage(state));
		if state.ISKF(f0), c = 'g'; else, c = 'y'; end;
		set(state.TH,'color',c,'string',sprintf('%s',FmtFrame(f0,state.FRATE)));
		drawnow;

		
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
% RESAMPLE:  resample anchor points along current contour

	case 'RESAMPLE',
		state = get(gcbf,'userData');
		nAnchors = GetNAP(state.NAPRESAMP);
		if isempty(nAnchors), return; end;
		state.NAPRESAMP = nAnchors;
		state.PREVAP = state.ANCHORS;
		k = [0 ; cumsum(sqrt(sum(diff(state.XY).^2,2)))];
		state.ANCHORS = interp1(k,state.XY,linspace(0,k(end),nAnchors),'linear');
		delete(state.ALH);
		state.ALH = [];
		set(gcbf,'userData',UpdateContour(state));
		
		
%-----------------------------------------------------------------------------
% RESET:  reset to original image

	case 'RESET',
		state = get(gcbf,'userData');
		state.USEAVG = 0;
		set(state.DH,'checked','off');
		set(state.IH,'cdata',GetImage(state));
		set(gcbf,'userData',state);
		set(gca,'xdir','normal','ydir','reverse','clim',state.CLIM);
		set(state.TH,'HorizontalAlignment','left');
		GetContours('MAP','Gray');
		
		
%-----------------------------------------------------------------------------
% SNAKE:  apply snake to current frame

	case 'SNAKE',
		state = get(gcbf,'userData');
		nAnchors = size(state.ANCHORS,1);
		if nAnchors < 3, 
			fprintf('Snake requires at least three anchor points\n');
			return;
		end;
		state.PREVAP = state.ANCHORS;
		set(gcbf,'pointer','watch'); drawnow;
		img = im2double(GetImage(state));
		xy = state.XY;
		k = [0 ; cumsum(sqrt(sum(diff(xy).^2,2)))];
		try,
			eg = ComputeImageForces(img,state.PARAMS.SIGMA);
			useBE = 1;		% use band energy
			xy = interp1(k, xy, linspace(0, k(end), state.PARAMS.NSPTS)', 'pchip');
			xy = make_snake(img', eg', xy, state.PARAMS.DELTA*ones(state.PARAMS.NSPTS,1), state.PARAMS.BPEN, state.PARAMS.ALPHA, state.PARAMS.LAMBDA, useBE);
		catch,
			fprintf('failed to apply snake\n');
			return;
		end;
		set(gcbf,'pointer','arrow');
		k = [0 ; cumsum(sqrt(sum(diff(xy).^2,2)))];		
		state.XY = interp1(k, xy, linspace(0, k(end), state.NPOINTS)', 'pchip');
		set(state.CLH,'xdata',state.XY(:,1),'ydata',state.XY(:,2));
		state.ANCHORS = interp1(k, xy, linspace(0, k(end), nAnchors)', 'pchip');
		for k = 1 : nAnchors,				% redistribute anchors
			set(state.ALH(k),'xdata',state.ANCHORS(k,1),'ydata',state.ANCHORS(k,2));
		end;
		drawnow;
		set(gcbf,'userData',state);
		
	
%-----------------------------------------------------------------------------
% TRACK:  process automatic tracking

	case 'TRACK',
		state = get(gcbf,'userData');
		
		switch varargin{2},
			case 'CLEAR',		% clear tracker
				state.TRACKER = [];
				state.TPAR = [];
				state.TRKRES = [];
				set(state.TMH(1),'label','<None>','enable','off');
				set(state.TMH(2:end),'enable','off');
				delete(findobj(gca,'tag','CONTOUR'));
				state.CLH = []; state.ALH = []; state.ANCHORS = [];
				
			case 'SELECT',		% select tracker
				[fn,p] = uigetfile([fileparts(which('GetContours')),filesep,'gct*.m'], 'Select tracking procedure...');
				if fn == 0; return; end;		% cancelled
				[~,fn] = fileparts(fn);
				d = pwd; Tracker = str2func(fn); cd(d);			% function handle
				state.TRACKER = Tracker;
				set(state.TMH(1),'label',[upper(fn), ':  configure'],'enable','on');		% set name into menu
				set(state.TMH(2:end),'enable','on');
				state.TPAR = Tracker('CONFIG', state, 1);		% configure
				state.TRKRES = [];
				
			case 'EXPORT', 		% export tracker data to file
				state.TRACKER('EXPORT', state);
				
			case 'CONFIG',		% configure tracker
				cfg = state.TRACKER('CONFIG', state, 0);
				if ~isempty(cfg), state.TPAR = cfg; end;
				
			case 'APPLY',		% apply tracker to current frame
				if isempty(state.TPAR),
					cfg = state.TRACKER('CONFIG', state, 1);
					if ~isempty(cfg), state.TPAR = cfg; end;
					set(gcbf,'userData',state);
					return;
				end;
				state.PREVAP = state.ANCHORS;
				set(gcbf,'pointer','watch'); drawnow;
				firstTime = 1;
				ns = state.TRACKER('TRACK', state, firstTime);
				if isempty(ns), 		% tracking failure
					set(gcbf,'pointer','arrow','userData',state); 
					drawnow;
					return; 
				end;				
				state = ns;
				firstTime = 0;

% update output variable state
				v = evalin('base',state.VNAME);
				k = find(state.CURFRAME == cell2mat({v.FRAME}));
				v(k).XY = state.XY;
				v(k).ANCHORS = state.ANCHORS;
				v(k).NOTE = get(state.LH,'string');
				v(k).TRKRES = state.TRKRES;
				assignin('base',state.VNAME,v);
				
% track thru sequence
				kf = state.KEYFRAMES;
				if isempty(kf), kf = [1 state.NFRAMES]; end;
				while strcmpi('on',get(state.TMH(3),'checked')),
					f = state.CURFRAME + 1;
					k = find(f>=kf(:,1) & f<=kf(:,2));
					if isempty(k), break; end;			% end of interval

% set new frame
					state = NewFrame(state, f, 1);
					set(gcbf,'userData',state);

% track
					ns = state.TRACKER('TRACK', state, firstTime);
					drawnow;
					if isempty(ns), break; end;			% tracking failure
					state = ns;

% update output new frame state
					v = evalin('base',state.VNAME);
					k = find(state.CURFRAME == cell2mat({v.FRAME}));
					v(k).XY = state.XY;
					v(k).ANCHORS = state.ANCHORS;
					v(k).NOTE = get(state.LH,'string');
					v(k).TRKRES = state.TRKRES;
					assignin('base',state.VNAME,v);						
				end;
				set(gcbf,'pointer','arrow'); drawnow;
    
			case 'SAVE', 		% save tracker data to base ws (note: save option must be added to menu explicitly by tracker; cf. gct_lines)
				state.TRACKER('SAVE', state);
				
			case 'DIAGNOSTICS', % show tracker diagnostics (note: save option must be added to menu explicitly by tracker; cf. gct_SLURP)
				state.TRACKER('DIAGS', state);
				
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
		set(lh, 'xdata',cp(1), 'ydata',cp(2));
		state = get(gcbf,'userData');
		state.TRKRES = [];					% invalidate tracker result
		k = find(lh == state.ALH);
		state.ANCHORS(k,:) = cp;
		set(gcbf,'userData',UpdateContour(state));
		
		
%-----------------------------------------------------------------------------
% INITialize

	otherwise,
		fh = findobj('tag','GETCONTOURS');
		if isempty(fh),
			Initialize(varargin{:});
		else,
			varargin{end+1} = 'LOCK';		% only first object can edit
			varargin{end+1} = true;
			Initialize(varargin{:});
		end;
		if isunix, [s,r] = unix('osascript -e ''tell application "MATLAB" to activate'''); end;
end;


%=============================================================================
% COEFTOSHAPE  - map Fourier coefficients to tongue shape

function [xy,d,mag,phi] = CoefToShape(C, S, DC, glx, gly)

C = C(:); S = S(:);
N = size(glx,2);					% # grid lines
gLen = abs(glx(1,end));				% grid line length
n = [0 : N-1];
k = [1 : length(C)]';

% distance function (lips to larynx)
dd = C*ones(1,N).*cos(2*pi*k*n./N) + S*ones(1,N).*sin(2*pi*k*n./N);		
d = DC + sum(dd,1);		% composite
dd = DC + dd;			% contribution by coef

% semi-polar gridlines
idx = find(diff(gly),1,'last');		% index of first pharyngeal gridline
th = cart2pol(glx(1,1:idx),gly(1,1:idx));
[x,y] = pol2cart(th,d(1:idx));
[glx(2,1:idx),gly(2,1:idx)] = pol2cart(th,.5);

% pharyngeal gridlines
y(idx:N) = gly(1,idx:N);
x(idx:N) = d(idx:N);
xy = [x(:),y(:)];
flipped = diff(glx(:,end)) > 0;
if flipped,			% faces right
	xy(idx:end,1) = -xy(idx:end,1); 
	glx(2,idx:end) = -.5;
end;	

% magnitude (constriction degree) and phase (constriction location)
mag = sqrt(C.^2 + S.^2);
phi = (180/pi)*atan(S./C)./k;

% return distance function as [N x nCoef]
d = dd';


%=============================================================================
% COMPUTEIMAGEFORCES  - computes Gaussian-filtered image derivatives

function img = ComputeImageForces(img, sigma)

[x,y] = ndgrid(floor(-3*sigma):ceil(3*sigma),floor(-3*sigma):ceil(3*sigma));
g = exp(-(x.^2 + y.^2)/(2*sigma^2));
d = 2*pi*sigma^4;
ix = imfilter(img, -(x./d).*g);
iy = imfilter(img, -(y./d).*g);
img = 1 - sqrt(ix.*ix + iy.*iy);


%=============================================================================
% CURSORADJUST  - set current frame based on click in audio panel

function CursorAdjust(ah, evt)

imgH = ah.UserData;
if imgH.UserData, imgH.UserData = 0; return; end;	% if cycling exit loop and update there
cp = get(ah, 'currentPoint');
x = cp(1,1);
axes(imgH);			% reset image axis as foreground
fh = ah.Parent;
state = get(fh,'userData');
f = floor(((x-1)/state.AUDIO.SRATE)*state.FRATE) + 1;	% new frame
state = NewFrame(state,f,0);
set(fh,'userData',state);


%=============================================================================
% DISTTOCOEF  - compute Fourier coefficients from VT distance function

function [DC,C,S,mag,phase] = DistToCoef(d, nCoef)

d = d(:)';
N = length(d);
n = [0:N-1];
if nCoef > floor(N/2), nCoef = floor(N/2); end;
k = [1 : nCoef]';

DC = 1/N * nansum(d);
C = 2/N * nansum(cos(2*pi*k*n./N) .* (ones(length(k),1)*d) , 2);
S = 2/N * nansum(sin(2*pi*k*n./N) .* (ones(length(k),1)*d) , 2);

mag = sqrt(C.^2 + S.^2);
phase = atan(S./C);


%=============================================================================
% DOCONFIG  - configure parameters

function par = DoConfig(par, dPar)

figPos = get(0, 'ScreenSize');
width = 290; height = 420;
figPos = [figPos(1)+(figPos(3)-width)/2, figPos(2)+(figPos(4)-height)/2, width, height];

cfg = dialog('Name', 'Configure Parameters', ...
	'tag', 'GETCONTOURS', ...
	'menubar', 'none', ...
	'Position', figPos, ...
	'KeyPressFcn', 'set(gcbf,''UserData'',1);uiresume', ...
	'UserData', 0);

% export filename
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Export Filename:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0109 0.8929 0.5469 0.0464]);
fn = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', par.FNAME, ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.5687 0.8929 0.3281 0.0536]);

% save images
si = uicontrol(cfg, ...
	'Style', 'checkbox', ...
	'HorizontalAlignment', 'left', ...
	'String', 'Save Images', ...
	'Value', par.SAVEIMG, ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.3937 0.8107 0.5469 0.0536]);

% # of contour points
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','# of Contour Points:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0109 0.7286 0.5469 0.0464]);
cp = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', num2str(par.NPOINTS), ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.5687 0.7286 0.3281 0.0536]);

% play interval
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Play Interval (ms):', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0109 0.6536 0.5469 0.0464]);
pint = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%d %d',par.PLAYINT), ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.5687 0.6536 0.3281 0.0536]);

% mpp
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','mm/pixel:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0109 0.5786 0.2188 0.0464]);
mpp = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%.3f',par.MPP), ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.2406 0.5786 0.21 0.0536]);

% origin
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Origin:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.455 0.5786 0.1531 0.0464]);
ogn = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%.0f %.0f',par.ORIGIN), ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.6125 0.5786 0.2778 0.0536]);

% image forces Gaussian Sigma
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Image Forces Sigma:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0109 0.5036 0.5469 0.0464]);
gs = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%d %d',par.SIGMA), ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.5687 0.5036 0.3281 0.0536]);

% snake frame
ph = uipanel('parent',cfg,'title','Snake Parameters','position',[.05 .12 .905 .35]);

% delta
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Delta:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0766 0.3750 0.4812 0.0464]);
ds = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%.1f',par.DELTA), ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.5687 0.3750 0.3281 0.0536]);
	
% band penalty
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Band Penalty:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0766 0.3000 0.4812 0.0464]);
bps = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%.1f',par.BPEN), ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.5687 0.3000 0.3281 0.0536]);
	
% alpha
cbs1='set(get(gcbo,''userdata''),''string'',sprintf(''%.2f'',get(gcbo,''value'')))';
cbs2='if ~isempty(str2num(get(gcbo,''string'')))&&str2num(get(gcbo,''string''))>=0&&str2num(get(gcbo,''string''))<=1,set(get(gcbo,''userdata''),''value'',str2num(get(gcbo,''string'')));else,set(gcbo,''string'',sprintf(''%0.2f'',get(get(gcbo,''userdata''),''value'')));end';
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Alpha:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0656 0.2250 0.1850 0.0464]);
als = uicontrol(cfg, ...
	'Style', 'slider', ...
	'Value', par.ALPHA, ...
	'Min', 0, ...
	'Max', 1, ...
	'Callback', cbs1, ...
	'Units', 'normalized', ...
	'Position', [0.2675 0.2253 0.4325 0.0436]);
ale = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%0.2f',par.ALPHA), ...
	'FontSize', 10, ...
	'Callback', cbs2, ...
	'UserData', als, ...
	'Units', 'normalized', ...
	'Position', [0.7219 0.2270 0.175 0.0536]);
als.UserData = ale;
	
% lambda
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Lambda:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0656 0.1500 0.1850 0.0464]);
las = uicontrol(cfg, ...
	'Style', 'slider', ...
	'Value', par.LAMBDA, ...
	'Min', 0, ...
	'Max', 1, ...
	'Callback', cbs1, ...
	'Units', 'normalized', ...
	'Position', [0.2675 0.1503 0.4325 0.0436]);
lae = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%0.2f',par.LAMBDA), ...
	'FontSize', 10, ...
	'Callback', cbs2, ...
	'UserData', las, ...
	'Units', 'normalized', ...
	'Position', [0.7219 0.1500 0.175 0.0536]);
las.UserData = lae;
	

% OK, Defaults, cancel buttons
uicontrol(cfg, ...
	'String','OK', ...
	'FontSize', 10, ...
	'Callback','set(gcbf,''UserData'',1);uiresume', ...
	'Units', 'normalized', ...
	'Position', [0.1844 0.0333 0.1875 0.0595]);
uicontrol(cfg, ...
	'String','Defaults', ...
	'FontSize', 10, ...
	'Callback','set(gcbf,''UserData'',2);uiresume', ...
	'Units', 'normalized', ...
	'Position', [0.4031 0.0333 0.1875 0.0595]);
uicontrol(cfg, ...
	'String','Cancel', ...
	'FontSize', 10, ...
	'Callback','set(gcbf,''UserData'',0);uiresume', ...
	'Units', 'normalized', ...
	'Position', [0.6219 0.0333 0.1875 0.0595]);

% wait for input
while 1,
	uiwait(cfg);
	if ~ishandle(cfg),						% window closed
		par = []; 
		break;
	end;
	
	fName = strtok(get(fn,'string'));
	if isempty(fName), 
		fName = dPar.FNAME;
		continue;
	end;
	
	nPoints = str2num(get(cp,'string'));
	if isempty(nPoints),
		set(cp,'string',num2str(dPar.NPOINTS));
		continue;
	elseif nPoints<10,
		set(cp,'string','10');
		continue;
	elseif nPoints>100,
		set(cp,'string','100');
		continue;
	end;
	
	playInt = str2num(get(pint,'string'));
	if isempty(playInt),
		set(pint,'string',sprintf('%d %d',dPar.PLAYINT));
		continue;
	elseif length(playInt) == 1,
		playInt = [-1 1]*playInt;
		set(pint,'string',sprintf('%d %d',playInt));
		continue;
	elseif length(playInt) > 2,
		playInt = playInt(1:2);
		set(pint,'string',sprintf('%d %d',playInt));
		continue;
	end;
	playInt = sort(playInt);
	
	mppVal = str2num(get(mpp,'string'));
	
	origin = str2num(get(ogn,'string'));
	if ~isempty(origin) && length(origin)~=2, origin = dPar.ORIGIN; end;
	
	sigma = str2num(get(gs,'string'));
	if isempty(sigma),
		set(ds,'string',num2str(dPar.SIGMA));
		continue;
	end;

	delta = str2num(get(ds,'string'));
	if isempty(delta),
		set(ds,'string',sprintf('0.1f',dPar.DELTA));
		continue;
	end;

	bpen = str2num(get(bps,'string'));
	if isempty(bpen),
		set(bps,'string',sprintf('0.1f',dPar.BPEN));
		continue;
	end;

	switch get(cfg,'userdata'),
		case 0, 		% cancel
			par = []; 
			break;
		case 1,			% ok
			par.FNAME = fName;
			par.NPOINTS = nPoints;
			par.PLAYINT = playInt;
			par.MPP = mppVal;
			par.ORIGIN = origin;
			par.SIGMA = sigma;
			par.SAVEIMG = get(si,'value');
			par.DELTA = delta;
			par.BPEN = bpen;
			par.ALPHA = get(als,'value');
			par.LAMBDA = get(las,'value');
			break;
		case 2,			% defaults
			set(fn,'string',dPar.FNAME);
			set(cp,'string',num2str(dPar.NPOINTS));
			set(pint,'string',sprintf('%d %d',dPar.PLAYINT));
			set(mpp,'string',sprintf('%.3f',dPar.MPP));
			set(ogn,'string',sprintf('%.0f %.0f',dPar.ORIGIN));
			set(gs,'string',num2str(dPar.SIGMA));
			set(si,'value',dPar.SAVEIMG);
			set(ds,'string',sprintf('%0.1f',dPar.DELTA));
			set(bps,'string',sprintf('%0.1f',dPar.BPEN));
			set(als,'value',dPar.ALPHA); set(ale,'string',sprintf('%.2f',dPar.ALPHA));
			set(las,'value',dPar.LAMBDA); set(lae,'string',sprintf('%.2f',dPar.LAMBDA));
			continue;
	end;
end;
if ishandle(cfg), delete(cfg); end;


%=============================================================================
% FMTFRAME  - format displayed frame/time string

function s = FmtFrame(f, frate)

t = (f-1) / frate;		% time in seconds
min = floor(t/60);
sec = mod(floor(t),60);
cs = mod(floor(t*100),100);
s = sprintf('%04d  %02.0f:%02.0f:%02.0f', f, min, sec, cs);


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
		img = uint8(mean(GetMovieFrame(state.MH,[pre post],state.CROP,state.RESIZE,state.IMGMODP,state.IMGMODA{:}),3));
	catch
		fail = 1;
	end
elseif strcmpi('on',get(state.DH,'checked')),	% frame differencing enabled (curFrame+1 - curFrame)
	f = state.CURFRAME;
	pre = f;
	post = f + 1;
	if post > state.NFRAMES, pre = f-1; post = f; end;
	try,
		img = diff(GetMovieFrame(state.MH,[pre post],state.CROP,state.RESIZE,state.IMGMODP,state.IMGMODA{:}),[],3);
	catch
		fail = 1;
	end
else,
	try,
		img = GetMovieFrame(state.MH,state.CURFRAME,state.CROP,state.RESIZE,state.IMGMODP,state.IMGMODA{:});
	catch
		fail = 1;
	end
end;
if fail,
	img = zeros(get(state.MH,'Height'),get(state.MH,'Width'),'uint8');
	fprintf('unable to load frame %d from %s\n', state.CURFRAME, get(state.MH,'name'));
end


%=============================================================================
% GETMOVIEFRAME  - get frame(s) from open video object
%
% frame can be scalar or a [first,last] vector specifying frames to return
% result is averaged across R,G,B planes; i.e. [nRows x nCols x nFrames]
% and returned as uint8
%
% optional CROP is [Xmin Ymin width height]
% optional RESIZE is scaling factor
% optional IMGMOD is a function pointer with optional args applied to the image

function img = GetMovieFrame(mh, frame, crop, reSize, ImgMod, varargin)

if nargin < 3, crop = []; end;
if nargin < 4, reSize = []; end;
if nargin < 5, ImgMod = []; end;

if ischar(mh),		% load from DICOM
	if length(frame) > 1, frame = [frame(1):frame(2)]; end;
	try,
		img = dicomread(mh, 'frames', frame);
	catch,
		error('unable to read frame %d from specified DICOM file (%s)',frame(1),mh); 
	end;
else,
	try, 
		img = read(mh, frame); 
	catch, 
		error('unable to read frame %d from specified movie handle',frame(1)); 
	end;
end;

% average across RGB
img = squeeze(mean(img,3));	
nImg = size(img,3);

% crop
if ~isempty(crop), 
	if nImg == 1, 
		img = imcrop(img, crop);
	else,
		img0 = img;
		q = imcrop(img(:,:,1), crop);
		img = zeros(size(q,1),size(q,2),nImg);
		for k = 1 : size(img,3), img(:,:,k) = imcrop(img0(:,:,k), crop); end;
		clear img0;
	end;
end;

% resize
if ~isempty(reSize) && reSize~=1,

% filter when scale < 1 to avoid aliasing
	if reSize < 1,
		[ih,iw] = size(img(:,:,1));
		mih = floor(ih*reSize); miw = floor(iw*reSize);
		h = fir1(10, mih/ih)' * fir1(10, miw/iw);
		for k = 1 : nImg, img(:,:,k) = filter2(h, img(:,:,k)); end;
	end;

% build interpolation indices
	[ih,iw,~] = size(img);
	mih = floor(ih*reSize); miw = floor(iw*reSize);
	[x,y] = meshgrid(1:(iw-1)/(miw-1):iw, 1:(ih-1)/(mih-1):ih);

% interpolate
	if nImg > 1,
		img0 = img;
		img = zeros(mih,miw,nImg);
		for k = 1 : nImg, img(:,:,k) = interp2(img0(:,:,k), x, y, 'cubic'); end;
		clear img0;
	else,
		img = interp2(img, x, y, 'cubic');
	end;
end;

% apply image modifier
if ~isempty(ImgMod),
	if nImg > 1,
		q = ImgMod(img(:,:,1), varargin{:});
		[ih,iw] = size(q);
		img0 = img;
		img = zeros(ih,iw,nImg);
		img(:,:,1) = q;
		for k = 2 : nImg, img(:,:,k) = ImgMod(img0(:,:,k), varargin{:}); end;
		clear img0;
	else,
		img = ImgMod(img, varargin{:}); 
	end;
end;

% output as uint8
img = uint8(img);


%=============================================================================
% GETMPP  - get MPP cm units choice
%
% returns [] on cancel

function cm = GetMPP

width = 260;
height = 110;
fPos = get(gcf, 'position');
pos = [fPos(1)+(fPos(3)-width)/2 , fPos(2)+(fPos(4)-height)/2 , width , height];

cfg = dialog('name', 'Find mm/pixel ratio...', ...
	'menubar', 'none', ...
	'position', pos, ...
	'keyPressFcn', 'set(gcbf,''UserData'',0);uiresume', ...
	'UserData', 0);

% blurb
uicontrol(cfg, ...
	'position', [10 height-45 width-20 40], ...
	'fontsize', 12, ...
	'style', 'text', ...
	'string', 'Use the interactive rectangular selection tool to enclose at least five calibration markings');

% cm choice radio buttons
cb = 'set(gcbo,''value'',1);set(get(gcbo,''userdata''),''value'',0)';
uicontrol(cfg, ...
	'Style', 'text', ...
	'HorizontalAlignment', 'right', ...
	'String', 'Markings show', ...
	'FontSize', 10, ...
	'Position', [15 height-65 80 20]);
r1 = uicontrol(cfg, ...
	'Style', 'radiobutton', ...
	'String','1.0 cm', ...
	'Value', 0, ...
	'Position',[100 height-65 60 24], ...
	'callback', cb);
r2 = uicontrol(cfg, ...
	'Style', 'radiobutton', ...
	'String','0.5 cm', ...
	'Value', 1, ...
	'Position',[165 height-65 60 24], ...
	'userdata', r1, ...
	'callback', cb);
r1.UserData = r2;

% OK, cancel buttons
uicontrol(cfg, ...		% buttons
	'Position',[width/2-70 15 60 25], ...
	'String','Proceed', ...
	'Callback','set(gcbf,''UserData'',1);uiresume');
uicontrol(cfg, ...
	'Position',[width/2+10 15 60 25], ...
	'String','Cancel', ...
	'Callback','set(gcbf,''UserData'',0);uiresume');

% wait for input
cm = [];
while 1,
	uiwait(cfg);
	if ishandle(cfg) && get(cfg, 'UserData'),
		if get(r1,'value'),
			cm = 1;
		else,
			cm = .5;
		end;
		break;
	else,
		break;
	end;
end;
if ishandle(cfg), delete(cfg); end;


%=============================================================================
% GETNAP  - get # anchor points
%
% returns [] on cancel

function nap = GetNAP(nAnchors)

width = 250;
height = 100;
fPos = get(gcf, 'position');
pos = [fPos(1)+(fPos(3)-width)/2 , fPos(2)+(fPos(4)-height)/2 , width , height];

cfg = dialog('name', 'Specify # Anchor Points...', ...
	'menubar', 'none', ...
	'position', pos, ...
	'keyPressFcn', 'set(gcbf,''UserData'',1);uiresume', ...
	'UserData', 0);

% entry field
uicontrol(cfg, ...
	'position', [20 60 120 20], ...
	'style', 'text', ...
	'horizontalAlignment', 'right', ...
	'string', 'Resample to # Points:');
eh = uicontrol(cfg, ...
	'position', [150 60 40 25], ...
	'style', 'edit', ...
	'string', num2str(nAnchors));

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
while 1,
	uiwait(cfg);
	if ishandle(cfg) && get(cfg, 'UserData'),
		nap = round(str2num(get(eh, 'string')));
		if isempty(nap), continue; end;
		if nap < 3, nap = 3; elseif nap > 50, nap = 50; else, break; end;
		set(eh,'string',num2str(nap));
	else,
		nap = [];
		break;
	end;
end;
if ishandle(cfg), delete(cfg); end;


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
anchors = [];
audio = 1;
clean = 0;
crop = [];
inherit = 0;
keyFrames = [];
first = [];
mpp = [];
nPoints = 100;
origin = [];
reSize = [];
textgrid = [];
ImageMod = [];
Tracker = [];
tState = [];
cfg = [];
lock = false;
defCfg = struct('DOTSIZE',25,'DOTCOLOR','r','LINEWIDTH',1,'LINECOLOR','y');
if ispc, defCfg.DOTSIZE = 15; end;

[~,vName] = fileparts(fName);
for ai = 2 : 2 : length(varargin),
	switch upper(varargin{ai-1}),
		case 'ANCHORS', anchors = varargin{ai};
		case 'AUDIO', audio = varargin{ai};
		case 'CLEAN', clean = varargin{ai};
		case 'CONFIG', cfg = varargin{ai}; 
		case 'CROP', crop = varargin{ai}; 
		case 'IMGMOD', ImageMod = varargin{ai}; 
		case 'FRAME', first = varargin{ai};
		case 'KEYFRAME', keyFrames = varargin{ai};
		case 'LOCK', lock = varargin{ai};
		case 'MPP', mpp = varargin{ai};
		case 'NPOINTS', nPoints = varargin{ai};
		case 'ORIGIN', origin = varargin{ai};
		case 'RESIZE', reSize = varargin{ai};
		case 'TEXTGRID',textgrid = varargin{ai};
		case 'TRACKER', Tracker = varargin{ai};
		case 'VNAME',vName = varargin{ai};
		otherwise, error('unsupported option (%s)', varargin{ai-1});
	end;
end;

% override display defaults
if ~isempty(cfg),
	fnc = fieldnames(cfg);
	fnd = fieldnames(defCfg);
	for fi = 1 : length(fnc),
		k = find(strcmpi(fnc{fi},fnd));
		if ~isempty(k), defCfg.(fnd{k}) = cfg.(fnc{fi}); end;
	end;
end;
cfg = defCfg;

% load existing data:  
%	if base ws variable exists it has priority
%	if no ws variable but mat file exists reload variable from it
if evalin('base',sprintf('exist(''%s'',''var'')',vName)),
	v = evalin('base',vName);
elseif exist([vName,'.mat']) == 2,
	v = load([vName,'.mat']);
	v = v.(vName);
else,
	v = [];
end;
if lock, v = []; Tracker = []; end;

% verify overwrite if either exist
if ~isempty(v),
	if isunix, [s,r] = unix('osascript -e ''tell application "MATLAB" to activate'''); end;
	if strcmp(questdlg(sprintf('%s exists; overwrite its values?\n(%s_old backup variable will be created)',vName,vName), ... 
			'Verify...', 'Yes', 'No', 'Yes'), 'No'), return; end;

% make backup copy of original variable in ws
	assignin('base',[vName,'_old'],v);		% any previous variable renamed to <VNAME>_old
end;

% open the movie
fNameExt = fName;
[p,fName,e] = fileparts(fName);
try,
	info = dicominfo(fNameExt);
catch,
	info = [];
end;
if ~isempty(info),			% DICOM
	frameRate = info.RecommendedDisplayFrameRate;
	nFrames = info.NumberOfFrames;
	mh = fNameExt;
	audio = false;
else,						% movie
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
	if audio,
		try,
			[s,sr] = audioread(fNameExt);
		catch,
			audio = false;
		end;
	end;
end;

% parse the TextGrid (if no explicit keyFrames)
labs = []; isPraatInt = 0;
if isempty(keyFrames) & ~isempty(textgrid), 
	[keyFrames,labs] = ParseTextGrid(textgrid,frameRate); 
	isPraatInt = (size(keyFrames,2) == 2);	% Praat intervals
end;

% set up keyframes
if isempty(keyFrames), 
	frame = 1; ts = ''; enableS = 'off';
else,
	if ~isPraatInt,
		keyFrames = unique(keyFrames(:));		% format as [nExtents x first,last] array
		idx = find(diff([-1;keyFrames]) > 1);
		len = diff([idx;length(keyFrames)+1]);
		keyFrames = [keyFrames(idx) , keyFrames(idx)+len-1];
	end;
	frame = keyFrames(1,1); 
	if isempty(labs), 
		ts = '<keyframe>'; 
		labs = repmat({ts},size(keyFrames,1),1);
	else, 
		ts = labs{1}; 
	end; 
	enableS = 'on';
end;
if ~isempty(first) && first > 0 && first <= nFrames, frame = round(first); end;
kf = []; for k = 1 : size(keyFrames,1), kf = [kf , keyFrames(k,1):keyFrames(k,2)]; end;
isKF = ismember([1:nFrames], kf);		% true for keyframes

% load first image
p = []; a = {};
if ~isempty(ImageMod),
	if iscell(ImageMod),
		p = ImageMod{1};			% function
		a = ImageMod(2:end);		% args
	else,
		p = ImageMod;
	end;
	if ischar(p), p = str2func(p); end;		% coerce to pointer
end;
ImageModP = p; ImageModA = a;
try,
	img = GetMovieFrame(mh,frame,crop,reSize,p,a{:});
catch,
	error('unable to load frame %d in %s (%d frames available)', frame, fNameExt, nFrames);
end;
[rows,cols] = size(img);

% display image
pos = get(0, 'defaultfigureposition');
fh = figure('menubar','none', 'units','pixels', 'resize','off', 'tag','GETCONTOURS', 'position',[pos(1) , pos(2)+pos(4)-rows , cols, rows], 'visible','off');
set(fh, 'colormap', gray(256));
ih = imagesc(img);
imgH = gca;							% image axis
clim = get(imgH,'clim');

% gca userdata tracks cycling: 0 off, -1 back, +1 fwd
set(imgH, 'position',[0 0 1 1], 'xtick',[],'ytick',[],'userdata',0);

% make room for slider and audio
pos = get(fh,'position');
set(imgH,'units','pixels');
if audio, dy = 85; else, dy = 25; end;
set(fh,'position',[pos(1) pos(2)-dy pos(3) pos(4)+dy]);
set(imgH,'position',[1 dy+1 pos(3) pos(4)]);

% init frame slider controls
if ispc,
	fs = 10;
	pp = [.1,.4,10,1.2;10.6,.1,10,1.5;112,3,30,20;145,3,pos(3)-180,19;pos(3)-32,3,30,20];
else,
	fs = 12;
	pp = [.8,.25,6,1.2;7.4,.15,7.6,1.2;113 5 28 16;145,.1,pos(3)-180,19;pos(3)-32,5,28,16];
end;
uicontrol(fh, ...
	'style', 'text', ...
	'horizontalAlignment', 'right', ...
	'string', 'Frame:', ...
	'fontSize', fs, ...
	'units', 'characters', ...
	'backgroundColor', get(fh, 'color'), ...
	'foregroundColor', 'k', ...
	'position', pp(1,:));
frameH = uicontrol(fh, ...
	'style', 'edit', ...
	'horizontalAlignment', 'left', ...
	'units', 'characters', ...
	'string', num2str(frame), ...
	'fontSize', fs, ...
	'position', pp(2,:), ...
	'callback',{@GetContours,'FRAME','SPECIFY'});
uicontrol(fh, ...
	'style', 'pushbutton', ...
	'string', ' << ', ...
	'fontSize', 10, ...
	'units', 'pixels', ...
	'position', pp(3,:), ...
	'callback',{@GetContours,'CYCLE','BCK'});
scrollerH = uicontrol(fh, ...
	'style', 'slider', ...
	'min', 1, 'max', nFrames, ...
	'sliderStep', [1 50]/nFrames, ...
	'value', frame, ...
	'units', 'pixels', ...
	'position', pp(4,:), ...
	'callback',{@GetContours,'FRAME','SCROLL'});
uicontrol(fh, ...
	'style', 'pushbutton', ...
	'string', '>>', ...
	'fontSize', 10, ...
	'units', 'pixels', ...
	'position', pp(5,:), ...
	'callback',{@GetContours,'CYCLE','FWD'});

% display frame number
if isKF(frame), c = 'g'; else, c = 'y'; end;
if ispc, fs = 14; else, fs = 18; end;
th = text(50,80,sprintf('%s',FmtFrame(frame,frameRate)),'fontsize',fs,'color',c);

% display frame label
lh = text(50,100,ts, 'fontsize',fs,'color','c','interpreter','none');

% init menu
cbs = 'if strcmp(''on'',get(gcbo,''checked'')),set(gcbo,''checked'',''off'');else,set(gcbo,''checked'',''on'');end';
menu = uimenu(fh, 'label','GETCONTOURS');
uimenu(menu,'label','About...', ...
			'callback',{@GetContours,'ABOUT'});
uimenu(menu,'label','Configure...', ...
			'accelerator', '9', ...
			'callback',{@GetContours,'CONFIG'});
uimenu(menu,'label','Frame Information', ...
			'accelerator', 'I', ...
			'callback',{@GetContours,'INFO'});
uimenu(menu,'label','Export Contours...', ...
			'accelerator', 'E', ...
			'callback',{@GetContours,'EXPORT'});
uimenu(menu,'label','Find mm/pixel...', ...
			'callback',{@GetContours,'MPP'});
uimenu(menu,'label','Find Origin...', ...
			'callback',{@GetContours,'ORIGIN'});
if ~isempty(mpp) && ~isempty(origin), fce = 'on'; else, fce = 'off'; end;
fch = uimenu(menu,'label','Fourier Coefficients', ...
			'enable', fce, ...
			'callback',{@GetContours,'COEF'});

uimenu(menu,'label','Previous Frame', ...
			'separator','on', ...
			'accelerator', 'B', ...
			'busyAction', 'cancel', ...
			'callback',{@GetContours,'FRAME','PREV'});
uimenu(menu,'label','Next Frame', ...
			'accelerator', 'F', ...
			'busyAction', 'cancel', ...
			'callback',{@GetContours,'FRAME','NEXT'});
uimenu(menu,'label','Annotate Frame', ...
			'accelerator', 'A', ...
			'callback',{@GetContours,'ANNOTATE'});
uimenu(menu,'label','Prev Frame With Data', ...
			'busyAction', 'cancel', ...
			'callback',{@GetContours,'FRAME','DPREV'});
uimenu(menu,'label','Next Frame With Data', ...
			'busyAction', 'cancel', ...
			'callback',{@GetContours,'FRAME','DNEXT'});
uimenu(menu,'label','Previous Key Frame', ...
			'enable', enableS, ...
			'accelerator', '1', ...
			'busyAction', 'cancel', ...
			'callback',{@GetContours,'FRAME','KPREV'});
uimenu(menu,'label','Next Key Frame', ...
			'enable', enableS, ...
			'accelerator', '2', ...
			'busyAction', 'cancel', ...
			'callback',{@GetContours,'FRAME','KNEXT'});

uimenu(menu,'label','Stop Cycling', ...
			'separator','on', ...
			'accelerator', 'X', ...
			'callback',{@GetContours,'CYCLE','STOP'});
uimenu(menu,'label','Play @ Current Frame', ...
			'accelerator', 'P', ...
			'callback',{@GetContours,'PLAY'});

nh = uimenu(menu,'label','Inherit Anchors', ...
			'separator','on', ...
			'accelerator', '1', ...
			'callback',cbs);
uimenu(menu,'label','Redistribute Anchors', ...
			'accelerator','R', ...
			'callback',{@GetContours,'REDISTRIBUTE',1});
uimenu(menu,'label','Resample Anchors', ...
			'accelerator','8', ...
			'callback',{@GetContours,'RESAMPLE'});
uimenu(menu,'label','Delete Anchors', ...
			'accelerator','Y', ...
			'callback',{@GetContours,'DELETE'});
uimenu(menu,'label','Undo Last Change', ...
			'accelerator','Z',...
			'callback',{@GetContours,'UNDO',1});
if exist('make_snake') == 3, e = 'on'; else, e = 'off'; end;
uimenu(menu,'label','Apply Snake', ...
			'separator','on', ...
			'enable', e, ...
			'accelerator', 'N', ...
			'callback',{@GetContours,'SNAKE'});
uimenu(menu,'label','Image Forces', ...
			'accelerator', 'L', ...
			'callback',{@GetContours,'FILTER'});
dh = uimenu(menu,'label','Frame Differencing', ...
			'accelerator', 'D', ...
			'callback',{@GetContours,'DIFF'});
oh = uimenu(menu,'label','Reset Original Image', ...
			'accelerator', 'G', ...
			'callback',{@GetContours,'RESET'});
uimenu(menu,'label','Averaging...', ...
			'callback',{@GetContours,'AVERAGING'});

h = uimenu(menu,'label','Contrast','separator','on');
uimenu(h,'label','Full Range','callback',{@GetContours,'RANGE','FULL'});
uimenu(h,'label','Adjust Range','callback',{@GetContours,'RANGE','ADJUST'});

h = uimenu(menu,'label','Colormap');
gh = uimenu(h,'label','Gray','checked','on','callback',{@GetContours,'MAP','Gray'});
oh.UserData = gh;
uimenu(h,'label','Inv Gray','callback',{@GetContours,'MAP','Inv Gray'});
cMaps = {'Autumn','Bone','Colorcube','Cool','Copper','Flag','Hot','HSV','Jet','Lines','Pink','Prism','Spring','Summer','Winter'};
for k = 1 : length(cMaps),
	if k == 1, sep = 'on'; else, sep = 'off'; end;
	uimenu(h,'label',cMaps{k},'separator',sep,'callback',{@GetContours,'MAP',cMaps{k}});
end;

h = uimenu(menu,'label','Flip');
uimenu(h,'label','Horizontal', ...
			'userdata', th, ...
			'callback',{@GetContours,'FLIP','HORIZONTAL'});
uimenu(h,'label','Vertical', ...
			'callback',{@GetContours,'FLIP','VERTICAL'});

menu = uimenu(fh, 'label','TRACKING');			
e = 'on';
d = 'off';
if isempty(Tracker),
	n = '<None>'; e = 'off';
elseif ischar(Tracker), 
	n = Tracker;
	Tracker = str2func(Tracker);
else,
	n = func2str(Tracker);
end;
if exist(n),
	n = [upper(n), ':  configure'];
elseif ~isempty(Tracker),
	fprintf('Warning: tracking function %s not found\n', n);
	n = '<None>'; e = 'off';
	Tracker = [];
end;

tmh = uimenu(menu,'label',n, ...
			'enable',e, ...
			'callback',{@GetContours,'TRACK','CONFIG'});
uimenu(menu,'label','Clear', ...
			'separator','on', ...
			'callback',{@GetContours,'TRACK','CLEAR'});
if lock, es = 'off'; else, es = 'on'; end;
uimenu(menu,'label','Select...', ...
			'enable', es, ...
			'callback',{@GetContours,'TRACK','SELECT'});
tmh(2) = uimenu(menu,'label','Export Tracker Data...', ...
			'enable',e, ...
			'callback',{@GetContours,'TRACK','EXPORT'});
if isempty(keyFrames), ee = 'off'; else, ee = 'on'; end;
tmh(3) = uimenu(menu,'label','Track thru Sequence', ...
			'enable',e, ...
			'checked', ee, ...
			'accelerator','K', ...
			'callback',cbs);
tmh(4) = uimenu(menu,'label','Apply Tracking', ...
			'accelerator','T', ...
			'enable',e, ...
			'callback',{@GetContours,'TRACK','APPLY'});
		
% init audio display
if audio,
	p = [5 26 pos(3)-9 dy-29];
	axes('units','pixels','position',p);
	plot(s,'HitTest','off');
	q = max(abs(s),[],'all') * 1.1;
	set(gca,'xtick',[],'ytick',[],'xlim',[1 size(s,1)],'ylim',q*[-1 1],'ButtonDownFcn',@CursorAdjust,'userData',imgH);
	if ~isempty(keyFrames),
		d = diff(keyFrames,[],2);
		k = find(~d);		% points
		f = keyFrames(k,1);
		f = floor(((f-1)/frameRate)*sr)+1;
		if ~isempty(f), line([f,f]',get(gca,'ylim'),'color',[0 .6 0],'HitTest','off'); end;
		idx = find(d);		% intervals
		for k = 1 : length(idx),
			f = keyFrames(idx(k),:);
			x = floor(((f-1)/frameRate)*sr)+1;
			x = [x;x];
			x = x(:)';
			y = get(gca,'ylim');
			y = [y fliplr(y)];
			patch(x,y,[0 .6 0],'edgecolor',[0 .6 0],'facealpha',.1,'HitTest','off');
		end;
	end;
	f = floor(((frame-1)/frameRate)*sr) + 1;	
	ach = line([f;f],get(gca,'ylim'),'color','r','HitTest','off');	% cursor
	audio = struct('SRATE',sr, 'ACH',ach);			% package audio params here (srate, cursor handle)
else,
	audio = [];
end;

% init internal state
state = struct('IH', ih, ...				% image handle
				'IMGH', imgH, ...			% image axis
				'MH', mh, ...				% movie handle
				'TH', th, ...				% frame text handle
				'LH', lh, ...				% annotation text handle
				'FCH', fch, ...				% Fourier coefs menu handle
				'NH', nh, ...				% inherit anchors menu handle
				'DH', dh, ...				% frame differencing menu handle
				'TMH', tmh, ...				% tracking menu handles (1: name, 2: export, 3: track, 4: apply, 5: diag)
				'FRAMEH', frameH, ...		% frame box handle
				'SCROLLERH', scrollerH, ...	% scroller handle
				'CLIM', clim, ...			% original contrast limits
				'CONFIG', cfg, ...			% display configuration
				'CROP', crop, ...			% crop rect
				'RESIZE', reSize, ...		% resize factor
				'IMGMODP', ImageModP, ...	% image preprocessor 
				'IMGMODA', [], ...			% and its arguments
				'TRACKER', Tracker, ...		% automatic tracking procedure 
				'TPAR', [], ...				% and its parameters
				'TRKRES', [], ...			% and its result
				'CURFRAME', frame, ...		% currently displayed frame
				'KEYFRAMES', keyFrames, ...	% keyframe list
				'ISKF', isKF, ...			% true for keyframes
				'NFRAMES', nFrames, ...		% number of available frames
				'FRATE', frameRate, ...		% movie frame rate
				'LABELS', [], ...			% Praat keyframe labels
				'NPOINTS', nPoints, ...		% number of contour points
				'ANCHORS', [], ...			% current anchor points
				'NAPRESAMP', 11, ...		% resample to this number of anchor points
				'ALH', [], ...				% and their line handles
				'PREVAP', [], ...			% previous anchor points
				'XY', [], ...				% current contour points
				'CLH', [], ...				% and their line handle
				'USEAVG', 0, ...			% frame averaging enabled
				'AUDIO', audio, ...			% display audio params
				'AVG', 1, ...				% pre/post # frames to average
				'RH', [], ...				% contrast adjustment handle
				'FNAME', fName, ...			% movie name (w/o path and extension)
				'FNAMEEXT', fNameExt, ...	% full movie name
				'VNAME', vName, ...			% emit array name
				'MPP', mpp, ...				% mm per pixel factor
				'ORIGIN', origin, ...		% origin
				'LOCK', lock, ...			% true if editing permitted
				'DEFPAR', [], ...			% defaults
				'PARAMS', []);				% parameters
state.LABELS = labs;
state.IMGMODA = ImageModA;

% parameter defaults
dPar = struct('FNAME', vName, ...				% export filename
				'NPOINTS', 100, ...				% # of contour points
				'PLAYINT', [-500 500], ...		% play interval around current frame (ms)
				'MPP', mpp, ...					% mm per pixel factor
				'ORIGIN', origin, ...			% origin
				'SIGMA', 5, ...					% image forces Gaussian sigma
				'SAVEIMG', 0, ...				% save images
				'DELTA', 2, ...					% snake params:  delta
				'BPEN', 2, ...					%   band penalty
				'ALPHA', .8, ...				%   alpha
				'LAMBDA', .95, ...				%   lambda1
				'NSPTS', 39);					% # snake estimation points
state.DEFPAR = dPar;
state.PARAMS = dPar;

n = sprintf('%s  [%d of %d]',fName,frame,nFrames);
if ~lock, n = [n , ' (editable)']; end;
set(fh,'name', n, ...
		'closeRequestFcn',{@GetContours,'CLOSE'}, ...
		'tag','GETCONTOURS', ...
		'visible','on', ...
		'userData',state);
if ~lock, set(ih,'buttonDownFcn', {@GetContours,'DOWN'}); end		% trap buttonDown on image
axes(imgH);

% initialize ws state variable
if isempty(v) || clean,		% no pre-existing variable
	v = struct('XY',state.XY,'ANCHORS',[],'FRAME',frame,'NOTE','','TRKRES',[]);
else,
	if isstruct(v) && isfield(v,'ANCHORS'),
		if isfield(v,'TIME'), v = rmfield(v,{'TIME','IMAGE'}); end;
		if ~isfield(v,'TRKRES'), v(1).TRKRES = []; end;		% bw compatibility
		k = find(frame == cell2mat({v.FRAME}));
		if isempty(k),
			v(end+1) = struct('XY',[],'ANCHORS',[],'FRAME',frame,'NOTE','','TRKRES',[]);
		else,
			anchors = v(k).ANCHORS;
		end;
	end;	
end;
if ~lock, assignin('base',vName,v); end;

% configure tracker if any
if ~isempty(Tracker), 
	state.TPAR = Tracker('CONFIG', state, 1); 
	delete(findobj(gca,'tag','CONTOUR'));			% delete any accidentally created points
	state.CLH = []; state.ALH = []; state.ANCHORS = [];
end;

% initialize contour
if ~isempty(anchors),
	state.ANCHORS = anchors;
	state = UpdateContour(state); 
end;

% set internal state
set(fh, 'userdata',state);


%===============================================================================
% ISPC  - returns true for PCWIN versions of Matlab

function result = ispc

result = strncmp(computer,'PC',2);


%===============================================================================
% LINEINTERSECT  - returns intersection of lines determined by A and B

function [intersectPt, cross] = LineIntersect(A, B)

X = 1; Y = 2;
FixRoundOffErr = inline('round(v*1e10)*1e-10', 'v');

denom = A(1,X) * (B(2,Y) - B(1,Y)) + ...
		A(2,X) * (B(1,Y) - B(2,Y)) + ...
		B(2,X) * (A(2,Y) - A(1,Y)) + ...
		B(1,X) * (A(1,Y) - A(2,Y));

intersectPt = [];
cross = [];
if ~denom, return; end;				% lines are coincident or parallel

s = FixRoundOffErr( (	A(1,X) * (B(2,Y) - B(1,Y)) + ...
						B(1,X) * (A(1,Y) - B(2,Y)) + ...
						B(2,X) * (B(1,Y) - A(1,Y))) / denom );
	  
t = FixRoundOffErr(-( 	A(1,X) * (B(1,Y) - A(2,Y)) + ...
						A(2,X) * (A(1,Y) - B(1,Y)) + ...
						B(1,X) * (A(2,Y) - A(1,Y))) / denom);
		
intersectPt = FixRoundOffErr([(A(1,X) + s * (A(2,X) - A(1,X))) (A(1,Y) + s * (A(2,Y) - A(1,Y)))]);
cross = (0<=s & s<=1 & 0<=t & t<=1);


%===============================================================================
% MAKEGRID  - build semi-polar grid

function [glx,gly] = MakeGrid(N, len, res, flp)

% find polar gridlines with RES separation at GL/2
dth = 2*asin(res/len);
th = linspace(0,pi,pi/dth);
[glx,gly] = pol2cart(th,len);
glx(2,:) = zeros(1,length(glx));
gly(2,:) = zeros(1,length(gly));

% find vertical gridlines; adjust vertical resolution to distance between polar line midpoints
adjRes = mode(sqrt(sum(diff([mean(glx)' , mean(gly)']).^2,2)));
dy = -adjRes : -adjRes : -len;
gly = [gly , [dy;dy]];
glx = [glx , [-ones(1,length(dy))*len ; zeros(1,length(dy))]];

if flp, glx = -glx; end;


%===============================================================================
% MAKEPOINT  - create anchor point

function lh = MakePoint(x,y,cfg)

if ismac, ms = 25; else, ms = 15; end;
lh = line(x,y,'marker','.','markerSize',cfg.DOTSIZE,'color',cfg.DOTCOLOR,'tag','CONTOUR','buttonDownFcn',{@GetContours,'DOWN','POINT'});


%===============================================================================
% NEWFRAME  - initialize new frame

function state = NewFrame(state, newFrame, forceInherit, external)

if nargin < 4, external = 0; end;

% update output variable state before change
if state.LOCK,
	frames = [];
else,
	v = evalin('base',state.VNAME);			% frame data in base ws (VNAME)
	frames = cell2mat({v.FRAME});			% frames with data
	k = find(state.CURFRAME == frames);
	v(k).XY = state.XY;
	v(k).ANCHORS = state.ANCHORS;
	v(k).NOTE = get(state.LH,'string');
	v(k).TRKRES = state.TRKRES;
	assignin('base',state.VNAME,v);
end;

% change frame
state.CURFRAME = newFrame;

% update audio cursor
if ~isempty(state.AUDIO),
	f = floor(((newFrame-1)/state.FRATE)*state.AUDIO.SRATE) + 1;
	set(state.AUDIO.ACH,'xdata',[f;f]);
end;

% get current annotation if any
ts = '';
k = find(newFrame == frames);
if ~isempty(k) && ~state.LOCK, ts = v(k).NOTE; end;

% use Praat labels if no annotation
if isempty(ts) && ~isempty(state.KEYFRAMES),
	k = find(newFrame>=state.KEYFRAMES(:,1) & newFrame<=state.KEYFRAMES(:,2));
	if ~isempty(k), k = k(end); ts = state.LABELS{k}; end;
end;

% update frame
set(state.IH,'cdata',GetImage(state));				% image
ah = gca;
if external, axes(get(state.IH,'parent')); end;
if state.ISKF(newFrame), c = 'g'; else, c = 'y'; end;
set(state.TH,'color',c,'string',sprintf('%s',FmtFrame(newFrame,state.FRATE)));	% frame count, time
set(state.LH,'string',ts);							% annotation
n = sprintf('%s  [%d of %d]',state.FNAME,newFrame,state.NFRAMES);
if ~state.LOCK, n = [n , ' (editable)']; end;
set(gcbf,'name',n);
set(state.FRAMEH, 'string', num2str(newFrame));
set(state.SCROLLERH, 'value', newFrame);
if state.LOCK, return; end;

delete(findobj(gca,'tag','CONTOUR'));				% delete existing handles
state.ALH = []; state.CLH = [];
k = find(newFrame == frames);
if isempty(k),					% virgin frame
	if strcmpi('off',get(state.NH,'checked')) && ~forceInherit,
		state.ANCHORS = [];		% don't inherit anchors from last frame
		state.XY = [];
		state.TRKRES = [];
	end;
	vv = struct('XY',state.XY,'ANCHORS',state.ANCHORS,'FRAME',newFrame,'NOTE',ts,'TRKRES',[]);
	v(end+1) = vv;
else,							% update from existing anchors
	if forceInherit || (isempty(v(k).ANCHORS) && strcmpi('on',get(state.NH,'checked'))),
		v(k).XY = state.XY;		% inherit into empty frame
		v(k).ANCHORS = state.ANCHORS;
		state.PREVAP = [];
	else,
		state.PREVAP = state.ANCHORS;
	end;			
	state.ANCHORS = v(k).ANCHORS;
	state.XY = v(k).XY;
	v(k).NOTE = ts;
	state.TRKRES = v(k).TRKRES;
end;
assignin('base',state.VNAME,v);

% create new handles
for k = 1 : size(state.ANCHORS,1),
	state.ALH(k) = MakePoint(state.ANCHORS(k,1),state.ANCHORS(k,2),state.CONFIG);
end;
if ~isempty(state.XY),
	state.CLH = line(state.XY(:,1),state.XY(:,2),'color',state.CONFIG.LINECOLOR,'linewidth',state.CONFIG.LINEWIDTH,'tag','CONTOUR','hitTest','off');
	uistack(state.CLH,'bottom'); uistack(state.CLH,'up');

% update contour line using active Tracker handler
	if ~isempty(state.TRACKER),
		ns = state.TRACKER('PLOT', state);
		if ~isempty(ns), state = ns; end;
	end;

end;
if external, axes(ah); end;


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
	labs(k) = [];
end


%===============================================================================
% SAVEVAR  - save variable as VNAME to specified MAT file

function SaveVar(var, fName)

vName = fName;
vs = vName;
eval([vName '=var;save ' fName ' ' vs]);
fprintf('wrote %s.mat\n', fName);


%===============================================================================
% SHAPETODIST  - map tongue contour to unwrapped VT distance function

function [d,xyi] = ShapeToDist(xy, glx, gly)

res = 100;				% comparison resolution (nPoints)

% resample contour to 100 equally spaced points
k = [0 ; cumsum(sqrt(sum(diff(xy).^2,2)))];
xy = interp1(k,xy,linspace(0,k(end),res),'pchip');

% find intersection of each gridline with contour
nLines = size(glx,2);
xyi = NaN(nLines,2);
d = NaN(nLines,1);
n = ones(1,res);
for gi = 1 : nLines,
	x = linspace(glx(1,gi),glx(2,gi),res);
	y = interp1(glx(:,gi),gly(:,gi),x);
	
% distance between every point of contour (rows) with every point of interpolated gridline (cols)
	dist = sqrt((xy(:,1)*n-n'*x).^2+(xy(:,2)*n-n'*y).^2);
	[~,k] = min(dist(:));
	[r,c] = ind2sub(size(dist),k);			% minimum distance contour point

% indices of bracketing contour point
	if r == 1,
		idx = [1 2];
	elseif r == res,
		idx = [res-1 res];
	else,
		idx = [r-1 r+1];
	end;

% compute intersection
	[p,cross] = LineIntersect(xy(idx,:),[glx(:,gi),gly(:,gi)]);
	if ~cross, continue; end;			% no intersection with this gridline
	xyi(gi,:) = p;

% compute distance along line
	if gly(1,gi) && diff(gly(:,gi)) == 0,	% pharyngeal section
		d(gi) = abs(p(1));
	else,
		d(gi) = sqrt(sum(p.^2));			% polar section
	end;

end;


%===============================================================================
% UPDATECONTOUR  - update displayed contour
%
% non-empty MODE passed to tracker plot fn on move

function state = UpdateContour(state, mode)

if nargin < 2, mode = []; end;

% create anchor points
if isempty(state.ALH) && ~isempty(state.ANCHORS),
	for k = 1 : size(state.ANCHORS,1),
		state.ALH(k) = MakePoint(state.ANCHORS(k,1),state.ANCHORS(k,2),state.CONFIG);
	end;
end;

% create contour line
if isempty(state.CLH) && ~isempty(state.ANCHORS),
	state.CLH = line(state.ANCHORS(1,1),state.ANCHORS(1,2),'color',state.CONFIG.LINECOLOR,'linewidth',state.CONFIG.LINEWIDTH,'tag','CONTOUR','hitTest','off');
	if ~verLessThan('matlab','8.5.0'), set(state.CLH,'PickableParts','none'); end;
	uistack(state.CLH,'bottom'); uistack(state.CLH,'up');
end;

% update contour line using active Tracker handler
if ~isempty(state.TRACKER),
	ns = state.TRACKER('PLOT', state, mode);
	if ~isempty(ns),
		state = ns;
		return;
	end;
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

