function [p,ih] = DotsPlace(src, varargin)
%DOTSPLACE  - place tracking locations on image
%
%	usage:  [p,ih] = DotsPlace(src, ...)
%
% use this procedure in tandem with DOTSTRACK to track image features of interest
% through an image sequence
%
% for a new set of points:
%   click to add a new point (origin is ULC)
%   click and drag an existing point to change its position
%   right click on an existing point to change its name
%   double-click on an existing point to delete
%
% for an existing set of points:
%   click and drag on the existing points to adjust their location
%   red points have failed tracking threshold
%
% the procedure blocks the command line until the "CONTINUE" button is clicked
%
% image SRC can be any one of the following:
%   movie filename (e.g. 'foo.mov')
%   movieHandle specifying an open movie handle
%   [nRows x nCols (xRGB)] image array
%
% optional supported 'NAME',VALUE pairs:
%   EDIT   - allow editing of existing point array if 'T' (default 'F')
%   IH     - image handle to reuse (default opens a new figure sized to SRC)
%   FRAME  - 1-based frame number to load from a movie sequence (default 1)
%   P0     - array of existing points ([1 x nPoints]; default none);
%   XFORM  - conditioning function applied to each frame (e.g., @(img) 255-uint8(mean(img,3)))
%
% returns updated points array
% optionally returns image handle IH
%
% points array is an array-of-structs [1 x nPoints] with fields
%   POS    - [X,Y] (pixel coordinates, relative to ULC origin)
%   LABEL  - point id string
%   FRAME  - 1-based frame number of fitted image
%   TIME   - 0-based offset of frame from beginning of movie in seconds
%   STATUS - 0: valid, 1: failed tracking step; 2: user flagged invalid point (tracking skipped)
%   CONF   - tracking confidence value
%
% see also DOTSTRACK, DOTSPLOT

% mkt 07/15
% mkt 10/15 support for names and structured points arrays

%----- branch by action

if nargin > 2 && ~ischar(varargin{1}),
	action = varargin{2};
else,
	action = 'INIT';
end;

switch action,

%----- DOWN  - mouse down handler

	case 'DOWN',
		newPts = get(gca,'UserData');		% nonzero:  permit creation, deletion

% clicked on image:  make new point
		if strcmpi(get(src,'Type'),'image'),	
			if ~newPts, return; end;		% if only point moves allowed
			xy = get(gca,'CurrentPoint');		
			xy = floor(xy(1,1:2));
			lh = line(xy(1),xy(2),'marker','+','color','g','tag','DPLC','userdata',[],'buttonDownFcn',{@DotsPlace,'DOWN'});
			
% clicked on existing point
		else,
			if strcmp(get(src,'type'),'text'),
				th = src;
				lh = get(th,'userdata');
			else,
				lh = src;
				th = get(lh,'userdata');
			end;
			switch get(gcbf,'SelectionType'),
				case 'open',					% double-click
					if newPts, 
						if ~isempty(th), delete(th); end;
						delete(lh);								% delete existing point
					end;
					return;
				case 'normal',					% no mod -- fall thru
				otherwise,						% change label name
					if newPts,
						name = get(th,'string');
						name(1:2) = [];
						set(th,'string',['  ',GetName(name)]);
						return;
					end;
			end;
		end;
		set(gcbf, ...						% set motion handlers
			'WindowButtonMotionFcn', {@DotsPlace,'MOVE',lh}, ...	
			'WindowButtonUpFcn', {@DotsPlace,'UP',lh}, ...	
			'Pointer', 'crosshair');


%----- INITialze

	case 'INIT',
		if nargin < 1, eval('help DotsPlace'); return; end

% parse args
		ih = []; frame = []; p0 = []; allowEdit = 0; xform = [];
		for ai = 2 : 2 : length(varargin),
			switch upper(varargin{ai-1}),
				case 'EDIT', allowEdit = strcmpi(varargin{ai}(1),'T');
				case 'IH', ih = varargin{ai};
				case 'FRAME', frame = varargin{ai};
				case 'P0', p0 = varargin{ai};
				case 'XFORM', xform = varargin{ai};
				otherwise, error('unrecognized parameter (%s)', varargin{ai-1});
			end;
		end;
		if isempty(frame), frame = 1; end;

% load and display image
		if ischar(src),				% movie filename
			try, mh = VideoReader(src); catch, error('unable to open movie file %s', src); end;
			img = GetMovieFrame(mh, frame);
			if ~isempty(xform), f = xform{1}; a = xform(2:end); img = f(img,a{:}); end;
			t = (frame-1)/mh.FrameRate;
		elseif ishandle(src),		% open movie handle
			img = GetMovieFrame(src, frame); 
			if ~isempty(xform), f = xform{1}; a = xform(2:end); img = f(img,a{:}); end;
			t = (frame-1)/mh.FrameRate;
		else,						% explicit image
			img = src;
			if ~isempty(xform), f = xform{1}; a = xform(2:end); img = f(img,a{:}); end;
			t = 0;
		end;

		if isempty(ih),				% existing image handle
			[~,ih] = implot(img);
		else,						% construct new figure
			set(ih,'cdata',img);
		end;
		ah = get(ih,'Parent');		% "get" provides backwards compatibility
		fh = get(ah,'Parent');

% add buttons
		bh = uicontrol(fh,'style','pushbutton','position',[5 5 60 20],'String','CONTINUE','callback','uiresume(gcbf)');

% init buttonDown event handling
		if isempty(p0),				% no previous points:  allow create, move, delete, rename
			allowEdit = 1;
			set(ih,'buttonDownFcn',{@DotsPlace,'DOWN'});
			set(ah,'userData',1);			% flag point creation, deletion allowed
			bh(2) = uicontrol(fh,'style','pushbutton','position',[5 25 60 20],'String','CLRALL','callback','delete(findobj(gca,''tag'',''DPLC''))');
		else,						% previous points
			set(ah,'userData',allowEdit);	% flag mvt only unless editing permitted
			if allowEdit, set(ih,'buttonDownFcn',{@DotsPlace,'DOWN'}); end;
			col = 'grkm';					% good, bad track, invalid, bad+invalid
			for li = 1 : length(p0),
				if isempty(p0(li).STATUS), c = 1; else, c = p0(li).STATUS+1; end;
				lh = line(p0(li).POS(1),p0(li).POS(2),'color',col(c),'marker','+','tag','DPLC','buttonDownFcn',{@DotsPlace,'DOWN'});
				if allowEdit, s = sprintf('  %s',p0(li).LABEL); else, s = sprintf('  %d %s',li,p0(li).LABEL); end;
				th = text(p0(li).POS(1),p0(li).POS(2),s,'color','w','userdata',lh,'tag','DPLC','interpreter','none','buttonDownFcn',{@DotsPlace,'DOWN'});
				set(lh,'userdata',th);
			end;
		end;

% block command line until "CONTINUE" clicked
		uiwait(fh);
		if ~ishandle(fh), p = p0; return; end;		% closing window returns original points

% update points
		lh = flipud(findobj(gca,'marker','+','tag','DPLC'));
		if length(lh) && (isempty(p0) || allowEdit),
			p(1,length(lh)) = struct('POS',[],'LABEL',[],'FRAME',[],'TIME',t,'STATUS',[],'CONF',[]);
		else,
			p = p0;
		end;
		for li = 1 : length(lh),
			x = get(lh(li),'XData');
			y = get(lh(li),'YData');
			p(li).POS = [x,y];
			if allowEdit,
				s = get(get(lh(li),'userdata'),'string'); 
				if length(s)>2, s(1:2) = []; end;
				p(li).LABEL = s;
			end;
		end;

% clean up and exit
		delete(bh);
		delete(findobj(gca,'tag','DPLC'));
		return;

%----- MOVE  - mouse mvt handler

	case 'MOVE',
		xy = get(gca, 'CurrentPoint');
		xy = floor(xy(1,1:2));
		lh = varargin{3};
		if ~ishandle(lh), return; end;
		x = get(lh,'XData');
		y = get(lh,'YData');
		[~,k] = min(sqrt((x-xy(1)).^2+(y-xy(2)).^2));
		if x(k) ~= xy(1) || y(k) ~= xy(2), 
			x(k) = xy(1); y(k) = xy(2);
			set(lh, 'XData',x, 'YData',y);
			th = get(lh,'userdata');
			if ~isempty(th), set(th,'position',[x y 0]); end;
		end;
	
%----- UP  - mouse up handler

	case 'UP',
		set(gcbf, ...						% clear motion handlers
			'WindowButtonMotionFcn', '', ...
			'WindowButtonUpFcn', '', ...
			'Pointer', 'arrow');
		lh = varargin{3};
		th = get(lh,'userdata');
		if get(gca,'UserData') && isempty(th),		% add label if new point
			x = get(lh,'XData'); y = get(lh,'YData');
			name = GetName('');
			th = text(x,y,['  ',name],'color','w','tag','DPLC','userdata',lh,'buttonDownFcn',{@DotsPlace,'DOWN'});
			set(lh,'userdata',th);
		end;

%----- error

	otherwise,
		error('DOTSPLACE:  unknown parameter (%s)', action);
		
end;


%=============================================================================
% GETNAME  - edit dot name

function name = GetName(name)

width = 250;
height = 100;
ts = 'Enter Point Name';

fh = gcf;
units = get(fh, 'units');
if strcmpi(units,'pixels'),
	units = [];
else,
	set(fh, 'units','pixels');
end;	
fPos = get(fh, 'position');
if ~isempty(units), set(fh, 'units',units); end;

pos = [fPos(1)+(fPos(3)-width)/2 , fPos(2)+(fPos(4)-height)/2 , width , height];

cfg = dialog('Name', ts, ...
	'menubar', 'none', ...
	'Position', pos, ...
	'UserData', 0);

% name field
eh = uicontrol(cfg, ...
	'Position', [20 60 width-40 25], ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', name);

% OK, cancel buttons
uicontrol(cfg, ...		% buttons
	'Position',[width/2-70 15 60 25], ...
	'String','OK', ...
	'Callback','set(gcbf,''UserData'',3);uiresume');
uicontrol(cfg, ...
	'Position',[width/2+10 15 60 25], ...
	'String','Cancel', ...
	'Callback','uiresume');

% wait for input
drawnow;
uicontrol(eh);
uiwait(cfg);
if get(cfg, 'UserData'),
	name = get(eh, 'string');
end;
delete(cfg);


