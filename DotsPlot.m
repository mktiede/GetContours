function DotsPlot(fName, p, frames)
%DOTSPLOT - cycle through image sequence displaying tracked dots fit
%
%	usage:  DotsPlot(fName, p, frames)
%
% this procedure cycles through the image sequence specified by movie FNAME
% displaying the fit of the tracked features specified by P [nFrames x nPoints]
% as returned by DOTSTRACK (if specified)
%
% optional FRAMES defaults to [1:nFrames]
%
% see also DOTSPLACE, DOTSTRACK

% mkt 07/15
% mkt 10/15 support for names and structured points arrays

if nargin < 1, eval('help DotsPlot'); return; end;
if nargin < 2, p = []; end;
if nargin < 3, frames = []; end;

% handle toggle button status change
if nargin == 3 && ischar(frames) && strcmp(frames,'TOGGLE'),
	set(gcbf,'userData',get(p.NewValue,'userData'));
	uiresume;
	return;
end;

% open the movie
try, mh = VideoReader(fName); catch, error('unable to open movie file %s', fName); end;
if verLessThan('matlab','8.5.0'),
	nFrames = get(mh,'NumberOfFrames');
	frameRate = get(mh,'FrameRate');
else,
	frameRate = mh.FrameRate;
	nFrames = floor(mh.Duration*frameRate);
end;

% validate frames
if isempty(p),
	if isempty(frames), frames = [1:nFrames]; end;
else,
	pFrames = cell2mat({p(:,1).FRAME});
	if isempty(frames),
		frames = pFrames;
	else,
		unknownFrames = setdiff(frames,pFrames);
		if ~isempty(unknownFrames),
			unknownFrames
			error('FRAMES not found in P');
		end;
	end;
end;

% display the first image
img = GetMovieFrame(mh,frames(1));
[fh,ih] = implot(img);
set(fh,'userData',0,'name',fName);

% display the markers
colors = 'grkm';		% good, bad track, invalid, bad+invalid
if ~isempty(p), 
	fIdx = find(frames(1) == pFrames);
	xy = cell2mat({p(fIdx,:).POS}');
	labels = {p(fIdx,:).LABEL};
	for k = 1 : length(labels),
		if isempty(p(fIdx,k).STATUS), c = 1; else, c = p(fIdx,k).STATUS+1; end;
		lh(k) = line(xy(k,1),xy(k,2),'color',colors(c),'marker','+');
		th(k) = text(xy(k,1),xy(k,2),['  ',labels{k}],'color','w');
	end;
end;

% add controls
bh = uibuttongroup(fh,'units','pixels','position',[5 5 104 20]);
uicontrol(bh,'style','togglebutton','position',[2 2 18 15],'string','<<','userdata',-3);
uicontrol(bh,'style','pushbutton','position',[22 2 18 15],'string','<','callback','set(gcbf,''userdata'',-1);uiresume');
b = uicontrol(bh,'style','togglebutton','position',[42 2 18 15],'string','O','userdata',0);
uicontrol(bh,'style','pushbutton','position',[62 2 18 15],'string','>','callback','set(gcbf,''userdata'',1);uiresume');
uicontrol(bh,'style','togglebutton','position',[82 2 18 15],'string','>>','userdata',3);
set(bh,'SelectedObject',b,'SelectionChangeFcn',{@DotsPlot,'TOGGLE'});

uicontrol(fh,'style','edit','position',[5 28 50 18],'horizontalAlignment','right','enable','inactive','string','Frame: ');
eh = uicontrol(fh,'style','edit','position',[54 28 54 18],'horizontalAlignment','left', ...
			'string',sprintf(' %04d',frames(1)),'callback','set(gcbf,''userdata'',2);uiresume');

% display loop
idx = 1;
while 1,
	state = get(fh,'userData');
	switch state,
		case {-3,-1}, idx = idx - 1; if idx < 1, idx = length(frames); end;
		case 0,
		case {1,3}, idx = idx + 1; if idx > length(frames), idx = 1; end;
		case 2,		% frame field
			f = round(str2num(get(eh,'string')));
			if isempty(f), f = frames(idx); end;
			[~,idx] = min(abs(frames-f));
	end;

	img = GetMovieFrame(mh,frames(idx));
	set(ih,'cdata',img);
	set(eh,'string',sprintf(' %04d',frames(idx)));
	set(fh,'name',sprintf('%s  FRAME %04d / %04d  (%.3f secs)', fName, frames(idx), length(frames), (frames(idx)-1)/frameRate));
	if ~isempty(p), 
		fIdx = find(frames(idx) == pFrames);
		xy = cell2mat({p(fIdx,:).POS}');
		for k = 1 : size(xy,1), 
			if isempty(p(fIdx,k).STATUS), c = 1; else, c = p(fIdx,k).STATUS+1; end;
			set(lh(k),'XData',xy(k,1),'YData',xy(k,2),'color',colors(c));
			set(th(k),'position',xy(k,:)); 
		end;
	end;
	
	if abs(state) < 3,
		set(fh,'userdata',0);
		set(bh,'SelectedObject',b);
		uiwait;
	else,
		drawnow;
	end;
	if ~ishandle(fh), break; end;
end;
