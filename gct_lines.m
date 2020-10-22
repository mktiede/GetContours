function r = gct_lines(action, state, varargin)
%GCT_LINES  - GetContours tracker that follows intersection of lines with contour
%
%	usage:  r = gct_lines(action, state, ...)
%
% this procedure assumes that paired anchors define line segments intersecting 
% the visible contour and tracks the intensity shift along them (from green to red anchor pts)
%
% configuration dialog parameters:
%   Algorithm  - one of "Overall Maximum" (highest intensity value along line), 
%                "First Max above Thresh" (highest intensity value in first values
%                exceeding current Threshold), or "First Edge above Thresh" (first
%                intensity value found above current Threshold)
%   Threshold  - value used to delimit multiple intensity peaks along line;
%                entering an empty value causes the current intensity maximum along 
%                the last defined line segment to be displayed
%   Plot       - plots the current intensity values along the last defined line
%   Med. Filt. - applies a 3rd order median filter to the intensity values before detection
%   Stretch    - stretches contrast to range along line before detection
%   Wt by Last - applies a weighting function determined by the last frame's intersection point
%   Window Len - length of weighting window
%
% the intersection for a given anchor line pair is stored as a marker whose handle is
% stored in the line connecting that pair; clicking on that line allows updating its position
%
% data saved to workspace have format:
% [nFrames x X,Y, % along green to red line of intersection, intensity @intersection, frame# x nAnchorPairs]
%
% see also GCT_FAN

% mkt 06/18
% mkt 07/18 fix empty TRKVAL save-to-ws issue  
% mkt 12/19 facelift & fixes
% mkt 10/20 UltraFest IX release

% trap line click callback
if ~ischar(action),
	HandleClick(action, state);
	return;
end;

%	branch by ACTION

switch upper(action),
		
%-----------------------------------------------------------------------------
% ADDPT:  add new point flag (called by GetContours:DOWN handler)
%
%	returns 1 indicating new point should always be appended to anchor list

	case 'ADDPT',
		r = 1;


%-----------------------------------------------------------------------------
% CONFIG:  handle configuration
%
% 	returns TPAR = new tracker params
%
%	varargin{1} set to 1 on first entry; 0 on subsequent entries
%
% r = [] flags user cancelled

	case 'CONFIG',
		r = DoConfig(mfilename, state, varargin{1});


%-----------------------------------------------------------------------------
% DIAGS:  show diagnostics
%
% 	return ignored

	case 'DIAGS',
		PlotStripe(state);


%-----------------------------------------------------------------------------
% EXPORT:  output data to file in tab-delimited format

	case 'EXPORT',
		v = evalin('base',state.VNAME);			% frame data in base ws (VNAME)
		f = cell2mat({v.FRAME});				% frames with data
		[f,k] = sort(f);
		v = v(k);								% sorted	

% open the file (append to existing)
		fName = state.TPAR.FNAME;
		if strcmp(fName,'<CMD_LINE>'),
			isNew = 1;
			fid = 1;
		else,
			isNew = ~exist(fName,'file');
			fid = fopen(fName, 'at');
			if fid == -1,
				error('error attempting to open %s', fName);
			end;
		end;

% write headers if necessary
		if isNew,
			fprintf(fid, 'FRAME\tNOTE\tANCHOR\tX\tY\tPCT\tINT\n');
		end;
		
% append data
		for k = 1 : length(v),
			d = v(k).TRKRES;
			for a = 1 : size(d,1),
				fprintf('%d\t%s\t%d\t%.0f\t%.0f\t%.3f\t%.3f\n', f(k), v(k).NOTE, a, d(a,:));
			end;
		end;

% clean up
		if fid > 1, 
			fclose(fid);
			fprintf('\nTracking data from %s appended to %s\n', state.VNAME, fName);
		end;
		
		
%-----------------------------------------------------------------------------
% PLOT:  post-tracking plotting handler (called after CLH creation)
%
% overrides contour plotting to obtain line segments mapped onto anchors
%
%	returns updated STATE variable 

	case 'PLOT',
		nAnchors = size(state.ANCHORS,1);
		xy = state.ANCHORS;
		if mod(nAnchors,2),
			set(state.ALH(end),'color','c');
			xy(end,:) = []; 				% drop unpaired
			nAnchors = nAnchors - 1; 
		end;	
		if nAnchors < 2,
			r = [];
			return;
		end;

% modify anchor markers such that first in pair is green and second is red
		h = reshape(state.ALH(1:nAnchors),2,nAnchors/2)';
		set(h(:,1),'color','g');
		set(h(:,2),'color','r');

% plot line segments between them:  the last tracked intersection is stored as the userdata of each line
		xx = reshape(xy(:,1),2,nAnchors/2)';
		yy = reshape(xy(:,2),2,nAnchors/2)';
		for k = 1 : length(state.CLH), delete(state.CLH(k).UserData); end;
		delete(state.CLH);
		for k = 1 : size(xx,1),
			lh(k) = line(xx(k,:),yy(k,:),'color',state.CONFIG.LINECOLOR,'linewidth',state.CONFIG.LINEWIDTH,'tag','CONTOUR','buttonDownFcn',@gct_lines);
			uistack(lh(k),'bottom'); uistack(lh(k),'up');
			if ~isempty(state.TRKRES) && size(state.TRKRES,1)>=k,
				lh(k).UserData = line(state.TRKRES(k,1),state.TRKRES(k,2),'color',state.CONFIG.LINECOLOR, 'marker','o', 'linestyle','none', ...
							'markersize',15, 'hitTest','off', 'tag','CONTOUR');				
				if ~verLessThan('matlab','8.5.0'), set(lh(k).UserData,'PickableParts','none'); end;
			end;
		end;
		state.CLH = lh;

% define "contour"
		xx = [xx , NaN(size(xx,1),1)]';
		yy = [yy , NaN(size(yy,1),1)]';
		state.XY = [xx(:) , yy(:)];
		
		r = state;
		
		
%-----------------------------------------------------------------------------
% SAVE:  save tracker data as variable in base ws
%
%	returns nothing
%
% this handler saves TRKRES data as VNAME_trk, formatted as [nFrames x X,Y,%,I,frame x nAnchorPairs]
% tracking data:  [nAnchors x ,Y, percent of line segment & intensity at intersection, frame#]
% assumes consistent number of anchors each frame

	case 'SAVE',
		v = evalin('base',state.VNAME);			% frame data in base ws (VNAME)
		f = cell2mat({v.FRAME});				% frames with data
		[f,k] = sort(f);
		v = v(k);								% sorted	

		nAP = floor(size(state.ANCHORS,1)/2);
		d = {v.TRKRES};
		if nAP < 1 || isempty(cell2mat(d)), return; end;
		idx = find(cellfun(@isempty,d));
		for k = idx, d{k} = NaN(nAP,5); end;		% fill untracked frames
		d = permute(reshape(cell2mat(d')',[5, nAP, length(v)]),[3 1 2]);
		
		vName = sprintf('%s_trk', state.VNAME);
		if evalin('base',sprintf('exist(''%s'',''var'')',vName)),
			if strcmp(questdlg(sprintf('%s exists; overwrite it?',vName), 'Verify...', 'Yes', 'No', 'Yes'), 'No'), return; end;
		end;
		
		assignin('base', vName, d);
		fprintf('%s created in base workspace\n', vName);
		
		
%-----------------------------------------------------------------------------
% TRACK:  track current frame handler
%
%	returns updated STATE variable on success, [] on failure
%
% 	varargin{1} is nonzero on first track of a sequence
%
% state.TRKRES set to [nAnchors x X,Y, percent of line & intensity at intersection, frame]

	case 'TRACK',
		nAnchors = size(state.ANCHORS,1);
		if nAnchors < 2,
			r = [];
			return;
		end;
		xy = state.ANCHORS;
		if mod(nAnchors,2), xy(end,:) = []; nAnchors = nAnchors - 1; end;
		img = GetContours('GETMOVIEFRAME',state.CURFRAME);
		xx = reshape(xy(:,1),2,nAnchors/2);
		yy = reshape(xy(:,2),2,nAnchors/2);
		d = round(sqrt(diff(xx).^2 + diff(yy).^2));		% line seg lengths
		np = length(d);									% number of line segs
		data = NaN(np,5);
%		idx = NaN(1,np);								% indices of intersections
		lh = state.CLH;
		for k = 1 : np,
			x = round(linspace(xx(1,k),xx(2,k),d(k)));
			y = round(linspace(yy(1,k),yy(2,k),d(k)));
			stripe = double(img(sub2ind(size(img),y,x))) / 255;		% scale 0:1
			[n,v] = DoTrack(stripe, state.TPAR, lh(k));		% track intensity boundary
			if isempty(n),
				data(k,:) = NaN(1,5);						% no intersection
				continue;
			end;

% update marker at intersection
			if isempty(lh(k).UserData),			% make new intersection
				lh(k).UserData = line(x(n),y(n),'color',state.CONFIG.LINECOLOR, 'marker','o', ...
						'markersize',15, 'hitTest','off', 'tag','CONTOUR');
				if ~verLessThan('matlab','8.5.0'), set(lh(k).UserData,'PickableParts','none'); end;
			else,
				set(lh(k).UserData,'xdata',x(n), 'ydata',y(n));
			end;
				
			data(k,:) = [x(n),y(n),n/d(k),v,state.CURFRAME,];	% result
%			idx(k) = n;							% index of intersection
		end;
		drawnow;
		state.TRKRES = data;					% save data
		r = state;
		
		
%-----------------------------------------------------------------------------
% error

	otherwise,
		error('GCT_LINES:  unrecognized action (%s)', action);
	
end;


%=============================================================================
% DEFCFG  - set default configuration
%
%	returns default values for tPar

function tPar = DefCfg(idString)

tPar = struct('ID', idString, ...
					'FNAME', '<CMD_LINE>', ...		% export filename
					'DOFILT', 1, ...				% apply median filter
					'STRETCH', 1, ...				% stretch contrast
					'THRESH', .5, ...				% threshold for grouping (0:1)
					'ALG', 1, ...					% 1: overall max; 2: max of first thresholded group; 3: first edge above threshold
					'WTBYLAST', .5, ...				% weight by proximity to last detected intersection
					'WINLEN', 30);					% wtByLast window length
				

%=============================================================================
% DOCONFIG  - config handler
%
%   returns non-empty tPar on OK, [] on cancel

function tPar = DoConfig(idString, state, firstTime)

dfn = '<CMD_LINE>';

figPos = get(0, 'ScreenSize');
width = 290; height = 380;
figPos = [figPos(1)+(figPos(3)-width)/2, figPos(2)+(figPos(4)-height)/2, width, height];

% initialize
tPar = state.TPAR;
dPar = DefCfg(idString);
if isempty(tPar) || ~strcmp(tPar.ID, idString), 
	tPar = dPar; 
	if firstTime,
		state.TMH(end+1) = uimenu(get(state.TMH(1),'parent'),'label','Show Diagnostics','callback',{@GetContours,'TRACK','DIAGNOSTICS'});
		state.TMH(end+1) = uimenu(get(state.TMH(1),'parent'),'label','Save Tracker Data to base WS','callback',{@GetContours,'TRACK','SAVE'});
	end;
end;

cfg = dialog('Name', idString, ...
	'tag', 'GETCONTOURS', ...
	'menubar', 'none', ...
	'Position', figPos, ...
	'KeyPressFcn', 'set(gcbf,''UserData'',1);uiresume', ...
	'UserData', 0);

% about
blurb = ['This procedure assumes that paired anchors define line segments intersecting ', ...
		'the visible contour and tracks the intensity shift along them'];
if ismac, fs = 12; else; fs = 9; end;

uicontrol(cfg, ...
	'Style', 'frame', ...
	'Units', 'normalized', ...
	'Position', [0.0310 0.8132 0.9345 0.1579]);
uicontrol(cfg, ...
	'Style', 'text', ...
	'HorizontalAlignment', 'left', ...
	'String', blurb, ...
	'FontSize', fs, ...
	'Units', 'normalized', ...
	'Position', [0.0414 0.8211 0.9103 0.1421]);

% label filename
h = 17.5;
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Export Filename:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0241 0.6908 0.3862 0.0513]);
fn = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', tPar.FNAME, ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.4345 0.6908 0.4345 0.0592]);

% threshold slider and label
cbs1='set(get(gcbo,''userdata''),''string'',sprintf(''%.2f'',get(gcbo,''value'')))';
cbs2='if ~isempty(str2num(get(gcbo,''string'')))&&str2num(get(gcbo,''string''))>=0&&str2num(get(gcbo,''string''))<=1,set(get(gcbo,''userdata''),''value'',str2num(get(gcbo,''string'')));else,set(gcbo,''string'',sprintf(''%0.2f'',get(get(gcbo,''userdata''),''value'')));end';
h = h - 2.3;
if ismac, dy = 0; else, dy = .007; end;
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Threshold:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0241 0.6000 0.3138 0.0513]);
ths = uicontrol(cfg, ...
	'Style', 'slider', ...
	'Value', tPar.THRESH, ...
	'Min', 0, ...
	'Max', 1, ...
	'Callback', cbs1, ...
	'Units', 'normalized', ...
	'Position', [0.3621 0.5882+dy 0.3621 0.0592]);
the = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%0.2f',tPar.THRESH), ...
	'FontSize', 10, ...
	'Callback', cbs2, ...
	'UserData', ths, ...
	'Units', 'normalized', ...
	'Position', [0.7483 0.6000 0.1207 0.0592]);
ths.UserData = the;

% median filter checkbox and label
h = h - 2.3;
mf = uicontrol(cfg, ...
	'Style', 'checkbox', ...
	'HorizontalAlignment', 'left', ...
	'String', 'Apply Median Filter', ...
	'FontSize', 10, ...
	'Value', tPar.DOFILT, ...
	'Units', 'normalized', ...
	'Position', [0.3138 0.5092 0.4828 0.0592]);

% stretch contrast checkbox and label
h = h - 2;
str = uicontrol(cfg, ...
	'Style', 'checkbox', ...
	'HorizontalAlignment', 'left', ...
	'String', 'Stretch Contrast', ...
	'FontSize', 10, ...
	'Value', tPar.STRETCH, ...
	'Units', 'normalized', ...
	'Position', [0.3138 0.4303 0.4828 0.0592]);

% weight by last slider and label
h = h - 2.3;
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Weight by Last:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0241 0.3395 0.3138 0.0513]);
wts = uicontrol(cfg, ...
	'Style', 'slider', ...
	'Value', tPar.WTBYLAST, ...
	'Min', 0, ...
	'Max', 1, ...
	'Callback', cbs1, ...
	'Units', 'normalized', ...
	'Position', [0.3621 0.3276+dy 0.3621 0.0592]);
wte = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%0.2f',tPar.WTBYLAST), ...
	'FontSize', 10, ...
	'Callback', cbs2, ...
	'UserData', wts, ...
	'Units', 'normalized', ...
	'Position', [0.7483 0.3395 0.1207 0.0592]);
wts.UserData = wte;

% win length slider and label
cbs3='set(get(gcbo,''userdata''),''string'',sprintf(''%.0f'',get(gcbo,''value'')))';
cbs4='if ~isempty(str2num(get(gcbo,''string'')))&&str2num(get(gcbo,''string''))>=2&&str2num(get(gcbo,''string''))<=52,set(get(gcbo,''userdata''),''value'',str2num(get(gcbo,''string'')));else,set(gcbo,''string'',sprintf(''%.0f'',get(get(gcbo,''userdata''),''value'')));end';
h = h - 2.3;
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Weighting Window:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0241 0.2487 0.3138 0.0513]);
wls = uicontrol(cfg, ...
	'Style', 'slider', ...
	'Value', tPar.WINLEN, ...
	'Min', 2, ...
	'Max', 52, ...
	'Sliderstep', [.02 .2], ...
	'Callback', cbs3, ...
	'Units', 'normalized', ...
	'Position', [0.3621 0.2368 0.3621 0.0592]);
wle = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%.0f',tPar.WINLEN), ...
	'FontSize', 10, ...
	'Callback', cbs4, ...
	'UserData', wls, ...
	'Units', 'normalized', ...
	'Position', [0.7483 0.2487 0.1207 0.0592]);
wls.UserData = wle;

% algorithm popup and label
h = h - 2.3;
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Algorithm:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0241 0.1618 0.3138 0.0513]);
alg = uicontrol(cfg, ...
	'Style', 'popupmenu', ...
	'String', 'Overall Maximum|First Max above Thresh|First Edge above Thresh', ...
	'FontSize', 10, ...
	'Value', tPar.ALG, ...
	'Units', 'normalized', ...
	'Position', [0.3379 0.1579 0.5552 0.0592]);

% OK, Defaults, cancel buttons
uicontrol(cfg, ...
	'String','OK', ...
	'FontSize', 10, ...
	'Callback','set(gcbf,''UserData'',1);uiresume', ...
	'Units', 'normalized', ...
	'Position', [0.1517 0.0368 0.2069 0.0658]);
uicontrol(cfg, ...
	'String','Defaults', ...
	'FontSize', 10, ...
	'Callback','set(gcbf,''UserData'',2);uiresume', ...
	'Units', 'normalized', ...
	'Position', [0.3931 0.0368 0.2069 0.0641]);
if firstTime, es = 'off'; else, es = 'on'; end;		% cancel not permitted first time (initialization)
uicontrol(cfg, ...
	'enable', es, ...
	'String','Cancel', ...
	'FontSize', 10, ...
	'Callback','set(gcbf,''UserData'',0);uiresume', ...
	'Units', 'normalized', ...
	'Position', [0.6345 0.0368 0.2069 0.0641]);

% wait for input
while 1,
	uiwait(cfg);
	if ~ishandle(cfg),						% window closed
		tPar = []; 
		break;
	end;
	
	fName = strtok(get(fn,'string'));
	if isempty(fName) || strcmp(fName,dPar.FNAME), 
		fName = dPar.FNAME;
	else,
		[p,f,e] = fileparts(fName);
		if isempty(e), e = '.txt'; end;
		fName = fullfile(p,[f,e]);
	end;

	switch get(cfg,'userdata'),
		case 0, 		% cancel
			tPar = []; 
			break;
		case 1,			% ok
			tPar.FNAME = fName;
			tPar.THRESH = get(ths,'value');
			tPar.DOFILT = get(mf,'value');
			tPar.STRETCH = get(str,'value');
			tPar.ALG = get(alg,'value');
			tPar.WTBYLAST = get(wts,'value');
			tPar.WINLEN = round(get(wls,'value'));
			if mod(tPar.WINLEN,2), tPar.WINLEN = tPar.WINLEN-1; end;
			break;
		case 2,			% defaults
			set(fn,'string', dPar.FNAME);
			set(ths,'value', dPar.THRESH);
			set(the,'string', sprintf('%0.2f',dPar.THRESH));
			set(mf,'value', dPar.DOFILT);
			set(str,'value', dPar.STRETCH);
			set(wts,'value', dPar.WTBYLAST);
			set(wte,'string', sprintf('%0.2f',dPar.WTBYLAST));
			set(wls,'value', dPar.WINLEN);
			set(wle,'string', sprintf('%.0f',dPar.WINLEN));
			set(alg,'value', dPar.ALG);
			continue;
	end;
end;
if ishandle(cfg), delete(cfg); end;


%=============================================================================
% DOTRACK  - detect intensity boundary
%
%	lh is handle of current anchor pair line segment
%
%   returns index of point along stripe, [] on fail and associated intensity

function [n,v] = DoTrack(stripe, tPar, lh)

% filter if necessary
if tPar.DOFILT, stripe = medfilt1(stripe); end;		% default 3-order
urStripe = stripe;
if tPar.STRETCH, stripe = (stripe - min(stripe)) ./ range(stripe); end;

% weight by last intersection 
if tPar.WTBYLAST > 0,
	if ishandle(lh) & ishandle(lh.UserData),
		x = [lh.XData(1) lh.UserData.XData];
		y = [lh.YData(1) lh.UserData.YData];
		idx = round(sqrt(diff(x).^2+diff(y).^2));	% index of intersection along line segment
		if idx < 1, idx = 1; elseif idx>length(stripe), idx = length(stripe); end;
		if ~isnan(idx),
			winLen = tPar.WINLEN;				% weighted neighborhood
			winBump = tPar.WTBYLAST;			% weight magnitude
			win = hamming(winLen*2+1);
			win = win - min(win);
			win = win ./ max(win) * winBump + 1-winBump;			% window 1 at idx, 1-winBump at tails
			wt = (1-winBump) * ones(1,length(stripe)+2*winLen);
			ht = [-1 1]*winLen + idx + winLen;
			wt(ht(1):ht(2)) = win;
			wt([1:winLen end-winLen+1:end]) = [];
			stripe = stripe .* wt;
		end;
	end;
end;

% apply detection algorithm:  find index of intersection N associated with max intensity V
switch tPar.ALG,

	case 1, 		% overall maximum
		[~,n] = max(stripe);
		
	case 2,			% max of first group exceeding threshold
		stripe(stripe < tPar.THRESH) = 0;			% clip sub-threshold values
		h = find(stripe > 0, 1);					% start of first group exceeding threshold
		if isempty(h),								% nothing above threshold
			n = []; v = [];
			return;
		end;
		t = find(stripe(h:end) == 0, 1) + h - 2;	% end of first group
		[~,n] = max(stripe(h:t));
		n = n + h - 1;

	case 3,			% first edge above threshold
		stripe(stripe < tPar.THRESH) = 0;			% clip sub-threshold values
		n = find(stripe > 0, 1);
		v = stripe(n);
		
end;

v = urStripe(n);		% unstretched intensity


%=============================================================================
% HANDLECLICK  - process click on anchor pair line
%
%	used to update preferred intersection point
%
% passed line handle and event properties

function HandleClick(lh, evt)

xy = round(evt.IntersectionPoint);
state = get(gcf, 'userdata');

% update graphics
if isempty(lh.UserData),		% make new intersection
	lh.UserData = line(xy(1),xy(2),'color',state.CONFIG.LINECOLOR, 'marker','o', ...
				'markersize',15, 'hitTest','off', 'tag','CONTOUR');
	if ~verLessThan('matlab','8.5.0'), set(lh.UserData,'PickableParts','none'); end;
else,							% update existing intersection
	set(lh.UserData,'xdata',xy(1), 'ydata',xy(2));
end;


%=============================================================================
% PLOTSTRIPE  - display stripe threshold

function PlotStripe(state)

switch size(state.ANCHORS,1),
	case {0,1}, 
		fprintf('no anchor pairs defined\n');
		return;
	case {2,3}
		xy = state.ANCHORS(1:2,:);
	otherwise,
		fprintf('last defined anchor pair\n');
		xy = state.ANCHORS;
		if mod(size(xy,1),2), xy(end,:) = []; end;
		xy = xy(end-1:end,:);
end;
tPar = state.TPAR;
img = get(state.IH,'cdata');
d = round(sqrt(sum(diff(xy).^2,2)));
x = round(linspace(xy(1,1),xy(2,1),d));
y = round(linspace(xy(1,2),xy(2,2),d));
stripe = double(img(sub2ind(size(img),y,x))) / 255;
if tPar.DOFILT, stripe = medfilt1(stripe); end;
maxInt = max(stripe);

% plot stripe
n = state.TRKRES(end,3) * length(stripe);	% intersection point
figure;
stem(stripe);
xl = [0 length(x)+1];
set(gca,'xlim',xl);
h = line(xl',[1;1]*tPar.THRESH,'color','r','linestyle','--');
h(2) = line([n,n],get(gca,'ylim'),'color',[0 .6 0],'linewidth',2);
legend(h, 'Current threshold','Last intersection');
ylabel('intensity');
xlabel('pixels along stripe')
title('Image intensity along stripe defined by last defined anchor pair');
uicontrol('style', 'pushbutton', ...
			'units', 'characters', ...
			'string', 'Cursor', ...
			'FontSize', 10, ...
			'position', [2 1 10 1.5], ...
			'callback', 'ginput(1)');


