function r = gct_Tahltan(action, state, varargin)
%GCT_TAHLTAN  - GetContours tracker specific to the DHW Tahltan project
%
%	usage:  r = gct_Tahltan(action, state, ...)
%
% The Tahltan project collects data in X-plane mode, with midsagittal tongue left of image midline, 
% and a coronal cross-section right of image midline.  This procedure assumes default behavior for  
% the midsagittal anchors left of image midline, and uses the first three defined anchors on the
% image right of midline to define the left/right parasagittal tongue peaks, and the midsagittal
% groove trough (anchors are sorted based on horizontal position, which is assumed to be subject
% left/right.

% mkt 09/20 for UltraFest IX

%	branch by ACTION

switch upper(action),
		
%-----------------------------------------------------------------------------
% ADDPT:  add new point flag (called by GetContours:DOWN handler)
%
% returns 
%   0 indicating new point should be appended using default processing (sagittal side)
%   1 indicating new point should be appended to end (coronal side)
%  -1 indicating new point should be ignored (corornal side with three pts already)

	case 'ADDPT',
		if isempty(state.ANCHORS),					% normal processing
			r = 0;
			return;
		end;
		cp = varargin{1};							% new point
		
% split into left (sagittal) and right (coronal) anchors
		[h,w] = size(get(state.IH,'CDATA'));
		split = round(w/2) + 20;					% account for greater left margin (approximately)
		if cp(1) < split, r = 0; return; end;		% sagittal side, normal processing
		corA = find(state.ANCHORS(:,1) >= split);	% coronal anchors
		if length(corA) < 3, r = 1; return; end;	% coronal side, append point
		r = -1;										% coronal side, already have 3 pts – ignore
		
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


%-----------------------------------------------------------------------------
% EXPORT:  output tracker data to file in tab-delimited format

	case 'EXPORT',
		fprintf('GCT_TAHLTAN:  export not available\n');
		
		
%-----------------------------------------------------------------------------
% PLOT:  post-tracking plotting handler (called after CLH creation)
%
% overrides contour plotting to distinguish between sagittal (left) and coronal (right) anchors
%
%	returns updated STATE variable 

	case 'PLOT',
		nAnchors = size(state.ANCHORS,1);

% split into left (sagittal) and right (coronal) anchors
		anchors = state.ANCHORS;
		[h,w] = size(get(state.IH,'CDATA'));
		split = round(w/2) + 20;					% account for greater left margin (approximately)
		if nAnchors > 0,
			sagA = find(anchors(:,1) < split);				% sagittal anchors
			corA = find(anchors(:,1) >= split);				% coronal anchors
		else,
			sagA = [];
			corA = [];
		end;
		
% update sagittal contour line from anchor points
		switch length(sagA),
			case 0,
				;
			case 1,
				state.XY = ones(state.NPOINTS,1) * anchors(sagA,:);
				set(state.CLH,'xdata',anchors(sagA,1),'ydata',anchors(sagA,2));
			case 2,
				k = [0 ; sqrt(sum(diff(anchors(sagA)).^2,2))];
				state.XY = interp1(k,anchors(sagA,:),linspace(0,k(end),state.NPOINTS),'linear');
				set(state.CLH,'xdata',state.XY(:,1),'ydata',state.XY(:,2));
			otherwise,
				k = [0 ; cumsum(sqrt(sum(diff(anchors(sagA)).^2,2)))];
				state.XY = interp1(k,anchors(sagA,:),linspace(0,k(end),state.NPOINTS),'pchip');
				set(state.CLH,'xdata',state.XY(:,1),'ydata',state.XY(:,2));
		end;

% delete any existing coronal lines, recreate handles
		delete([state.TPAR.PLH , state.TPAR.MLH]);
		if length(corA) > 0,
			PLH = line(anchors(corA(1),1),anchors(corA(1),2),'color',state.CONFIG.LINECOLOR,'linewidth',state.CONFIG.LINEWIDTH,'tag','CONTOUR','hitTest','off');
			MLH = line(anchors(corA(1),1),anchors(corA(1),2),'color',state.CONFIG.LINECOLOR,'linewidth',state.CONFIG.LINEWIDTH,'tag','CONTOUR','hitTest','off');
			if ~verLessThan('matlab','8.5.0'), set([PLH,MLH],'PickableParts','none'); end;
			uistack([PLH,MLH],'bottom'); uistack([PLH,MLH],'up');
			state.TPAR.PLH = PLH;
			state.TPAR.MLH = MLH;
		end;

% update coronal lines from anchor points
		switch length(corA),
			case 0,
				state.TPAR.PLH = []; 
				state.TPAR.MLH = []; 
			case {1,2},
				set([PLH,MLH],'xdata',anchors(corA,1),'ydata',anchors(corA,2));
			otherwise,
				axy = anchors(corA,:);
				axy = axy(1:3,:);				% ignore additional points
				[~,k] = sort(axy(:,1));			% sort left to right
				axy = axy(k,:);
				ixy = FindNormal(axy);			% intersection point of line connecting parasagittal peaks to midsagittal trough
				set(PLH,'xdata',axy([1 3],1),'ydata',axy([1 3],2));
				set(MLH,'xdata',[ixy(1);axy(2,1)],'ydata',[ixy(2);axy(2,2)]);

% always evaluate
				state.TRKRES = sqrt(sum((ixy - axy(2,:)).^2));
				if ixy(2) > axy(2,2), state.TRKRES = -state.TRKRES; end;	% negative if point above line
		end;

		r = state;		
		
%-----------------------------------------------------------------------------
% SAVE:  save tracker data as variable in base ws
%
%	returns nothing
%
% this handler saves TRKRES data as VNAME_trk, formatted as [nFrames x F,T,GLEN]
% where F is frame number, T is offset in secs, GLEN is groove length value

	case 'SAVE',
		v = evalin('base',state.VNAME);			% frame data in base ws (VNAME)
		f = cell2mat({v.FRAME});				% frames with data
		k = (state.CURFRAME == f);
		v(k).XY = state.XY;						% update output variabe
		v(k).ANCHORS = state.ANCHORS;
		v(k).TRKRES = state.TRKRES;
		[~,k] = sort(f);
		v = v(k);								% sorted	
		v(find(cellfun(@isempty,{v.TRKRES}))) = [];	% drop frame w/o TRKRES
		if isempty(v), 
			fprintf('nothing to save\n');
			return;
		end;
		d = zeros(length(v),3);
		d(:,1) = cell2mat({v.FRAME})';
		d(:,2) = (d(:,1)-1)/state.FRATE;
		d(:,3) = cell2mat({v.TRKRES})';
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
% state.TRKRES set to length of normal connecting midsagittal point to line connecting parasagittal peaks
% if that point is above the parasagittal line the distance is returned as negative

	case 'TRACK',
		PLH = state.TPAR.PLH;
		MLH = state.TPAR.MLH;
		if isempty(PLH), r = []; return; end;		% no coronal points
		pxy = [PLH.XData ; PLH.YData]';
		mxy = [MLH.XData ; MLH.YData]';
		if max(abs(pxy-mxy),[],'all') == 0, r = []; return; end;	% no midsagittal point
		state.TRKRES = sqrt(sum(diff(mxy).^2,2));	% distance along normal
		if mxy(2,2) < mxy(1,2), state.TRKRES = -state.TRKRES; end;	% negative if point above line
		r = state;
		
		
%-----------------------------------------------------------------------------
% error

	otherwise,
		error('GCT_TAHLTAN:  unrecognized action (%s)', action);
	
end;


%=============================================================================
% DEFCFG  - set default configuration
%
%	returns default values for tPar

function tPar = DefCfg(idString)

tPar = struct('ID', idString, ...
				'PLH', [], ...			% line handle connecting parasagittal peaks
				'MLH', []);				% line handle connecting midsagittal groove to its intersection with PLH
				

%=============================================================================
% DOCONFIG  - config handler
%
%   returns tPar

function tPar = DoConfig(idString, state, firstTime)

tPar = state.TPAR;
dPar = DefCfg(idString);

figPos = get(0, 'ScreenSize');
width = 300; height = 240;
figPos = [figPos(1)+(figPos(3)-width)/2, figPos(2)+(figPos(4)-height)/2, width, height];
if isempty(tPar) || ~strcmp(tPar.ID, idString), 
	tPar = dPar; 
	if firstTime,
		state.TMH(end+1) = uimenu(get(state.TMH(1),'parent'),'label','Save Tracker Data to base WS','callback',{@GetContours,'TRACK','SAVE'});
	end;
end;

cfg = dialog('Name', idString, ...
	'tag', 'GETCONTOURS', ...
	'menubar', 'none', ...
	'Position', figPos, ...
	'KeyPressFcn', 'uiresume');
	
% about
blurb = ['This procedure is intended for the DHW Tahltan project. Images are assumed to represent ', ...
		'X-plane mode representations with midsagittal views left of midline and coronal views ', ...
		'right of midline. Anchors left of image midline use default processing. The first three ', ...
		'anchors right of midline are assumed to designate parasagittal left/right tongue peaks ', ...
		'and the midsagittal groove (distinguished by left/right position).'];

uicontrol(cfg, ...
	'Style', 'frame', ...
	'Position', [10 height-180 width-19 160]);
uicontrol(cfg, ...
	'Style', 'text', ...
	'HorizontalAlignment', 'left', ...
	'String', blurb, ...
	'FontSize', 12, ...
	'Position', [14 height-177 width-27 154]);

% OK buttons
uicontrol(cfg, ...
	'Position',[width/2-30 15 60 25], ...
	'String','OK', ...
	'FontSize', 10, ...
	'Callback','uiresume');

% wait for input
uiwait(cfg);
if ishandle(cfg), delete(cfg); end;


%=============================================================================
% FINDNORMAL  - find point on AB normal to point C
%
% assumes xy is ordered [A (left), C (mid), B (right)]; finds P intersecting A:B on its normal to C

function P = FindNormal(xy)

X = 1; Y = 2;
A = xy(1,:);
B = xy(3,:);
C = xy(2,:);
r = ((A(Y)-C(Y))*(A(Y)-B(Y)) - (A(X)-C(X))*(B(X)-A(X))) / ((B(X)-A(X))^2 + (B(Y)-A(Y))^2);
P = A + r*(B-A);
