function r = gct_template(action, state, varargin)
%GCT_TEMPLATE  - GetContours tracking procedure
%
%	usage:  r = gct_template(action, state, ...)
%
% This is a template GetContours tracking procedure.  
%
% User labelling procedures must supply these action handlers:
%	ADDPT	- add new point flag
%	CONFIG	- process any configuration necessary
%	DIAGS   - show tracker-specific diagnostics
%	EXPORT  - export data to tab-delimited file
%	PLOT	- post-tracking plotting handler (called by GetContours:UpdateContour)
%	SAVE	- save tracker data as base ws variable
%	TRACK	- track current frame handler
%
% if the TRACK handler returns state the value of state.TRKRES it is attached to 
% the frame values (here it is set to the s.d. of the jitter across all anchors)

% mkt 06/18
% mkt 12/19 v2.3

%	branch by ACTION

switch upper(action),
		
%-----------------------------------------------------------------------------
% ADDPT:  add new point flag (called by GetContours:DOWN handler)
%
%	returns 1 indicating new point should always be appended to anchor list
%	    or  0 indicating new point should be appended using default processing
%       or -1 indicating new point should be ignored

	case 'ADDPT',
		r = 0;
		
		
%-----------------------------------------------------------------------------
% CONFIG:  handle configuration
%
% 	returns TPAR = new tracker params
%
%	varargin{1} set to 1 on first entry; 0 on subsequent entries
%
% for a minimal handler just set r = 1 and return
% r = [] flags user cancelled

	case 'CONFIG',
		r = DoConfig(mfilename, state, varargin{1});


%-----------------------------------------------------------------------------
% DIAGS:  show diagnostics
%
% 	return ignored

	case 'DIAGS',
		fprintf('diagnostics for %s\n', state.TPAR.ID);


%-----------------------------------------------------------------------------
% EXPORT:  output tracker data to file in tab-delimited format
%
%	returns nothing; minimal handler does nothing
%
% this handler returns a count of frames with non-empty TRKRES

	case 'EXPORT',
		v = evalin('base',state.VNAME);			% frame data in base ws (VNAME)
		f = cell2mat({v.FRAME});				% frames with tracking data
		n = length(f) - sum(cellfun(@isempty, {v.TRKRES}));
		fprintf('%d frames have tracking data\n', n);
		
		
%-----------------------------------------------------------------------------
% PLOT:  post-tracking plotting handler (called after CLH creation)
%
%	returns updated STATE variable if plotting handled locally, 
%	[] to use GetContours default plotting

	case 'PLOT',
		r = [];
		
		
%-----------------------------------------------------------------------------
% SAVE:  save tracker data as variable in base ws
%
%	returns nothing
%
% this handler saves TRKRES data as VNAME_trk

	case 'SAVE',
		v = evalin('base',state.VNAME);			% frame data in base ws (VNAME)
		f = cell2mat({v.FRAME});				% frames with data
		[f,k] = sort(f);
		v = v(k);								% sorted	
		d = [f ; cell2mat({v.TRKRES})]';
		vName = sprintf('%s_trk', state.VNAME);
		assignin('base', vName, d);
		fprintf('%s created in base workspace\n', vName);
		
		
%-----------------------------------------------------------------------------
% TRACK:  track current frame handler
%
%	returns updated STATE variable on success, [] on failure
%
% 	varargin{1} is nonzero on first track of a sequence
%
% for this example all anchors are jittered randomly with Gaussian noise
% and state.TRKRES is set to the s.d. of the jitter

	case 'TRACK',
		if nargin < 3 || isempty(varargin{1}), firstTime = 1; else, firstTime = 0; end;
		nAnchors = size(state.ANCHORS,1);
		if nAnchors < 1,
			r = [];
			return;
		end;

% nothing is done with the image here, but this shows how to obtain it based on current state parameters
		img = GetContours('GETMOVIEFRAME',state.CURFRAME);
		
		orig = state.ANCHORS;
		state.ANCHORS = state.ANCHORS + randn(size(state.ANCHORS));
		d = sqrt(sum((orig - state.ANCHORS).^2,2));
		state.TRES = std(d);
		for k = 1 : nAnchors,
			set(state.ALH(k),'xdata',state.ANCHORS(k,1),'ydata',state.ANCHORS(k,2));
		end;

% update contour line from anchor points
		switch size(state.ANCHORS,1),
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
		state.TRKRES = std(sqrt(sum((orig - state.ANCHORS).^2,2)));
		r = state;
		
		
%-----------------------------------------------------------------------------
% error

	otherwise,
		error('GCT_TEMPLATE:  unrecognized action (%s)', action);
	
end;


%=============================================================================
% DEFCFG  - set default configuration
%
%	returns default values for tPar

function tPar = DefCfg(idString)

tPar = struct('ID', idString, ...
					'VALUE', 0);
				

%=============================================================================
% DOCONFIG  - config handler
%
%   returns non-empty tPar on OK, [] on cancel

function tPar = DoConfig(idString, state, firstTime)

figPos = get(0, 'ScreenSize');
width = 300; height = 170;
figPos = [figPos(1)+(figPos(3)-width)/2, figPos(2)+(figPos(4)-height)/2, width, height];

% initialize if necessary
tPar = state.TPAR;
if isempty(tPar) || ~strcmp(tPar.ID, idString),
	tPar = DefCfg(idString);
end;

cfg = dialog('Name', idString, ...
	'tag', 'GETCONTOURS', ...
	'menubar', 'none', ...
	'Position', figPos, ...
	'KeyPressFcn', 'set(gcbf,''UserData'',1);uiresume', ...
	'UserData', 0);

% about
blurb = ['This is a sample GetContours tracking procedure.  When called ', ...
		'for tracking, it jitters all anchors randomly with Gaussian noise.'];

uicontrol(cfg, ...
	'Style', 'frame', ...
	'Position', [10 height-70 width-20 60]);
uicontrol(cfg, ...
	'Style', 'text', ...
	'Fontsize', 12, ...
	'HorizontalAlignment', 'left', ...
	'String', blurb, ...
	'Position', [13 height-67 width-26 54]);

% editable text field and label
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Parameter:', ...
	'Units', 'characters', ...
	'Position', [1 4 21 1.5]);
eh = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', num2str(tPar.VALUE), ...
	'Units', 'characters', ...
	'Position', [23 4.2 8 1.5], ...
	'Callback', 'set(gcbf,''UserData'',1);uiresume');

% OK, cancel buttons
uicontrol(cfg, ...		% buttons
	'Position',[width/2-70 15 60 25], ...
	'String','OK', ...
	'Callback','set(gcbf,''UserData'',1);uiresume');
if firstTime, es = 'off'; else, es = 'on'; end;		% cancel not permitted first time (initialization)
uicontrol(cfg, ...
	'enable', es, ...
	'Position',[width/2+10 15 60 25], ...
	'String','Cancel', ...
	'Callback','uiresume');

% wait for input
uiwait(cfg);
if ~ishandle(cfg), 
	tPar = []; 
	return;
end;
if get(cfg, 'UserData'),
	v = str2num(get(eh, 'string'));
	if ~isempty(v),
		tPar.VALUE = v;
	end;
else,
	tPar = [];
end;
if ishandle(cfg), delete(cfg); end;

