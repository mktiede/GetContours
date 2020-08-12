function r = gct_snake(action, state, varargin)
%GCT_SNAKE  - GetContours tracker that applies the default snake algorithm to successive frames
%
%	usage:  r = gct_snake(action, state, ...)
%
% This procedure uses the default snake algorithm (the EdgeTrak implementation of C. Laporte) with
% the relevant parameters configured using the main GetContours configuration dialog.  After 
% seeding of the first frame with a minimum of three points, subsequent frames are seeded using the 
% fit from the preceding frame.

% mkt 08/20 UltraFest IX release

%	branch by ACTION

switch upper(action),
		
%-----------------------------------------------------------------------------
% ADDPT:  add new point flag (called by GetContours:DOWN handler)
%
%	returns 0 indicating new point should be appended using default processing

	case 'ADDPT',
		r = 0;
		
		
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
		v = evalin('base',state.VNAME);		% load output variable
		frames = cell2mat({v.FRAME});
		[~,k] = sort(frames);
		v = v(k);
		v(cellfun(@isempty,{v.XY})) = [];	% delete empty frames
		frames = cell2mat({v.FRAME});
		if length(frames) < 2, return; end;
		ht = [frames(1) frames(end)];
		nf = diff(ht) + 1;
		energy = NaN(state.PARAMS.NSPTS,nf);
		for k = 1:length(frames),
			if isempty(v(k).TRKRES), continue; end;
			energy(:,frames(k)-ht(1)+1) = v(k).TRKRES;
		end;
		energy(end,:) = [];						% kill NaNs associated with last point
		e1 = energy(:,1) * ones(1,size(energy,2));
		energy = (energy - e1) ./ e1;			% normalize by energy of first frame
		energy = energy .* (energy>0) + 1;		% reserve 0 for no data
		energy(isnan(energy)) = 0;
		pos = get(0,'defaultFigurePosition');
		pos(1) = 5;
		
		fh = findobj('TAG','GCT_SNAKE');
		if isempty(fh),
			fh = figure('name','SLURP Diagnostics','position',pos,'tag','GCT_SNAKE','numbertitle','off','menubar','none');
		else,
			clf(fh);
			figure(fh);
		end;
		set(fh,'userData',1);	% used to avoid frame switch on double call to update fcn
		yy = linspace(1,size(state.ANCHORS,1),state.PARAMS.NSPTS);
		imagesc(flipud(energy),[0 4]);
		x = str2num(char(get(gca,'xticklabel')));
		x = cellstr(num2str(x + ht(1) - 1));
		y = str2num(char(get(gca,'yticklabel')));
		yy = flipud(linspace(1,size(state.ANCHORS,1),state.PARAMS.NSPTS)');
		y = cellstr(num2str(round(yy(y))));
		
		c = [[1 1 1];eval(sprintf('%s(63)',state.TPAR.CMAP))];
		colormap(c);
		xlabel('frames'); ylabel('left <- points -> right'); title('Energy')
		set(gca,'userdata',{ht(1),yy},'xticklabel',x,'yticklabel',y);

		h = datacursormode(fh);
		set(h,'enable','on','updateFcn',@SetFrame);


%-----------------------------------------------------------------------------
% EXPORT:  output tracker data to file in tab-delimited format (use default GC behavior)

	case 'EXPORT',
		GetContours('EXPORT');	
		
		
%-----------------------------------------------------------------------------
% PLOT:  post-tracking plotting handler (called after CLH creation)
%
%	returns [] to use GetContours default plotting

	case 'PLOT',
		r = [];
		
		
%-----------------------------------------------------------------------------
% SAVE:  save tracker data as variable in base ws (UNUSED)
%


%-----------------------------------------------------------------------------
% TRACK:  track current frame handler
%
%	returns updated STATE variable on success, [] on failure
%
% 	varargin{1} is nonzero on first track of a sequence
%
% TRKRES returns energy for frame fit

	case 'TRACK',
		firstTime = varargin{1};
		nAnchors = size(state.ANCHORS,1);
		if nAnchors < 3,
			fprintf('SNAKE requires at least three anchor points for initialization\n');
			r = [];
			return;
		end;
		img = im2double(GetContours('GETMOVIEFRAME',state.CURFRAME));
		xy = state.XY;
		k = [0 ; cumsum(sqrt(sum(diff(xy).^2,2)))];
		try,
			eg = ComputeImageForces(img,state.PARAMS.SIGMA);
			useBE = 1;		% use band energy
			xy = interp1(k, xy, linspace(0, k(end), state.PARAMS.NSPTS)', 'pchip');
			[xy,state.TRKRES] = make_snake(img', eg', xy, state.PARAMS.DELTA*ones(state.PARAMS.NSPTS,1), state.PARAMS.BPEN, state.PARAMS.ALPHA, state.PARAMS.LAMBDA, useBE);
		catch,
			fprintf('failed to apply snake\n');
			return;
		end;
		k = [0 ; cumsum(sqrt(sum(diff(xy).^2,2)))];		
		state.XY = interp1(k, xy, linspace(0, k(end), state.NPOINTS)', 'pchip');
		set(state.CLH,'xdata',state.XY(:,1),'ydata',state.XY(:,2));
		state.ANCHORS = interp1(k, xy, linspace(0, k(end), nAnchors)', 'pchip');
		for k = 1 : nAnchors,				% redistribute anchors
			set(state.ALH(k),'xdata',state.ANCHORS(k,1),'ydata',state.ANCHORS(k,2));
		end;
		drawnow;
		r = state;
		
		
%-----------------------------------------------------------------------------
% error

	otherwise,
		error('GCT_SNAKE:  unrecognized action (%s)', action);
	
end;


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
% DEFCFG  - set default configuration
%
%	returns default values for tPar

function dPar = DefCfg(idString)

dPar = struct('ID', idString, ...
				'CMAP', 'gray');					% diagnostic map colormap name


%=============================================================================
% DOCONFIG  - config handler
%
%   returns non-empty tPar on OK, [] on cancel

function tPar = DoConfig(idString, state, firstTime)

tPar = state.TPAR;

% initialize
dPar = DefCfg(idString);
if isempty(tPar) || ~strcmp(tPar.ID, idString),
	tPar = dPar;
	if firstTime,
		state.TMH(end+1) = uimenu(get(state.TMH(1),'parent'),'label','Show Diagnostics','callback',{@GetContours,'TRACK','DIAGNOSTICS'});
	end;
end;

figPos = get(0, 'ScreenSize');
width = 290; height = 230;
figPos = [figPos(1)+(figPos(3)-width)/2, figPos(2)+(figPos(4)-height)/2, width, height];

cfg = dialog('Name', idString, ...
	'tag', 'GETCONTOURS', ...
	'menubar', 'none', ...
	'Position', figPos, ...
	'KeyPressFcn', 'set(gcbf,''UserData'',1);uiresume', ...
	'UserData', 0);

% about
blurb = ['This procedure uses the default GC snake algorithm (the EdgeTrak implementation of C. Laporte) ', ...
		'with the relevant parameters configured using the main GetContours configuration dialog.  After ', ...
		'seeding of the first frame with a minimum of three points, subsequent frames are seeded using ', ...
		'the fit from the preceding frame.'];
if ispc, fs = 8; else; fs = 12; end;

uicontrol(cfg, ...
	'Style', 'frame', ...
	'Units', 'normalized', ...
	'Position', [0.0310 0.4304 0.9345 0.5217]);
uicontrol(cfg, ...
	'Style', 'text', ...
	'HorizontalAlignment', 'left', ...
	'String', blurb, ...
	'FontSize', fs, ...
	'Units', 'normalized', ...
	'Position', [0.0448 0.4435 0.9069 0.4957]);

% colormap
h = 4;
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Colormap:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0241 0.2609 0.3862 0.0848]);
cm = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', tPar.CMAP, ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.4345 0.2609 0.4828 0.0978]);

% OK, Defaults, cancel buttons
uicontrol(cfg, ...
	'String','OK', ...
	'FontSize', 10, ...
	'Callback','set(gcbf,''UserData'',1);uiresume', ...
	'Units', 'normalized', ...
	'Position', [0.1517 0.0609 0.2069 0.1087]);
uicontrol(cfg, ...
	'String','Defaults', ...
	'FontSize', 10, ...
	'Callback','set(gcbf,''UserData'',2);uiresume', ...
	'Units', 'normalized', ...
	'Position', [0.3931 0.0609 0.2069 0.1087]);
if firstTime, es = 'off'; else, es = 'on'; end;		% cancel not permitted first time (initialization)
uicontrol(cfg, ...
	'enable', es, ...
	'String','Cancel', ...
	'FontSize', 10, ...
	'Callback','set(gcbf,''UserData'',0);uiresume', ...
	'Units', 'normalized', ...
	'Position', [0.6345 0.0609 0.2069 0.1087]);

% wait for input
while 1,
	uiwait(cfg);
	if ~ishandle(cfg),						% window closed
		tPar = []; 
		break;
	end;
	
	cName = strtok(get(cm,'string'));
	try,
		eval(sprintf('%s(63);',cName));
	catch,
		cName = [];
	end;
	if isempty(cName), 
		cName = dPar.CMAP;
	end;

	switch get(cfg,'userdata'),
		case 0, 		% cancel
			tPar = []; 
			break;
		case 1,			% ok
			tPar.CMAP = cName;
			break;
		case 2,			% defaults
			set(cm,'string',dPar.CMAP);
			continue;
	end;
end;
if ishandle(cfg), delete(cfg); end;


%=============================================================================
% SETFRAME  - DataCursorManager UpdateFcn that sets current frame based on
%				clicked position in diagnostic map

function txt = SetFrame(h, evt)

pos = round(evt.Position);		% X,Y
ud = get(gca,'userdata');
fh = get(gca,'parent');
doFlag = get(fh,'userdata');
set(fh,'userData',~doFlag);
ih = get(gca,'children');
img = ih.CData;
energy = img(pos(2),pos(1)) - 1;	% report true energy value
frame = pos(1) + ud{1} - 1;			% add frame offset
yy = ud{2};
pt = round(yy(pos(2)));				% interpolate to # anchors
txt = {['Frame: ',num2str(frame)],['Point: ',num2str(pt),'  Energy: ',sprintf('%.2f',energy)]};
if doFlag, GetContours('FRAME','EXPLICIT',frame); end;
