function r = gct_SLURP(action, state, varargin)
%GCT_SLURP  - GetContours tracker implementing SLURP algorithm
%
%	usage:  r = gct_SLURP(action, state, ...)
%
% This procedure implements the SLURP contour tracking algorithm of Cathy Laporte.
% It assumes the availability of the following supporting files within the same
% directory:
%   make_snake (mex)          - CL implementation of Li et al. (2005) algorithm
%   SLURP_Default_Model (mat) - default models for snake particle initialization
%
% Publications which make use of this software should cite Laporte C & Ménard L. (2018). 
% Multi-hypothesis tracking of the tongue surface in ultrasound video recordings of normal 
% and impaired speech.  Medical Image Analysis, 44, 98-114.

% C. Laporte 11/19
% mkt 01/20 various tweaks
% mkt 08/20 UltraFest IX release
% mkt 01/21 fix max init issue

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
		nParticles = NaN(1,nf);
		for k = 1:length(frames),
			if isempty(v(k).TRKRES) || ~isfield(v(k).TRKRES,'ENERGY'), continue; end;
			energy(:,frames(k)-ht(1)+1) = v(k).TRKRES.ENERGY;
			nParticles(frames(k)-ht(1)+1) = v(k).TRKRES.NPARTICLES;
		end;
		energy(end,:) = [];						% kill NaNs associated with last point
		energy = energy .* (energy>0) + 1;		% reserve 0 for no data
		energy(isnan(energy)) = 0;
		pos = get(0,'defaultFigurePosition');
		pos(1) = 5;
		
		fh = findobj('TAG','GCT_SLURP');
		if isempty(fh),
			fh = figure('name','SLURP Diagnostics','position',pos,'tag','GCT_SLURP','numbertitle','off','menubar','none');
		else,
			clf(fh);
			figure(fh);
		end;
		set(fh,'userData',1);	% used to avoid frame switch on double call to update fcn
		subplot(211);
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
		set(gca,'userdata',{1,ht(1),yy},'xticklabel',x,'yticklabel',y);

		subplot(212);
		stem(ht(1):ht(2),nParticles);
		set(gca,'xlim',ht,'ylim',[1 500],'userdata',{2,0,0});
		xlabel('frames'); ylabel('# particles'); title('# of Particles')

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
% TRKRES returns normalized energy and particle count for frame fit

	case 'TRACK',
		firstTime = varargin{1};
		nAnchors = size(state.ANCHORS,1);
		if nAnchors < 3,
			fprintf('SLURP requires at least three anchor points for initialization\n');
			r = [];
			return;
		end;
		img = im2double(GetContours('GETMOVIEFRAME',state.CURFRAME));
		if firstTime,
			k = [0 ; cumsum(sqrt(sum(diff(state.XY).^2,2)))];
			init_pts = interp1(k,state.XY,linspace(0,k(end),state.TPAR.NPOINTS),'pchip');
		else,
			init_pts = [];
		end;
		[state.TPAR,xy,energy,nParticles] = DoTrack(img, state.TPAR, init_pts, firstTime);
		k = [0 ; cumsum(sqrt(sum(diff(xy).^2,2)))];		
		state.XY = interp1(k, xy, linspace(0, k(end), state.NPOINTS)', 'pchip');
		set(state.CLH,'xdata',state.XY(:,1),'ydata',state.XY(:,2));
		state.ANCHORS = interp1(k, xy, linspace(0, k(end), nAnchors)', 'pchip');
		for k = 1 : nAnchors,				% redistribute anchors
			set(state.ALH(k),'xdata',state.ANCHORS(k,1),'ydata',state.ANCHORS(k,2));
		end;
		drawnow;
		normEnergy = (energy - state.TPAR.STARTENERGY) ./ state.TPAR.STARTENERGY;
		state.TRKRES = struct('ENERGY',normEnergy,'NPARTICLES',nParticles);
		r = state;
		

%-----------------------------------------------------------------------------
% error

	otherwise,
		error('GCT_SLURP:  unrecognized action (%s)', action);
	
end;


%=============================================================================
% COMPUTEIMAGEFORCES  - computes Gaussian-filtered image derivatives
%
% N.B. unlike GC version does not return 1-result to delay normalization until after masking

function img = ComputeImageForces(img, sigma)

[x,y] = ndgrid(floor(-3*sigma):ceil(3*sigma),floor(-3*sigma):ceil(3*sigma));
g = exp(-(x.^2 + y.^2)/(2*sigma^2));
d = 2*pi*sigma^4;
ix = imfilter(img, -(x./d).*g);
iy = imfilter(img, -(y./d).*g);
img = sqrt(ix.*ix + iy.*iy);


%=============================================================================
% DEFCFG  - set default configuration
%
%	returns default values for tPar

function dPar = DefCfg(idString)

% load default models
defaultModelFname = 'SLURP_Default_Model';	
nPoints = 39;

dPar = struct('ID', idString, ...
				'MNAME', defaultModelFname, ...		% default models mat file
				'CMAP', 'gray', ...					% diagnostic map colormap name
				'SIGMA', 5, ...						% image forces Gaussian sigma
				'DELTA', 2, ...						% delta
				'BPEN', 2, ...						% band penalty
				'ALPHA', 0.8, ...					% alpha
				'LAMBDA', 0.95, ...					% lambda1
				'MINPART', 10, ...					% min # particles
				'MAXPART', 1000, ...				% max # particles
				'ADAPTIVE', 1, ...					% adaptive sampling 
				'NPOINTS', nPoints);


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
	try,
		dm = load(fullfile(fileparts(which(idString)),dPar.MNAME));
		nPoints = round(length(dm.ShapeData.x_mean)/2);	% 39
	catch,
		fprintf('missing required %s\n', dPar.MNAME);
		tPar = [];
		return;
	end;
	model = struct('Evectors', dm.ShapeData.Evectors,...
					'Evalues', dm.ShapeData.Evalues,...
					'x_mean', dm.ShapeData.x_mean,...
					'motion_cov', kron(dm.motion_model_var,dm.motion_model_var') .* dm.motion_model_corr_coef);  

% build masks
	maxFrames = 100;			% max # frames to search for mask construction
	mask = MakeMask(GetContours('GETMOVIEFRAME',[state.CURFRAME min([(state.CURFRAME+maxFrames-1) state.NFRAMES-1])]));
	se = strel('rectangle', [10 10]);
	mask = imopen(mask,se);
	gmask = ComputeImageForces(mask, dPar.SIGMA);
	gmask = gmask./max(gmask(:));	% mask gradient

% include state info
	tPar.MASK = mask;
	tPar.GMASK = gmask;
	tPar.MODEL = model;
	tPar.PFSTATE = [];
	tPar.STARTLENGTH = [];
	tPar.STARTENERGY = [];
end;

% init config
figPos = get(0, 'ScreenSize');
width = 290; height = 390;
figPos = [figPos(1)+(figPos(3)-width)/2, figPos(2)+(figPos(4)-height)/2, width, height];

cfg = dialog('Name', idString, ...
	'tag', 'GETCONTOURS', ...
	'menubar', 'none', ...
	'Position', figPos, ...
	'KeyPressFcn', 'set(gcbf,''UserData'',1);uiresume', ...
	'UserData', 0);

% about
blurb = ['This procedure implements SLURP, a multi- hypothesis snake contour tracking algorithm.'];
if ismac, fs = 12; else; fs = 9; end;
uicontrol(cfg, ...
	'Style', 'frame', ...
	'Units', 'normalized', ...
	'Position', [0.0483 0.8513 0.9000 0.1077]);
uicontrol(cfg, ...
	'Style', 'text', ...
	'HorizontalAlignment', 'left', ...
	'String', blurb, ...
	'FontSize', fs, ...
	'Units', 'normalized', ...
	'Position', [0.0759 0.8641 0.8690 0.0872]);

% colormap
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Colormap:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0241 0.7500 0.3862 0.0500]);
cm = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', tPar.CMAP, ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.4345 0.7500 0.4828 0.0577]);

% sigma
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Sigma:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0724 0.6615 0.3379 0.0500]);
gs = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%.1f',tPar.SIGMA), ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.4345 0.6615 0.4828 0.0577]);
	
% delta
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Delta:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0724 0.5731 0.3379 0.0500]);
ds = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%.1f',tPar.DELTA), ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.4345 0.5731 0.4828 0.0577]);
	
% band penalty
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Band Penalty:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0724 0.4846 0.3379 0.0500]);
bps = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%.1f',tPar.BPEN), ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.4345 0.4846 0.4828 0.0577]);
	
% alpha
cbs1='set(get(gcbo,''userdata''),''string'',sprintf(''%.2f'',get(gcbo,''value'')))';
cbs2='if ~isempty(str2num(get(gcbo,''string'')))&&str2num(get(gcbo,''string''))>=0&&str2num(get(gcbo,''string''))<=1,set(get(gcbo,''userdata''),''value'',str2num(get(gcbo,''string'')));else,set(gcbo,''string'',sprintf(''%0.2f'',get(get(gcbo,''userdata''),''value'')));end';
if ismac, dy = 0; else, dy = .007; end;
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Alpha:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0724 0.3962 0.1931 0.0500]);
als = uicontrol(cfg, ...
	'Style', 'slider', ...
	'Value', tPar.ALPHA, ...
	'Min', 0, ...
	'Max', 1, ...
	'Callback', cbs1, ...
	'Units', 'normalized', ...
	'Position', [0.2897 0.3846+dy 0.4345 0.0577]);
ale = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%0.2f',tPar.ALPHA), ...
	'FontSize', 10, ...
	'Callback', cbs2, ...
	'UserData', als, ...
	'Units', 'normalized', ...
	'Position', [0.7483 0.3962 0.1690 0.0577]);
als.UserData = ale;

% lambda
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Lambda:', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0724 0.3077 0.1931 0.0500]);
las = uicontrol(cfg, ...
	'Style', 'slider', ...
	'Value', tPar.LAMBDA, ...
	'Min', 0, ...
	'Max', 1, ...
	'Callback', cbs1, ...
	'Units', 'normalized', ...
	'Position', [0.2897 0.2962+dy 0.4345 0.0577]);
lae = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%0.2f',tPar.LAMBDA), ...
	'FontSize', 10, ...
	'Callback', cbs2, ...
	'UserData', las, ...
	'Units', 'normalized', ...
	'Position', [0.7483 0.3077 0.1690 0.0577]);
las.UserData = lae;

% adaptive sampling checkbox and label
as = uicontrol(cfg, ...
	'Style', 'checkbox', ...
	'HorizontalAlignment', 'left', ...
	'String', 'Adaptive Sampling', ...
	'FontSize', 10, ...
	'Value', tPar.ADAPTIVE, ...
	'Units', 'normalized', ...
	'Position', [0.3138 0.2192 0.4828 0.0577]);

% min/max particles
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Particles:    Min', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.0241 0.1346 0.4103 0.0500]);
minp = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%.0f',tPar.MINPART), ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.4466 0.1346 0.1690 0.0577]);
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Max', ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.6155 0.1346 0.1207 0.0500]);
maxp = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', sprintf('%.0f',tPar.MAXPART), ...
	'FontSize', 10, ...
	'Units', 'normalized', ...
	'Position', [0.7483 0.1346 0.1690 0.0577]);
	
% OK, Defaults, cancel buttons
if firstTime, es = 'off'; else, es = 'on'; end;		% cancel not permitted first time (initialization)
uicontrol(cfg, ...
	'String','OK', ...
	'FontSize', 10, ...
	'Callback','set(gcbf,''UserData'',1);uiresume', ...
	'Units', 'normalized', ...
	'Position', [0.1517 0.0359 0.2069 0.0641]);
uicontrol(cfg, ...
	'String','Defaults', ...
	'FontSize', 10, ...
	'Callback','set(gcbf,''UserData'',2);uiresume', ...
	'Units', 'normalized', ...
	'Position', [0.3931 0.0359 0.2069 0.0641]);
uicontrol(cfg, ...
	'enable', es, ...
	'String','Cancel', ...
	'FontSize', 10, ...
	'Callback','set(gcbf,''UserData'',0);uiresume', ...
	'Units', 'normalized', ...
	'Position', [0.6345 0.0359 0.2069 0.0641]);

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

	sigma = str2num(get(gs,'string'));
	if isempty(sigma),
		set(gs,'string',sprintf('%.0f',dPar.SIGMA));
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

	minPart = str2num(get(minp,'string'));
	if isempty(minPart),
		set(minp,'string',sprintf('%.0f',dPar.MINPART));
		continue;
	end;

	maxPart = str2num(get(maxp,'string'));
	if isempty(maxPart),
		set(maxp,'string',sprintf('%.0f',dPar.MAXPART));
		continue;
	end;

	switch get(cfg,'userdata'),
		case 0, 		% cancel
			tPar = []; 
			break;
		case 1,			% ok
			tPar.CMAP = cName;
			tPar.SIGMA = sigma;
			tPar.DELTA = delta;
			tPar.BPEN = bpen;
			tPar.ALPHA = get(als,'value');
			tPar.LAMBDA = get(las,'value');
			tPar.ADAPTIVE = get(as,'value');
			tPar.MINPART = minPart;
			tPar.MAXPART = maxPart;
			break;
		case 2,			% defaults
			set(cm,'string',dPar.CMAP);
			set(gs,'string',sprintf('%.0f',dPar.SIGMA));
			set(ds,'string',sprintf('%0.1f',dPar.DELTA));
			set(bps,'string',sprintf('%0.1f',dPar.BPEN));
			set(als,'value',dPar.ALPHA); set(ale,'string',sprintf('%.2f',dPar.ALPHA));
			set(las,'value',dPar.LAMBDA); set(lae,'string',sprintf('%.2f',dPar.LAMBDA));
			set(as,'value',dPar.ADAPTIVE);
			set(minp,'string',sprintf('%.0f',dPar.MINPART));
			set(maxp,'string',sprintf('%.0f',dPar.MAXPART));
			continue;
	end;
end;
if ishandle(cfg), delete(cfg); end;


%=============================================================================
% DOTRACK  - track contours using SLURP method

function [tPar,xy,energy,nParticles] = DoTrack(img, tPar, init_pts, firstTime)

% compute image forces
gimg = ComputeImageForces(img, tPar.SIGMA);

% convolve with mask
gimg = 1 - gimg .* tPar.MASK .* (1 - tPar.GMASK) ./ max(gimg(:));

% ----- initial frame
if firstTime,
	useBE = 1;		% use band energy
	[xy,energy] = make_snake(img', gimg', init_pts, tPar.DELTA*ones(tPar.NPOINTS,1), tPar.BPEN, tPar.ALPHA, tPar.LAMBDA, useBE);

% ensure XY defined increasing left to right (blows up below if right to left)
	if xy(1,1) > xy(end,1), xy = flipud(xy); end;
	
% create particles for subsequent calls
	nParticles = tPar.MAXPART;
	pfstate = zeros(3+size(tPar.MODEL.Evectors,2),nParticles);
	pfstate(1:2,:) = repmat(xy(1,:)',1,nParticles);
	start_length = sum(sqrt(sum(diff(xy).^2,2)));
	normalized_pts = reshape((xy - repmat(pfstate(1:2,1)', length(xy),1)),[],1)./start_length;
	pfstate(3,:) = ones(1,nParticles); 	% length scale wrt
												% original apparent
												% tongue length
	pfstate(4:end,:) = repmat(((normalized_pts - tPar.MODEL.x_mean)'*tPar.MODEL.Evectors)', 1, nParticles);
	tPar.STARTLENGTH = start_length;
	tPar.STARTENERGY = energy;
	
% ----- subsequent frames
else,
	pfstate = tPar.PFSTATE;
	
% evolve particles using state transition model
% CL check whether need to evolve max_particles or nParticles
	rv = mvnrnd(zeros(1,size(pfstate,1)), tPar.MODEL.motion_cov, tPar.MAXPART);
	pfstate(1:2,:) = pfstate(1:2,:) + tPar.STARTLENGTH*repmat(pfstate(3,:),[2 1]).*rv(:,1:2)';
	pfstate(3:end,:) = pfstate(3:end,:) + rv(:,3:end)';

% evaluate new particles using external energy terms in snake model
% CL: max particles or nParticles?
	pf_xy = zeros(tPar.NPOINTS, 2, tPar.MAXPART);
	pxy = zeros(size(pf_xy));
	p_energy = zeros(tPar.NPOINTS, tPar.MAXPART);
	p_total_energy = zeros(1, tPar.MAXPART);

	pp = 1;
	cumlike = 0;
	minlike = 7 * exp(-sum(tPar.STARTENERGY));
	useBE = 0;		% don't use band energy

% for adaptive sampling, cumulate enough particles to reach minimum likelihood threshold while 
% still keeping minimum number of particles and without exceeding maximum number allowed

	while pp < tPar.MINPART || (pp < tPar.MAXPART && (cumlike < minlike || ~tPar.ADAPTIVE)),

% get snake vertices from particle state
		pt_vec = (tPar.MODEL.x_mean + sum(repmat(pfstate(4:end,pp)', length(tPar.MODEL.Evectors),1).*tPar.MODEL.Evectors,2));
		pf_xy(:,:,pp) = reshape(pt_vec,tPar.NPOINTS,2).* pfstate(3,pp).*tPar.STARTLENGTH + repmat(pfstate(1:2,pp)',tPar.NPOINTS,1);  
		[pxy(:,:,pp),p_energy(:,pp)] = make_snake(img', gimg', pf_xy(:,:,pp), tPar.DELTA*ones(tPar.NPOINTS,1), tPar.BPEN, tPar.ALPHA, tPar.LAMBDA, useBE);

% compute particle likelihood
		p_total_energy(pp) = sum(p_energy(:,pp));
		len = sum(sqrt(sum(diff(pxy(:,:,pp)).^2,2)));
		lratio(pp) = max(len/tPar.STARTLENGTH, tPar.STARTLENGTH/len);
		like(pp) = exp(-p_total_energy(pp)*lratio(pp));
		cumlike = cumlike + like(pp);

		pp = pp + 1;
	end;

% add one particle to contain refined best snake result
	nParticles = pp + 1;        
	[max_w, best] = max(like);

% refine best snake using full objective function
	xy = pxy(:,:,best);
	useBE = 1;		% use band energy
	[xy,energy] = make_snake(img', gimg', xy, tPar.DELTA*ones(tPar.NPOINTS,1), tPar.BPEN, tPar.ALPHA, tPar.LAMBDA, useBE);
	len = sum(sqrt(sum(diff(xy).^2, 2)));
	lr = max(len/tPar.STARTLENGTH, tPar.STARTLENGTH/len);

% save updated best state to new particle
	pfstate(:,nParticles) = zeros(size(pfstate,1),1);
	pfstate(1:2,nParticles) = xy(1,:);
	pfstate(3,nParticles) = sum(sqrt(sum(diff(xy).^2,2)))./tPar.STARTLENGTH;
	normalized_pt = reshape((xy(:,:)-repmat(pfstate(1:2,nParticles)',length(xy(:,:)),1)),[],1)./(pfstate(3,nParticles)*tPar.STARTLENGTH);        
	pfstate(4:end,nParticles) = (normalized_pt - tPar.MODEL.x_mean)'*tPar.MODEL.Evectors;
	like(nParticles) = exp(-sum(energy*lr));

% sample new generation of particles to evolve based on particle likelihood (importance sampling)
	weight = like./sum(like);
	cdf = cumsum(weight);
	cdf_prev = circshift(cdf,[0 1]);
	cdf_prev(1) = 0.0;
	samples = rand(tPar.MAXPART, 1);
	for pp = 1 : tPar.MAXPART,           
		ix(pp) = find(cdf >= samples(pp) & cdf_prev <= samples(pp));
	end;
	pfstate = pfstate(:,ix(1:tPar.MAXPART));
	tPar.ENERGY = energy;
end;

% save for next round
tPar.PFSTATE = pfstate;    
tPar.ENERGY = energy;


%=============================================================================
% MAKEMASK  - compute image mask using non-varying regions through frame comparison
%
% find image mask by looking at parts where there is significant variation from one 
% image to the next

function mask = MakeMask(img)

var_img = var(double(img), 0, 3);
vmask = var_img > 2;

CC = bwconncomp(vmask);
numOfPixels = cellfun(@numel,CC.PixelIdxList);
[~,indexOfMax] = max(numOfPixels);
mask = zeros(size(vmask));
mask(CC.PixelIdxList{indexOfMax}) = 1; 


%=============================================================================
% SETFRAME  - DataCursorManager UpdateFcn that sets current frame based on
%				clicked position in diagnostic map

function txt = SetFrame(h, evt)

pos = round(evt.Position);		% X,Y
ud = get(gca,'userdata');
fh = get(gca,'parent');
doFlag = get(fh,'userdata');
set(fh,'userData',~doFlag);
if ud{1} == 1,
	ih = get(gca,'children');
	img = ih.CData;
	energy = img(pos(2),pos(1)) - 1;	% report true energy value
	frame = pos(1) + ud{2} - 1;		% add frame offset
	yy = ud{3};
	pt = round(yy(pos(2)));			% interpolate to # anchors
	txt = {['Frame: ',num2str(frame)],['Point: ',num2str(pt),'  Energy: ',sprintf('%.2f',energy)]};
else,
	frame = pos(1);
	txt = {['Frame: ',num2str(frame)],['nParticles: ',num2str(pos(2))]};
end;
if doFlag, GetContours('FRAME','EXPLICIT',frame); end;
