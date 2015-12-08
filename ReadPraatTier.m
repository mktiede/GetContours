function [segs,labs,tier,tierNames,dur] = ReadPraatTier(fName,tier)
%READPRAATTIER  - read tier data from PRAAT TextGrid file
%
%	usage:  [segs,labs,tier,tierNames,dur] = ReadPraatTier(fName, tier)
%
% use this procedure to load the specified point or interval TIER data 
% from PRAAT TextGrid format FNAME (both standard and short forms supported)
%
% default ".TextGrid" extension is optional
%
% if no return arguements are specified returns a list of available tier names
%
% TIER may be a tier name or numeric index; if none specified the first tier found is used
%
% returns SEGS offsets [nSegs x offs] and segment labels LABS [nSegs]
% for specified tier (for interval tiers offs == head,tail)
%
% optionally returns name of tier used (TIER), a list of all tiers (TIERNAMES),
% and DURation (secs)
%
% see also WRITEPRAATTEXTGRID

% mkt 03/09
% mkt 07/14 updated

% parse args
if nargin < 1,
	eval('help ReadPraatTier');
	return;
end;
if nargin<2, tier = []; end
if isempty(tier), tier = 1; end

[p,f,e] = fileparts(fName);
if isempty(e), fName = fullfile(p,[f,'.TextGrid']); end;

% vacuum file
try,
	fid = fopen(fName,'rt');
	lines = {};
	while 1,
		line = fgetl(fid);
		if ~ischar(line), break; end;
		lines{end+1,1} = line;
	end;
	fclose(fid);
catch,
	error('error attempting to load from %s', fName);
end;

% some rudimentary format checking
if length(lines)<15 || isempty(findstr(lines{1},'ooTextFile')) || isempty(findstr(lines{2},'"TextGrid"')) || isempty(findstr(lines{6},'exists')),
	error('%s has unrecognized file format', fName);
end;

% long or short?
try,
	data = str2num(char(lines([4 5 7])));
	isShort = 1;
	if isempty(data),
		data(1) = str2num(lines{4}(8:end));
		data(2) = str2num(lines{5}(8:end));
		data(3) = str2num(lines{7}(8:end));
		isShort = 0;
	end;
catch,
	error('%s has unrecognized file format', fName);
end;
head = data(1);
tail = data(2);
dur = tail;
nTiers = data(3);

% get tier names and types
pointTiers = []; intTiers = [];
if isShort,
	k = strmatch('"TextTier"',lines,'exact');
	if isempty(k),
		pointTierNames = {};
	else,
		pointTiers = k + 1;
		pointTierNames = lines(pointTiers);
		for k = 1 : length(pointTiers), pointTierNames{k} = pointTierNames{k}(2:end-1); end;
	end;
	k = strmatch('"IntervalTier"',lines,'exact');
	if isempty(k),
		intTierNames = {};
	else,
		intTiers = k + 1;
		intTierNames = lines(intTiers);
		for k = 1 : length(intTiers), intTierNames{k} = intTierNames{k}(2:end-1); end;
	end;
else,
	k = find(~cellfun(@isempty,regexp(lines,'class = "TextTier"')));
	if isempty(k),
		pointTierNames = {};
	else,
		pointTiers = k + 1;
		q = regexp(lines(k+1),'name = "(.+)"','tokens');
		for k = 1 : length(q), pointTierNames(k) = q{k}{1}; end;
	end;
	k = find(~cellfun(@isempty,regexp(lines,'class = "IntervalTier"')));
	if isempty(k),
		intTierNames = {};
	else,
		intTiers = k + 1;
		q = regexp(lines(k+1),'name = "(.+)"','tokens');
		for k = 1 : length(q), 
			if isempty(q{k}), intTierNames(k) = {''}; else, intTierNames(k) = q{k}{1}; end;
		end;
	end;
end;

% list names if no output args
tierNames = [pointTierNames , intTierNames];
if nTiers ~= length(tierNames),
	error('expecting %d tiers, found %d', nTiers, length(tierNames));
end;
if nargout < 1,
	fprintf('Tiers in %s:  "%s"', fName, tierNames{1});
	if nTiers>1, fprintf(', "%s"', tierNames{2:end}); end
	fprintf('\n');
	return;
end;

% format data
if isnumeric(tier),
	tiers = [];
	if ~isempty(pointTiers), tiers = [tiers, pointTiers]; end;
	if ~isempty(intTiers), tiers = [tiers , intTiers]; end;
	[tiers,k] = sort(tiers);
	tierNames = tierNames(k);
	tierStart = tiers(tier);
	tier = tierNames(tier);
	isPoint = ~isempty(strmatch(tier, pointTierNames));
else,
	k = strmatch(tier, tierNames, 'exact');
	if length(k) > 1,
		k = k(1);
		fprintf('more than one instance of %s found in %s; using first\n', tierNames{k}, fName);
	elseif isempty(k),
		k = strmatch(tier, tierNames);
		if isempty(k),
			error('%s not found in %s', tier, fName);
		else,
			fprintf('exact match to "%s" not found in %s; using "%s"\n', tier, fName, tierNames{1});
			tier = tierNames{1};
		end;
	end;
	k = strmatch(tier, pointTierNames, 'exact');
	if isempty(k),
		k = strmatch(tier, intTierNames, 'exact');
		tiers = intTiers;
		isPoint = 0;
	else,
		tiers = pointTiers;
		isPoint = 1;
	end;
	tierStart = tiers(k);
end;

% verify tier length against file length
if isShort,
	h = str2num(lines{tierStart+1});
	t = str2num(lines{tierStart+2});
else,
	h = str2num(lines{tierStart+1}(regexp(lines{tierStart+1},' \d'):end));
	t = str2num(lines{tierStart+2}(regexp(lines{tierStart+2},' \d'):end));
end;
if head ~= h || tail ~= t,
	error('mismatch between tier (%.1f:%.1f) and file lengths (%.1f:%.1f)', h,t,head,tail);
end;

% short format
urLines = lines;
if isShort,
	nSegs = str2num(lines{tierStart+3});
	M = (3 - isPoint);
	h = tierStart + 4;
	t = nSegs*M + h - 1;
	n = (t - h + 1) / M;
	lines = reshape(lines(h:t),M,n)';
	labs = lines(:,M);
	for k = 1 : length(labs),
		labs{k} = labs{k}(2:end-1);
	end;
	lines = reshape(lines(:,1:(M-1))',(M-1)*n,1);
	segs = reshape(str2num(char(lines)),(M-1),n)';

% long format
else,
	nSegs = str2num(lines{tierStart+3}(regexp(lines{tierStart+3},'\d'):end));
	M = (4 - isPoint);
	h = tierStart + 4;
	t = nSegs*M + h - 1;
	n = (t - h + 1) / M;
	lines = reshape(lines(h:t),M,n)';
	labs = lines(:,M);
	for k = 1 : length(labs),
		kk = findstr(labs{k},'"');
		if length(kk) ~= 2,
			n = strmatch(labs{k},urLines);
			error('tier name without closing " character in line %d of %s (unable to continue)', n, fName);
		end;
		labs{k} = labs{k}(kk(1)+1:kk(2)-1);
	end;
	lines = reshape(lines(:,2:(M-1))',(M-2)*n,1);
	idx = cell2mat(regexp(lines,' \d'));
	segs = zeros(length(lines),1);
	for k = 1 : length(lines),
		segs(k) = str2num(lines{k}(idx(k):end));
	end;
	segs = reshape(segs,(M-2),n)';
end;
