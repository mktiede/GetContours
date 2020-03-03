function ReshapeContours(fName, EdgeTrak)
%RESHAPECONTOURS  - convert long format GetContours output to wide format
%
%	usage:  ReshapeContours(fName, EdgeTrak)
%
% GETCONTOURS exports contour data to a tab-delimited text file in long format
% i.e., each line specifies for one coordinate of each contour the following
%
%   FRAME  TIME  NOTE  POINT  X  Y  (mmX  mmY)
%
% where FRAME is source movie frame number, TIME is the corresponding temporal offset (secs),
% NOTE is any annotation associated with the frame, POINT is the coordinate number, 
% X and Y are the POINT pixel coordinates (relative to image ULC), and optionally 
% (if origin and mm/pixel values defined) origin-based coordinates in mm
%
% this procedure loads file FNAME (default extension ".tsv") in this format and creates a
% corresponding wide format file output FNAME_wide.tsv in which each row has tab-delimited
% format FRAME TIME NOTE X1 Y1 X2 Y2 ... Xn Yn
%
% if mm data exist a second file is also created as FRAME_wide_mm.tsv with mm coordinates
%
% if optional EDGETRAK is nonzero output is in EdgeTrak format; i.e. FNAME_EdgeTrak.tsv
% and (if mm data exist) FNAME_EdgeTrak_mm.txt with format
%   FRAME1  FRAME2 ... FRAMEn
%   X   Y   X   Y      X   Y
%   (nPoints rows of contour data)

% mkt 02/20

if nargin < 1, eval('help ReshapeContours'); return; end;
if nargin < 2 || isempty(EdgeTrak), EdgeTrak = 0; end;

% slurp file
[p,fName,ext] = fileparts(fName);
if isempty(ext), ext = '.tsv'; end;
ifn = fullfile(p,[fName,ext]);

fid = fopen(ifn,'rt');
if fid == -1, error('unable to open %s', ifn); end;
lines = {};
while 1,
	lx  = fgetl(fid);
	if ~ischar(lx), break; end;
	lines{end+1,1} = lx;
end;
fclose(fid);

% is this a GC file?
hdr = lines{1};
lines(1) = [];
hdr = regexp(hdr,'\t','split');
if length(hdr) < 5 || ~strcmp(hdr{1},'FRAME') || ~strcmp(hdr{2},'TIME'),
	error('%s does not appear to be a GetContours export file', fName);
end;

% does it have mm coordinates?
gotmm = (length(hdr) == 8 && strcmp(hdr{8},'mmY'));

% init output files
if EdgeTrak,
	ofn = fullfile(p,sprintf('%s_EdgeTrak%s',fName,ext));
else,
	ofn = fullfile(p,sprintf('%s_wide%s',fName,ext));
end;
if exist(ofn) == 2 && strcmp(questdlg(sprintf('Overwrite %s?',ofn), 'File exists...', 'Yes', 'No', 'Yes'), 'No'),
	fprintf('Conversion of %s cancelled.\n', ifn);
	return;
end;
fid = fopen(ofn,'wt');
if fid == -1, error('unable to open %s for writing', ofn); end;
if gotmm, 
	if EdgeTrak,
		ofnmm = fullfile(p,sprintf('%s_EdgeTrak_mm%s',fName,ext));
	else,
		ofnmm = fullfile(p,sprintf('%s_wide_mm%s',fName,ext));
	end;
	fidmm = fopen(ofnmm,'wt');
end;

% find number of points per contour
q = regexp(lines{1},'\t','split');
f = str2num(q{1});
nPts = 1;
while 1,
	q = regexp(lines{nPts+1},'\t','split');
	n = str2num(q{1});
	if n > f, break; end;
	nPts = nPts + 1;
end;

% find number of frames
nFrames = length(lines) / nPts;

% load the contour data (for EdgeTrak option)
if EdgeTrak,
	c = zeros(nPts,2+gotmm*2,nFrames);
	f = zeros(nFrames,1);
	li = 1;
	for fi = 1 : nFrames,
		for pi = 1 : nPts,
			q = regexp(lines{li},'\t','split');
			if pi == 1, f(fi) = str2num(q{1}); end;
			c(pi,1,fi) = str2num(q{5});
			c(pi,2,fi) = str2num(q{6});
			if gotmm,
				c(pi,3,fi) = str2num(q{7});
				c(pi,4,fi) = str2num(q{8});
			end;
			li = li + 1;
		end;
	end;
end;

% write header
if EdgeTrak,
	fprintf(fid,'FRAME%d',f(1)); 
	if gotmm, fprintf(fidmm,'FRAME%d',f(1)); end;
	for k = 2 : nFrames, 
		fprintf(fid,'\t\tFRAME%d',f(k)); 
		if gotmm, fprintf(fidmm,'\t\tFRAME%d',f(k)); end;
	end;
	fprintf(fid,'\nX\tY');
	if gotmm, fprintf(fidmm,'\nX\tY'); end;
	for k = 2 : nFrames, 
		fprintf(fid,'\tX\tY');
		if gotmm, fprintf(fidmm,'\tX\tY'); end;
	end;
	fprintf(fid,'\n');
	if gotmm, fprintf(fidmm,'\n'); end;
else,
	hdr = hdr(1:3);
	hdrmm = hdr;
	for k = 1 : nPts, 
		hdr{end+1} = sprintf('X%03d',k); hdrmm{end+1} = sprintf('mmX%03d',k);
		hdr{end+1} = sprintf('Y%03d',k); hdrmm{end+1} = sprintf('mmY%03d',k);
	end;
	fprintf(fid,'%s\t',hdr{1:end-1});
	fprintf(fid,'%s\n',hdr{end});
	if gotmm,
		fprintf(fidmm,'%s\t',hdrmm{1:end-1});
		fprintf(fidmm,'%s\n',hdrmm{end});
	end;
end;

% write data
if EdgeTrak,
	for pi = 1 : nPts,
		fprintf(fid,'%.2f\t%.2f',c(pi,1:2,1));
		if gotmm, fprintf(fidmm,'%.2f\t%.2f',c(pi,3:4,1)); end;
		for fi = 2 : nFrames,
			fprintf(fid,'\t%.2f\t%.2f',c(pi,1:2,fi));
			if gotmm, fprintf(fidmm,'\t%.2f\t%.2f',c(pi,3:4,fi)); end;
		end;
		fprintf(fid,'\n');
		if gotmm, fprintf(fidmm,'\n'); end;
	end;
else,
	li = 1;
	while li < length(lines),
	
% get the frame-invariant data
		q = regexp(lines{li},'\t','split');
		frame = str2num(q{1});
		time = str2num(q{2});
		note = q{3};
		if gotmm,
			c = NaN(nPts,4);
		else,
			c = NaN(nPts,2);
		end;

% load the contour
		for k = 1 : nPts
			q = regexp(lines{li},'\t','split');
			c(k,1) = str2num(q{5});
			c(k,2) = str2num(q{6});
			if gotmm,
				c(k,3) = str2num(q{7});
				c(k,4) = str2num(q{8});
			end;
			li = li + 1;
		end;

% write the contour
		v = c(:,1:2)';
		fprintf(fid,'%d\t%f\t%s', frame, time, note);
		fprintf(fid,'\t%.2f',v(:));
		fprintf(fid,'\n');
		if gotmm,
			fprintf(fidmm,'%d\t%f\t%s', frame, time, note);
			v = c(:,3:4)';
			fprintf(fidmm,'\t%.2f',v(:));
			fprintf(fidmm,'\n');
		end;
	end;
end;

% clean up
fclose(fid);
if gotmm,
	fclose(fidmm);
	fprintf('%s converted to %s and %s\n', ifn, ofn, ofnmm);
else,
	fprintf('%s converted to %s\n', ifn, ofn);
end;
