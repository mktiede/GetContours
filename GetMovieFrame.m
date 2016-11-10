function img = GetMovieFrame(mh, frame)
% GETMOVIEFRAME  - get frame from open video object
%
%	usage:  img = GetMovieFrame(mh, frame)
%
% returns video image associated with 1-based FRAME from open movie handle MH

% mkt 09/15

if verLessThan('matlab','8.5.0'),
	try, img = read(mh, frame); catch, error('unable to read frame %d from specified movie handle',frame); end;
else,
	mh.CurrentTime = (frame-1)/mh.FrameRate;
	try, img = readFrame(mh); catch, error('unable to read frame %d from specified movie handle',frame); end;
end;
