function [xframe,st,ed]= chopsignal(frIdx,x,frLen,ovrLapFac,noFr);
% Extracts a frame of signal given a frame index.
%
% Usage:    [xframe,st,ed]= chopsignal(frIdx,x,frLen,ovrLapFac,noFr);
%
% Inputs:   frIdx       = frame index
%           x           = signal matrix (length-by-NoOfChannels)
%           frLen       = frame size (in samples)
%           ovrLapFac   = overlapping factor
%           noFr        = total no. of frames computed using compFrNum.m
% 
% Outputs:  xframe      = the extracted signal of a particular frame
%           st          = start of index from which signal is extracted
%           ed          = end of index from which signal is extracted
%
% Copyright (C) Andy Khong 2009
%
% Last modified 24th March 2009 1539 Hrs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   ftp://prep.ai.mit.edu/pub/gnu/COPYING-2.0 or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xframe = zeros(frLen,size(x,2));

shft = (1-ovrLapFac)*frLen;

if frIdx == 1;
    xframe = x(1:frLen,:);
    st = 1;
    ed = frLen;
else
    st = (frIdx-1)*shft+1;
    ed = st+frLen-1;
    xframe = x(st:ed,:);
end
    

 