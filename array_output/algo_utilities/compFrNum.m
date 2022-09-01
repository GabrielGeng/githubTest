function noFr= compFrNum(len,frLen,ovrLapFac);
% Computes the number of frames given the total length of the signal, frame size and overlapping factor.
%
% Usage:    noFr= compFrNum(len,frLen,ovrLapFac);
%
% inputs:   len         = length of the signal
%           frLen       = frame size (in samples)
%           ovrLapFac   = overlapping factor
% 
% outputs:  noFr        = number of frames
%
% Copyright (C) Andy Khong 2009
%
% Last modified 24th Jan 2009
%               24th Mar 2009- added in line 48: 1-
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


%test vec
% len = 10;
% ovrLapFac = 0.5;
% frSze= 4;

if ovrLapFac ==0;
    noFr = floor(len/frLen);
else
    if frLen==1;
        noFr = len;
    else
        shftSamp    = (1-ovrLapFac)*frLen;          % number of samples per shift
        lftOvr      = len-frLen;                % left over samples after first frame
        lftOvrBloc  = floor(lftOvr/shftSamp);   % there will be this no. of shifts for the remaining samples (each shift will be a frame)    
        noFr        = lftOvrBloc+1;             % now add the first frame
    end
end



