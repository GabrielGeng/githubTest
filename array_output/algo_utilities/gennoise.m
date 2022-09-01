function [corruptsignal, newnoise, confirmSNR]= gennoise(signal, noise, SNR);
%   Generates a new corrupted signal based on a given SNR.
%
%   Generates the new noise sequence as output.
%   New SNR is confirmed via confirmSNR.
%
%   Usage: [corruptsignal, newnoise, confirmSNR]= gennoise(signal, noise, SNR);
%
%   Inputs          signal  : clean signal
%                   noise   : random noise of unit variance
%                   SNR     : the desired SNR
%
%   Outputs         corruptsignal   : noisy signal
%                   newnoise        : amount of noise added to clean signal
%                   confirmSNR      : SNR (for checking)
%
% Copyright (C) Andy Khong 2003
%
% Last modified 19th Dec 2003
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



signalpow=(signal'*signal)/length(signal);
noisepow= (noise'*noise)/length(noise);

newnoise= sqrt((signalpow/((10^(SNR*0.1))*noisepow)))*noise;

corruptsignal= signal+newnoise;
confirmSNR= 10*log10(signalpow/((newnoise'*newnoise)/length(noise)));