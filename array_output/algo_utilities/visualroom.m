function visualroom(source, receiver, room);
%   Visualize a synthetic room, microphone and source placement.
%
%   Usage:  visualroom(source, receiver, room);
%
%   Inputs      source      : source positions (srcNum-by-3)
%               receiver    : receiver positions (recNum-by-3)
%               room        : room dimensions (1-by-3)
%
% Copyright (C) Andy Khong 2004
%
% Last modified 12th Jul 2003
%               24th May 2004
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


figure; hold on;
axis([0 room(1) 0 room(2) 0 room(3)]);
source_size     = size(source);
receiver_size   = size(receiver);
source_number   = source_size(1);
receiver_number = receiver_size(1);

for i=1: receiver_number;    
    len     = length(0:0.01:receiver(i,3));
    recX    = receiver(i,1)*ones(len,1);
    recY    = receiver(i,2)*ones(len,1);
    recH    = [0:0.01:receiver(i,3)];
    plot3(recX,recY,recH,'LineWidth',2,'Color','k');
    plot3(recX(end),recY(end),recH(end),'ko','MarkerSize',10);
end

for i=1: source_number;    
    len     = length(0:0.01:source(i,3));
    socX    = source(i,1)*ones(len,1);
    socY    = source(i,2)*ones(len,1);
    socH    = [0:0.01:source(i,3)];
    plot3(socX,socY,socH,'LineWidth',2,'Color','k');
    plot3(socX(end),socY(end),socH(end),'k*','MarkerSize',10);
end

box on;
grid on;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');








