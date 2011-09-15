%%                             imshowbig.m
% Jacques Bothma                                      
% Levine Lab, UC Berkeley                        
% Functionally complete                             Last Modified: 09/13/10
%
%% Attribution:
% Feel free to use, modify and distribute this code provided that you
% attribute Jacques Bothma for development.
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.
%
%  Important Notes:
% % This version written for Windows.
%
% Overview:
%
% This code simply plots images usiong imshow and makes them as big as the size of the screen allows. Uses maxwindow from matlab file exchange.
%
% Input:
%
% F -  image to be plotted compatible with imshow plotting


function imshowbig(F)
figure, imshow(F,'Border','tight','InitialMagnification',100), maxwindow
%set(gcf,'Visible','off')
end
