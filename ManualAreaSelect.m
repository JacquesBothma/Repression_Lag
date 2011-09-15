%%                              ManualAreaSelect.m
% Jacques Bothma                                    Last Modified: 09/13/11     
% Levine Lab, UC Berkeley                       
% Functionally complete                          
%
%
%% Attribution:
%  Feel free to use, modify and distribute this code provided that you
%  attribute Jacques Bothma for development.
%  This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
%  To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.
%%
%
% Overview:
%
% This code allows manual selection of region by user using a mouse of an
% arbirary number of regions. 
%
% Inputs:
%
% Rs -  Number of sepearte regions to select
% Nuc - Image that is the same size as the one being used to define regions
%
% Outputs:
% 
% Mask - Binary imnage that has the selected regions as ones and background region as zeros. 
%

function Mask=ManualAreaSelect(Rs,Nuc)

%figs = get(0, 'Children') %List of open figures 
%imshowbig(label2rgb(Nuc,'jet',[0,0,0]))
%figure(figs(2)); %Go to the figure with the highest handle, typically the one last active
imshowbig(Nuc);


hold on

for k=1:Rs % Loop, picking up the points.
        
xy = []; % initially the list of points is empty
n = 0; % initiate number of points

    but = 1;
    while but == 1
        [xi,yi,but] = ginput(1);
        plot(xi,yi,'wo')
        n = n+1;
        xy(:,n) = [xi;yi];
    end
    % Interpolate with a smooth (spline) curve and finer spacing.
    ts = linspace(1,n,150);
    xys = spline(1:n,xy,ts);

    % Plot the interpolated curve and export data
    bnd_x{k} = [xys(1,:), xys(1,1)]; 
    bnd_y{k} = [xys(2,:), xys(2,1)];
    bndrys(k).ConvexHull = [bnd_x{k}',bnd_y{k}'];
    
    
    
    plot(bnd_x{k},bnd_y{k},'m-');
    
end

[aa,bb,cc]=size(Nuc);

BW=zeros([aa,bb]);

for k=1:Rs
    
BW = poly2mask(bnd_x{k}, bnd_y{k}, aa, bb) + BW; % create mask that correspond to perimeter defined by the previous proc.

end

BW=logical(BW);

hold off

Mask=BW;
