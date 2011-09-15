%%                               PairmRNAs.m
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
%
% Overview:
%
% This code identifies mRNA pairs and isolted mRNAs from the input of two
% mRNA channels. It does this by rotating the matrices with the locations
% of the different mRNAs relative to each other and looking for overlap.
%
% The inputs are the logical array of mRNA1 and mRNA2 center pixels and the size of the neighborhood in pixels. The
% outputs are the logical arrays of the mRNA pairs, isolated mRNA1s and isolated mRNA2s.
%
% Input:
%
% mRNA1 - logical array of mRNA1 centers
% mRNA2 - logical array of mRNA2 centers
% n - number of pixel shells to consider
%
% Output:
%
% mRNA1o - logical array of isolated mRNA1 centers
% mRNA2o - logical array of isolated mRNA2 centers
% mRNAp - logical array of paired mRNA centers
%
%
%

function [mRNA1o, mRNA2o, mRNAp] =PairmRNAs(mRNA1,mRNA2,n)


%%% Generate an array that contains the xy shifts neccessary to pair
%%% neighbouring mRNAs.

S=[0 0];       % First element that corresponds to perfect colocalization
R=[0,-1;1,0];  % 90 degree rotation matrix

for i=1:n
    
t=i*ones(2*i,1);
t2=[-i+1:i]';
T=[t,t2]; 
ST=[T;T*R;T*R^2;T*R^3]; % Generate the shifts that cover the i'th shell
S=[S;ST];               % Combine the shifts of the i'th shell with those of all the smaller ones

end


%%%%%



[a,b] = size(mRNA1); % Size of matrix

mRNA1o=sparse(a,b); %defining arrays
mRNA2o=sparse(a,b);
mRNAp=sparse(a,b);
mRNApT=sparse(a,b);

mRNA1=sparse(mRNA1); % redefining inputs as sparse
mRNA2=sparse(mRNA2); %



ls=length(S);

for i=1:ls   
mRNApT=mRNA1.*circshift(mRNA2,S(i,:)); % Shift mRNA2 relative to mRNA1, where overlap exists define pairs
mRNA1=mRNA1-mRNApT;                    % Remove pixel in mRNA1 that is associated with pair
mRNA2=mRNA2-circshift(mRNApT,-S(i,:)); % Remove pixel in mRNA2 that is associated with pair
mRNAp=mRNAp+mRNApT;                    % Add newly found pairs to pair matrix
end


mRNA1o=mRNA1;
mRNA2o=mRNA2;

 mRNA1o=full(mRNA1o);
 mRNA2o=full(mRNA2o);
 mRNAp=full(mRNAp);
