function [I] = f_indicator(z_thresh, ObsA, ObsB)
%% function for calculating indicator
% -------------- Input -------------- 

% -------------- Output --------------
% - I            [1,1]                indicator value with classification
%                                     according to a threshold. (1 if 
%                                     observation (row), target (column)

% -------------- Version --------------
% - 2020/03/30 Stephanie Thiesen: intial version

% -------------- Script --------------
ObsA_I = double(ObsA <= z_thresh);
ObsB_I = double(ObsB <= z_thresh);

I = (ObsA_I - ObsB_I)^2;

end