function [z_sim] = f_montecarlosim(pmf, edges_z)
%% function to perform Monte Carlo Simulation from PMF
% randomly draw a p-value uniformly distributed between 0 and 1 and obtain the z value from the CCDF

% -------------- Input --------------
% - pmf         [1,n]     predicted z PMF for a target 
% - edges_z     [1,n+1]   bin edges of the z PMF

% -------------- Output --------------
% - z_sim       s         simulated z-value 

% -------------- Version --------------
% - 2020/06/16 Stephanie Thiesen: intial version

% -------------- Script --------------
    p_value = rand; %single uniformly distributed random number in the interval (0,1)
    ccdf = [0 cumsum(pmf)]; %cummulative PMF
    % z_sim = interp1(ccdf,edges_z,p_value);
    idx_ = find(ccdf >= p_value, 1); %find where ccdf >= p_value (index of the 
    z_sim = interp1(ccdf([idx_-1 idx_]),edges_z([idx_-1 idx_]),p_value);
    
end