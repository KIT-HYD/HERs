%% HER - Histogram via entropy reduction 
%   HER1: Spatial correlation structure ("her" structure cointaining model definitions)
%   HER2: Optimization ("her")
%   HER3: Prediction ("pmf_pred" structure containing the pmf predictions)
%   HERs4: Sequential simulation ("pmf_realiz" structure containing realizations)

% Version
% 03: Version using the DKL optimization for one single bin, continuous 
% weights for AND_OR optimization +  isotropic nugget
% 03.1: Version with steep slope for the 1st class asymptotically 
% increasing the weight towards infinity as the distance approaches zero (no nugget)

%% Load dataset 
clear all;
close all;
clc; 
addpath('functions');
addpath('1.0.raw data');

load Jura_dataset_Pb_log10.mat %1- 'Argovian'; 2- 'Kimmeridgian'; 3- 'Sequanian'; 4- 'Portlandian'; 5- 'Quaternary'
shp_basin = shaperead('jura_basin.shp');
her.txt = txt; %dataset identification

%to avoid numerical problems
z = round(z,5);
z_cal = round(z_cal,5);
z_val = round(z_val,5);
her.z_thresh = round(her.z_thresh,5);
%% Define infogram properties and aggregation method
her.lag_dist = 0.07; % DEFINE the class size (lag width in meters if the input x,y is given in meters)
her.binwidth = 0.015; % DEFINE delta z binwidth
her.binwidth_z = 0.015; % DEFINE binwidth for the z PMFs
her.n_neighbor = 9; % DEFINE the maximum number of neighbord to be considered in the interpolation and AND_OR optimization (n_neighbor-1) (min 1, max dim_cal)
her.z_thresh = log10(50); %DEFINE z threshold for the probability map. [] if no one is desired

%to be used in HER3 step (no need to run HER1 and HER2 parts when you change it)
her.aggregation_type = 'andor';  %DEFINE the aggregation method

%% HER1: Spatial Characterization
her = f_her1_infogram(x_cal, y_cal, z_cal, her);

% plot Infogram cloud, Infogram & Number of pairs, Histogram of delta_z classes, normalized infogram + variogram
f_plot_infogram(her); %todo bin center

%% HER2: Weight optimization 
% (leave-one-out cross validation on calibration dataset)
her = f_her2_weight_NOnugget(x_cal, y_cal, z_cal, her);

% plot optimum weight
f_plot_weights(her);

%% HER3: z PMF prediction | Local uncertainty
% GRID for ploting results
[pmf_pred_grid.x_target,pmf_pred_grid.y_target] = meshgrid(0:0.05:6, 0:0.05:6); %DEFINE a GRID for predicting and ploting results
pmf_pred_grid.x_target = round(pmf_pred_grid.x_target,2); pmf_pred_grid.y_target = round(pmf_pred_grid.y_target,2);
[pmf_pred_grid.pmf_pred, pmf_pred_grid.idx_zero_neigh_pred] = f_her3_predict_NOnugget(x_cal, y_cal, z_cal, pmf_pred_grid.x_target(:), pmf_pred_grid.y_target(:), her);

%% HERs4: Sequential simulation | Spatial uncertainty
% Monte Carlo simulation + sequential simulation using HER model

% GRID for ploting results
pmf_simul_grid.nsim = 1; %define the number of simulations
[pmf_simul_grid.x_target,pmf_simul_grid.y_target] = meshgrid(0:0.05:6, 0:0.05:6); %DEFINE a GRID for predicting and ploting results
pmf_simul_grid.x_target = round(pmf_simul_grid.x_target,2); pmf_simul_grid.y_target = round(pmf_simul_grid.y_target,2);
    for i = 1:pmf_simul_grid.nsim
        z_realiz_ = f_her4_HERs_NOnugget(x_cal, y_cal, z_cal, pmf_simul_grid.x_target(:), pmf_simul_grid.y_target(:), her);
        pmf_simul_grid.z_realiz(:,:,i) =reshape(z_realiz_, size(pmf_simul_grid.x_target,1), size(pmf_simul_grid.y_target,2));
        i
    end

clear z_realiz_ i

%% save
filename = sprintf('HERs_E03_1_lag%s_bw%s_%s_nn%s_%s_dataset_%s_Ralf.mat', num2str(round(her.lag_dist,2)), ...
    num2str(round(her.binwidth,4)), num2str(round(her.binwidth_z,4)), num2str(round(her.n_neighbor,0)),...
    her.aggregation_type, her.txt)

save(filename); 
