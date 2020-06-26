function [z_realiz] = f_her4_sHER_NOnugget(x_cal, y_cal,z_cal, x_target, y_target, her)
%% function to generate a realization using Sequential Simulation and HER

% -------------- Input --------------
% - x_cal       [N,1]   x coordidates of the calibration set
% - y_cal       [N,1]   y coordidates of the calibration set
% - z_cal       [N,1]   z values of the calibration set (variable under study)
% - x_target    [T,1]   x coordidates of the target set
% - y_target    [T,1]   y coordidates of the target set
% - her         struct  structure cointaining model definitions

% -------------- Output --------------
% - z_realiz    [1,T]    realization of the targets 

% -------------- Version --------------
% - 2020/06/16 Stephanie Thiesen: intial version

% -------------- Script --------------
    z_realiz = NaN(1,length(x_target));
    idx_ = randperm(length(x_target));
    for i = 1:length(idx_)
        [pmf, ~] = f_her3_predict_NOnugget(x_cal, y_cal, z_cal, x_target(idx_(i)), y_target(idx_(i)), her);
        z_realiz(idx_(i)) = f_montecarlosim(cell2mat(pmf), her.edges_z); %Monte carlo Simulation
        x_cal = [x_cal; x_target(idx_(i))]; %sequentially including the simulated value as data
        y_cal = [y_cal; y_target(idx_(i))];
        z_cal = [z_cal; z_realiz(idx_(i))];
    end  
   
end