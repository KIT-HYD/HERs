function [z_realiz] = f_her4_HERs_NOnugget(x_cal, y_cal,z_cal, x_target, y_target, her)
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
    x_cal_ = x_cal;
    y_cal_ = y_cal;
    z_cal_ = z_cal;
    for i = 1:length(idx_)
        for j = 1:length(x_cal)
            euc_distance_(j) = f_euclidean_dist(x_cal(j), y_cal(j), x_target(i), y_target(i)); %calculate the euclidean distance between cal. set and target
        end
        [dist_, idx_cal_] = min(euc_distance_);
        if dist_ <= 10e-6 %honoring calibration set
            z_realiz(idx_(i)) = z_cal(idx_cal_);
        else
        [pmf, ~] = f_her3_predict_NOnugget(x_cal_, y_cal_, z_cal_, x_target(idx_(i)), y_target(idx_(i)), her);
        z_realiz(idx_(i)) = f_montecarlosim(cell2mat(pmf), her.edges_z); %Monte carlo Simulation
        x_cal_ = [x_cal_; x_target(idx_(i))]; %sequentially including the simulated value as data
        y_cal_ = [y_cal_; y_target(idx_(i))];
        z_cal_ = [z_cal_; z_realiz(idx_(i))];
        end
    end  
   
end