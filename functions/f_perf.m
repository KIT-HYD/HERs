function [perf_val] = f_perf(z_val, z_pred_, prob_thres_pred_, z_thresh, marginal_prob_thres_cal, prob_thres, idx_marg)
    %% function for calculating performance metrics of non-HER models
    % - Using leave-one-out cross-validation, the goal is to calculate the 
    % mean DKL of the left-out targets using OR aggregation method.

    % -------------- Input -------------- 
    % - z_true_            [T,1]        true values of the z target 
    % - z_pred_            [T,1]        predicted value(mean, median, mode) of the z target  
    % - prob_thres_pred_   [T,1]        probability of the target exceeding the z threshold 

    % -------------- Output --------------
    % - perf_val          struct    Structure containing deterministic and probabilistic performance metrics

    % -------------- Version --------------
    % - 2020/04/13 Stephanie Thiesen: intial version

    % -------------- Script --------------
    perf_val.marginal_prob_thres_cal = marginal_prob_thres_cal;
    perf_val.prob_thres = prob_thres;
    perf_val.index_marg = idx_marg;
    % Root mean square deviation RMSD (RMSE)
    % Mean Error (ME)  
    % Mean Absolute Error (MAE) L1 norm, robust parameter estimator
    % Nash-Sutcliffe model efficiency (r2, coefficient of determination)
    [perf_val.error_sign, perf_val.RMSE, perf_val.ME, perf_val.MAE, perf_val.NSE, perf_val.correl_true_pred, perf_val.correl_true_residue] = f_performance_det(z_pred_,z_val);

% %     % scoring rule - DKL in relation to the mean value
% %     perf_val.DKL_score_mean = f_performance_prob(z_true_', PMF_pred_', ones(numel(PMF_pred_),1)', edges_z);

    % scoring rule - DKL of the probability map
    PMF_larger_thres_true = double(z_val > z_thresh);
    PMF_smaller_thres_true = double(z_val <= z_thresh);
    perf_val.DKL_score_largersmaller_thres = zeros(numel(z_val),1);
    for target_ = 1 : numel(z_val)
        if PMF_smaller_thres_true(target_) == 1
            perf_val.DKL_score_largersmaller_thres(target_,1) = (log2(PMF_smaller_thres_true(target_)) - log2(1-prob_thres_pred_(target_)))*PMF_smaller_thres_true(target_);
        else
            perf_val.DKL_score_largersmaller_thres(target_,1) = (log2(PMF_larger_thres_true(target_)) - log2(prob_thres_pred_(target_)))*PMF_larger_thres_true(target_);
        end
    end
    perf_val.DKL_score_largersmaller_thres_mean = mean(perf_val.DKL_score_largersmaller_thres);
    perf_val.DKL_score_larger_thres = zeros(numel(z_val),1);
    for target_ = 1 : numel(z_val)
        if PMF_larger_thres_true(target_) == 0
            perf_val.DKL_score_larger_thres(target_,1) = NaN;
        else
            perf_val.DKL_score_larger_thres(target_,1) = (log2(PMF_larger_thres_true(target_)) - log2(prob_thres_pred_(target_)))*PMF_larger_thres_true(target_);
        end
    end
    perf_val.DKL_score_larger_thres_mean = nanmean(perf_val.DKL_score_larger_thres);
    % [perf_HER_val.error_sign_thres, perf_HER_val.RMSE_thres, perf_HER_val.ME_thres, perf_HER_val.MAE_thres, perf_HER_val.NSE_thres,~,~] = f_performance_det(prob_thres_pred_, PMF_larger_thres_true);

    % misclassification (Goovaerts, 1997, p.365)
    perf_val.prob_thres = round([0:0.01:1]',2);

    for i = 1:length(perf_val.prob_thres)
        contamin_pred_ = double(prob_thres_pred_ > perf_val.prob_thres(i) ); %assigning as contaminated all places where the local risk of contamination exceeds the average risk of contamination over the region
        [perf_val.TP(i), perf_val.TN(i), perf_val.FP(i), perf_val.FN(i), perf_val.TPR(i), perf_val.TNR(i), perf_val.FPR(i), perf_val.Accuracy(i), perf_val.F1_score(i), perf_val.MCC(i), perf_val.DistanceROC(i)] = f_ROC(PMF_larger_thres_true, contamin_pred_);
        perf_val.misclassication(i,1) = 100 * (1 - perf_val.Accuracy(i)); %100 * ((sum(double(contamin_pred_ ~= PMF_thres_true_val))) ./ numel(PMF_pred_));
        perf_val.contamination(i,1) = 100 * (sum(contamin_pred_) ./ numel(z_val));
    end
    perf_val.misclass_mean = 100*(sum(double(z_val > z_thresh) ~= double(z_pred_> z_thresh)))/numel(z_val);
%     perf_val.misclass_median = 100*(sum(double(z_val > z_thresh) ~= (double(z_median_pred'> z_thresh))))/numel(z_val);

    perf_val.misclass_marginal = perf_val.misclassication(perf_val.index_marg);
    perf_val.Accuracy_marginal = 100*perf_val.Accuracy(perf_val.index_marg);
    perf_val.F1_score_marginal = perf_val.F1_score(perf_val.index_marg);
    perf_val.MCC_marginal = perf_val.MCC(perf_val.index_marg);
    perf_val.DistanceROC_marginal = perf_val.DistanceROC(perf_val.index_marg);

end