%% Extract PMF statistics, calculate performance and plot
clear
addpath('functions');
load(filename);

%% HER3 Extract PMF statistics and pot
% extract mean, mode, probability and plot

% GRID
[z_entropy_pred_, z_mean_pred_, z_mode_pred_, z_prob_pred_] = f_extract_pmf_statistics(pmf_pred_grid.pmf_pred, her.edges_z, her.bin_centers_edges_z, her.z_thresh);
    pmf_pred_grid.z_entropy_pred = reshape(z_entropy_pred_, size(pmf_pred_grid.x_target,1), size(pmf_pred_grid.y_target,2));
    pmf_pred_grid.z_mean_pred = reshape(z_mean_pred_, size(pmf_pred_grid.x_target,1), size(pmf_pred_grid.y_target,2));
    pmf_pred_grid.z_mode_pred = reshape(z_mode_pred_, size(pmf_pred_grid.x_target,1), size(pmf_pred_grid.y_target,2));
    pmf_pred_grid.z_prob_pred = reshape(z_prob_pred_, size(pmf_pred_grid.x_target,1), size(pmf_pred_grid.y_target,2));
x_plot = [3.6, 4.0]; y_plot = [4.4, 4.7]; %coordinates to plot predicted PMF
f_plot_prediction(pmf_pred_grid.z_mean_pred, pmf_pred_grid.z_entropy_pred, pmf_pred_grid.pmf_pred, pmf_pred_grid.x_target, pmf_pred_grid.y_target, x, y, z, idx_cal, idx_val, her, x_plot, y_plot, shp_basin);
f_plot_probabilitymap(pmf_pred_grid.z_prob_pred, her.z_thresh, txt, pmf_pred_grid.x_target, pmf_pred_grid.y_target, x, y, z, idx_cal, idx_val, shp_basin);

clear z_entropy_pred_ z_mean_pred_ z_median_pred_ z_mode_pred_ z_prob_pred_x_plot y_plot
clear h1 i idx_ ncols nrows randplot_set_ z_predicted_ target_ 

%% HER3: Calculate performance metrics (Ralf, this is not working for my dataset, but it will work for yours!)
% Define
z_pred_ = round(pmf_pred_grid.z_mean_pred(idx_val),5)';
PMF_pred_ = pmf_pred_grid.pmf_pred(idx_val)';
prob_thres_pred_ = round(pmf_pred_grid.z_prob_pred(idx_val),5)';%cpmf>threshold
z_true_ = z(idx_val);
perf_HER_val.marginal_prob_thres_cal = sum( double(z(idx_cal)) > her.z_thresh) / length(idx_cal); %marginal probability of contamination over the region

% Root mean square deviation RMSD (RMSE), Mean Error (ME), Mean Absolute Error (MAE) L1 norm, robust parameter estimator
% Nash-Sutcliffe model efficiency (r2, coefficient of determination)
[perf_HER_val.error_sign, perf_HER_val.RMSE, perf_HER_val.ME, perf_HER_val.MAE, perf_HER_val.NSE, perf_HER_val.correl_true_pred, perf_HER_val.correl_true_residue] = f_performance_det(z_pred_,z_true_);

% scoring rule - DKL in relation to the mean value | bin_true optimization HER03
perf_HER_val.DKL_score_mean = f_performance_prob(z_true_', PMF_pred_', ones(numel(PMF_pred_),1)', her.edges_z);

% scoring rule - DKL of the probability map | threshold optimization HER04
obs_larger_thres_ = double(z_true_ > her.z_thresh);
obs_smaller_thres_ = double(z_true_ <= her.z_thresh);
perf_HER_val.DKL_score_01_largersmaller_thres = zeros(numel(PMF_pred_),1); 
for target_ = 1 : numel(PMF_pred_)
    if obs_smaller_thres_(target_) == 1
        perf_HER_val.DKL_score_01_largersmaller_thres(target_,1) = (log2(obs_smaller_thres_(target_)) - log2(1-prob_thres_pred_(target_)))*obs_smaller_thres_(target_);
    else
        perf_HER_val.DKL_score_01_largersmaller_thres(target_,1) = (log2(obs_larger_thres_(target_)) - log2(prob_thres_pred_(target_)))*obs_larger_thres_(target_);
    end
end
perf_HER_val.DKL_score_01_largersmaller_thres_mean = mean(perf_HER_val.DKL_score_01_largersmaller_thres);
perf_HER_val.DKL_score_1_larger_thres = zeros(numel(PMF_pred_),1);
for target_ = 1 : numel(PMF_pred_)
    if obs_larger_thres_(target_) == 0
        perf_HER_val.DKL_score_1_larger_thres(target_,1) = NaN;
    else
        perf_HER_val.DKL_score_1_larger_thres(target_,1) = (log2(obs_larger_thres_(target_)) - log2(prob_thres_pred_(target_)))*obs_larger_thres_(target_);
    end
end
perf_HER_val.DKL_score_1_larger_thres_mean = nanmean(perf_HER_val.DKL_score_1_larger_thres);

% misclassification (Goovaerts, 1997, p.365)
perf_HER_val.prob_thres = round([0:0.01:1]',2);
perf_HER_val.index_marg = find(perf_HER_val.prob_thres == round(perf_HER_val.marginal_prob_thres_cal,2));

for i = 1:length(perf_HER_val.prob_thres)
    contamin_pred_ = double(prob_thres_pred_ > perf_HER_val.prob_thres(i) ); %assigning as contaminated all places where the local risk of contamination exceeds the average risk of contamination over the region
    [perf_HER_val.TP(i), perf_HER_val.TN(i), perf_HER_val.FP(i), perf_HER_val.FN(i), perf_HER_val.TPR(i), perf_HER_val.TNR(i), perf_HER_val.FPR(i), perf_HER_val.Accuracy(i), perf_HER_val.F1_score(i), perf_HER_val.MCC(i), perf_HER_val.DistanceROC(i)] = f_ROC(obs_larger_thres_, contamin_pred_);
    perf_HER_val.misclassication(i,1) = 100 * (1 - perf_HER_val.Accuracy(i)); %100 * ((sum(double(contamin_pred_ ~= PMF_thres_true_val))) ./ numel(PMF_pred_));
    perf_HER_val.contamination(i,1) = 100 * (sum(contamin_pred_) ./ numel(PMF_pred_));
end

% fraction of true values between symmetric probability interval (Goovaerts, 2001 & Meirvenne and Goovaerts, 2001)
[perf_HER_val.G, perf_HER_val.ksi_fraction, perf_HER_val.edges_interv, perf_HER_val.interv_size] = f_G(PMF_pred_, her.edges_z, z_true_);
perf_HER_val.ksi_fraction_uq_lq =  perf_HER_val.ksi_fraction(find(perf_HER_val.edges_interv(:,1)==0.25));% fraction of true values between lower and upper quantiles (Goovaerts, 2001 & Meirvenne and Goovaerts, 2001)

perf_HER_val.misclass_mean = 100*(sum(double(z_val > her.z_thresh) ~= (double(round(pmf_pred_fullset.z_mean_pred(idx_val),5)'> her.z_thresh))))/numel(PMF_pred_);
perf_HER_val.misclass_min = min(perf_HER_val.misclassication(:));
perf_HER_val.misclass_marginal = perf_HER_val.misclassication(perf_HER_val.index_marg);
perf_HER_val.Accuracy_marginal = 100*perf_HER_val.Accuracy(perf_HER_val.index_marg);
perf_HER_val.F1_score_marginal = perf_HER_val.F1_score(perf_HER_val.index_marg);
perf_HER_val.MCC_marginal = perf_HER_val.MCC(perf_HER_val.index_marg);
perf_HER_val.DistanceROC_marginal = perf_HER_val.DistanceROC(perf_HER_val.index_marg);

clear pmf_ z_true_ contamin_pred_ i target_ idx_ z_lq_pred_ z_uq_pred_ obs_larger_thres_ obs_smaller_thres_ xlim ylim

%plots

%residue
[z_true_sorted_, idx_] = sort(z_val);
z_predicted_sorted_ = z_pred_(idx_);

figure
subplot(1,2,1) %Correlation
scatter(z_true_sorted_,z_predicted_sorted_);
ylabel('Expected Z');
xlabel('True Z');
xlim([-2 2]);
ylim([-2 2]);
subplot(1,2,2) %Residue correlation
scatter(z_true_sorted_,z_predicted_sorted_ - z_true_sorted_);
xlabel('True Z');
ylabel('Residue');
xlim([-3 3]);
ylim([-3 3]);
sgtitle('Correlation');

%ROC curve
figure;
hold on;
plot(perf_HER_val.FPR,perf_HER_val.TPR,'o-');
ylabel('TPR');
xlabel('FPR');
title('ROC curve');

figure
subplot(1,2,1);
plot(perf_HER_val.prob_thres , perf_HER_val.contamination,'o-');
hold on
ylim = 1;
xlabel('Probability p');
ylabel('Locations declared contaminated [%]');
pbaspect([1 1 1]);
title('Proportion of test locations that are declared contaminated');

subplot(1,2,2);
plot(perf_HER_val.prob_thres , perf_HER_val.misclassication,'o-');
hold on
ylim = 1;
xlabel('Probability p');
ylabel('Misclassification [%]');
line([perf_HER_val.marginal_prob_thres_cal, perf_HER_val.marginal_prob_thres_cal], get(gca, 'ylim'),'Color','black','LineStyle','--');
text(perf_HER_val.marginal_prob_thres_cal,50,'Marginal Prob. \rightarrow','HorizontalAlignment','right');
pbaspect([1 1 1]);
title('Proportion of test locations that are wrongly classified');

sgtitle({strcat('Probability threshold p / ', txt);''});

clear prob_thres_pred_ PMF_pred_ z_pred_ ylim idx_ z_predicted_sorted_ z_true_sorted_ ans xlim ylim

%% sHER4 Extract PMF statistics and pot
% extract mean, mode, probability and plot

% GRID
k = 1;
for i = 1:size(pmf_simul_grid.z_realiz,1) %extract PMF
    for j = 1:size(pmf_simul_grid.z_realiz,2)
        realiz_ = reshape(pmf_simul_grid.z_realiz(i,j,:), [1, size(pmf_simul_grid.z_realiz,3)]); %flatten matrix
%         n_pairs_by_bin_ = histcounts(realiz_,her.edges_z);
        pmf_realiz_ = histcounts(realiz_,her.edges_z,'Normalization', 'probability');
%         pmf_realiz_ = pmf_realiz_ + her.prob_min; %avoiding empty bins
        pmf_simul_grid.pmf_realiz{1,k} = pmf_realiz_ ./ sum(pmf_realiz_(:)); %normalization
        pmf_simul_grid.realiz_mean(i,j) = mean(realiz_);
        k = k + 1;
    end
end
[z_entropy_realiz_, z_mean_realiz_, z_mode_realiz_, z_prob_realiz_] = f_extract_pmf_statistics(pmf_simul_grid.pmf_realiz, her.edges_z, her.bin_centers_edges_z, her.z_thresh);
    pmf_simul_grid.z_entropy_realiz = reshape(z_entropy_realiz_, size(pmf_simul_grid.x_target,1), size(pmf_simul_grid.y_target,2));
    pmf_simul_grid.z_mean_realiz = reshape(z_mean_realiz_, size(pmf_simul_grid.x_target,1), size(pmf_simul_grid.y_target,2));
    pmf_simul_grid.z_mode_realiz = reshape(z_mode_realiz_, size(pmf_simul_grid.x_target,1), size(pmf_simul_grid.y_target,2));
    pmf_simul_grid.z_prob_realiz = reshape(z_prob_realiz_, size(pmf_simul_grid.x_target,1), size(pmf_simul_grid.y_target,2));
% plot
x_plot = [3.6, 4.0]; y_plot = [4.4, 4.7]; %coordinates to plot predicted PMF
f_plot_prediction(pmf_simul_grid.z_mean_realiz, pmf_simul_grid.z_entropy_realiz, pmf_simul_grid.pmf_realiz, pmf_simul_grid.x_target, pmf_simul_grid.y_target, x, y, z, idx_cal, idx_val, her, x_plot, y_plot, shp_basin);
f_plot_probabilitymap(pmf_simul_grid.z_prob_realiz, her.z_thresh, txt, pmf_simul_grid.x_target, pmf_simul_grid.y_target, x, y, z, idx_cal, idx_val, shp_basin);

clear z_entropy_realiz_ z_mean_realiz_ z_mode_realiz_ z_prob_realiz_
clear z_mean_realiz_ z_prob_realiz_ z_mode_realiz_ z_entropy_realiz_ realiz_ i idx_ ncols nrows randplot_set_ z_predicted_ target_ x_plot y_plot i j k h1 ncols nrows pmf_realiz_ randplot_set_ 

% sHER4 plot realization
clims_logPb = [1.2 2.5];
clims_HlogPb = [1.4 2.2];
sz1 = 60;
sz2 = 35;
symspec_ = makesymbolspec('Polygon', {'ROCK', 'basin','FaceColor', [1 1 1], 'FaceAlpha',0,'LineWidth',2, 'EdgeColor', [0 0 0]});

i = 15; %realization number
z_plot_clip = pmf_simul_grid.z_realiz(:,:,i);
in_ = inpolygon(pmf_simul_grid.x_target, pmf_simul_grid.y_target, shp_basin.X, shp_basin.Y);
z_plot_clip(~in_) = nan;

figure
hold on
pcolor(pmf_simul_grid.x_target, pmf_simul_grid.y_target, z_plot_clip);
scatter(x(idx_cal), y(idx_cal), 1+70*normalize(z(idx_cal),'range'),'k+','LineWidth',1);
set(gca,'ColorScale','log')
colormap(gca,parula(8));
title({escape(filename);sprintf('sHER realization #%i',i);''});
h=colorbar;
caxis(clims_logPb);
shading flat;
pbaspect([1 1 1]);
mapshow(shp_basin,'SymbolSpec', symspec_);

clear h h1 i idx_ in_ symspec_ sz1 sz2 z_plot_clip clims_HlogPb clims_logPb 
%% sHER4: Calculate performance metrics (Ralf, this is not working for my dataset, but it will work for yours!)
% Define
z_realiz_ = round(pmf_simul_grid.z_mean_realiz(idx_val),5)';
PMF_realiz_ = pmf_simul_grid.pmf_realiz(idx_val)';
prob_thres_realiz_ = round(pmf_simul_grid.z_prob_realiz(idx_val),5)';%cpmf>threshold
z_true_ = z(idx_val);
perf_sHER_val.marginal_prob_thres_cal = sum( double(z(idx_cal)) > her.z_thresh) / length(idx_cal); %marginal probability of contamination over the region

% Root mean square deviation RMSD (RMSE)
% Mean Error (ME)  
% Mean Absolute Error (MAE) L1 norm, robust parameter estimator
% Nash-Sutcliffe model efficiency (r2, coefficient of determination)
[perf_sHER_val.error_sign, perf_sHER_val.RMSE, perf_sHER_val.ME, perf_sHER_val.MAE, perf_sHER_val.NSE, perf_sHER_val.correl_true_pred, perf_sHER_val.correl_true_residue] = f_performance_det(z_realiz_,z_true_);

% scoring rule - DKL in relation to the mean value | bin_true optimization HER03
perf_sHER_val.DKL_score_mean = f_performance_prob(z_true_', PMF_realiz_', ones(numel(PMF_realiz_),1)', her.edges_z);

% scoring rule - DKL of the probability map | threshold optimization HER04
obs_larger_thres_ = double(z_true_ > her.z_thresh);
obs_smaller_thres_ = double(z_true_ <= her.z_thresh);
perf_sHER_val.DKL_score_01_largersmaller_thres = zeros(numel(PMF_realiz_),1); 
for target_ = 1 : numel(PMF_realiz_)
    if obs_smaller_thres_(target_) == 1
        perf_sHER_val.DKL_score_01_largersmaller_thres(target_,1) = (log2(obs_smaller_thres_(target_)) - log2(1-prob_thres_realiz_(target_)))*obs_smaller_thres_(target_);
    else
        perf_sHER_val.DKL_score_01_largersmaller_thres(target_,1) = (log2(obs_larger_thres_(target_)) - log2(prob_thres_realiz_(target_)))*obs_larger_thres_(target_);
    end
end
perf_sHER_val.DKL_score_01_largersmaller_thres_mean = mean(perf_sHER_val.DKL_score_01_largersmaller_thres);
perf_sHER_val.DKL_score_1_larger_thres = zeros(numel(PMF_realiz_),1);
for target_ = 1 : numel(PMF_realiz_)
    if obs_larger_thres_(target_) == 0
        perf_sHER_val.DKL_score_1_larger_thres(target_,1) = NaN;
    else
        perf_sHER_val.DKL_score_1_larger_thres(target_,1) = (log2(obs_larger_thres_(target_)) - log2(prob_thres_realiz_(target_)))*obs_larger_thres_(target_);
    end
end
perf_sHER_val.DKL_score_1_larger_thres_mean = nanmean(perf_sHER_val.DKL_score_1_larger_thres);

% misclassification (Goovaerts, 1997, p.365)
perf_sHER_val.prob_thres = round([0:0.01:1]',2);
perf_sHER_val.index_marg = find(perf_sHER_val.prob_thres == round(perf_sHER_val.marginal_prob_thres_cal,2));
for i = 1:length(perf_sHER_val.prob_thres)
    contamin_pred_ = double(prob_thres_realiz_ > perf_sHER_val.prob_thres(i) ); %assigning as contaminated all places where the local risk of contamination exceeds the average risk of contamination over the region
    [perf_sHER_val.TP(i), perf_sHER_val.TN(i), perf_sHER_val.FP(i), perf_sHER_val.FN(i), perf_sHER_val.TPR(i), perf_sHER_val.TNR(i), perf_sHER_val.FPR(i), perf_sHER_val.Accuracy(i), perf_sHER_val.F1_score(i), perf_sHER_val.MCC(i), perf_sHER_val.DistanceROC(i)] = f_ROC(obs_larger_thres_, contamin_pred_);
    perf_sHER_val.misclassication(i,1) = 100 * (1 - perf_sHER_val.Accuracy(i)); %100 * ((sum(double(contamin_pred_ ~= PMF_thres_true_val))) ./ numel(PMF_pred_));
    perf_sHER_val.contamination(i,1) = 100 * (sum(contamin_pred_) ./ numel(PMF_realiz_));
end

% fraction of true values between symmetric probability interval (Goovaerts, 2001 & Meirvenne and Goovaerts, 2001)
[perf_sHER_val.G, perf_sHER_val.ksi_fraction, perf_sHER_val.edges_interv, perf_sHER_val.interv_size] = f_G(PMF_realiz_, her.edges_z, z_true_);
perf_sHER_val.ksi_fraction_uq_lq =  perf_sHER_val.ksi_fraction(find(perf_sHER_val.edges_interv(:,1)==0.25));% fraction of true values between lower and upper quantiles (Goovaerts, 2001 & Meirvenne and Goovaerts, 2001)

perf_sHER_val.misclass_mean = 100*(sum(double(z_val > her.z_thresh) ~= (double(round(pmf_simul_fullset.z_mean_realiz(idx_val),5)'> her.z_thresh))))/numel(PMF_realiz_);
perf_sHER_val.misclass_min = min(perf_sHER_val.misclassication(:));
% perf_HER_val.misclass_median = 100*(sum(double(z_val > her.z_thresh) ~= (double(round(pmf_pred_fullset.z_median_pred(idx_val),5)'> her.z_thresh))))/numel(PMF_pred_);
perf_sHER_val.misclass_marginal = perf_sHER_val.misclassication(perf_sHER_val.index_marg);
perf_sHER_val.Accuracy_marginal = 100*perf_sHER_val.Accuracy(perf_sHER_val.index_marg);
perf_sHER_val.F1_score_marginal = perf_sHER_val.F1_score(perf_sHER_val.index_marg);
perf_sHER_val.MCC_marginal = perf_sHER_val.MCC(perf_sHER_val.index_marg);
perf_sHER_val.DistanceROC_marginal = perf_sHER_val.DistanceROC(perf_sHER_val.index_marg);

clear pmf_ z_true_ contamin_pred_ i target_ idx_ z_lq_pred_ z_uq_pred_ obs_larger_thres_ obs_smaller_thres_ xlim ylim

% plots

%residue
[z_true_sorted_, idx_] = sort(z_val);
z_predicted_sorted_ = z_realiz_(idx_);

figure
subplot(1,2,1) %Correlation
scatter(z_true_sorted_,z_predicted_sorted_);
ylabel('Expected Z');
xlabel('True Z');
xlim([-2 2]);
ylim([-2 2]);
subplot(1,2,2) %Residue correlation
scatter(z_true_sorted_,z_predicted_sorted_ - z_true_sorted_);
xlabel('True Z');
ylabel('Residue');
xlim([-3 3]);
ylim([-3 3]);
sgtitle('Correlation');

%ROC curve
figure;
hold on;
plot(perf_sHER_val.FPR,perf_sHER_val.TPR,'o-');
ylabel('TPR');
xlabel('FPR');
title('ROC curve');

figure
subplot(1,2,1);
plot(perf_sHER_val.prob_thres , perf_sHER_val.contamination,'o-');
hold on
ylim = 1;
xlabel('Probability p');
ylabel('Locations declared contaminated [%]');
pbaspect([1 1 1]);
title('Proportion of test locations that are declared contaminated');

subplot(1,2,2);
plot(perf_sHER_val.prob_thres , perf_sHER_val.misclassication,'o-');
hold on
ylim = 1;
xlabel('Probability p');
ylabel('Misclassification [%]');
line([perf_sHER_val.marginal_prob_thres_cal, perf_sHER_val.marginal_prob_thres_cal], get(gca, 'ylim'),'Color','black','LineStyle','--');
text(perf_sHER_val.marginal_prob_thres_cal,50,'Marginal Prob. \rightarrow','HorizontalAlignment','right');
pbaspect([1 1 1]);
title('Proportion of test locations that are wrongly classified');

sgtitle({strcat('Probability threshold p / ', txt);''});

clear contamin_pred_ z_realiz_ z_uq_pred_ z_true_sorted_ z_true_ z_predicted_sorted_ ...
    z_lq_pred_ prob_thres_realiz_ PMF_larger_thres_true_ PMF_smaller_thres_true_ PMF_realiz_
clear h1 h i idx_ in_ clims symspec_ sz1 sz2 target_ z_plot_clip ylim

%% save
save(filename);

%% ERGODIC FLUCTUATIONS clip (Just for simulation, not sure if you are interested for now)
%     in_ = inpolygon(pmf_pred_grid.x_target, pmf_pred_grid.y_target, shp_basin.X, shp_basin.Y);
%     x_plot_clip = pmf_pred_grid.x_target;
%     x_plot_clip(~in_) = nan;
%     x_plot_clip = x_plot_clip(:);
%     x_plot_clip = x_plot_clip(~isnan(x_plot_clip));
%     y_plot_clip = pmf_pred_grid.y_target;
%     y_plot_clip(~in_) = nan;
%     y_plot_clip = y_plot_clip(:);
%     y_plot_clip = y_plot_clip(~isnan(y_plot_clip));  
% 
% % Plot CCDF
% figure
% hold on
%     % callibration
% h1 = cdfplot(z_cal);
% set( h1, 'LineStyle', '-', 'Color', 'r','LineWidth',2);
%     % E-type
% z_plot_clip = pmf_simul_grid.z_mean_realiz(:); %pmf_simul_grid.realiz_mean(:);
% z_plot_clip(~in_) = nan;
% z_plot_clip = z_plot_clip(:);
% z_plot_clip = z_plot_clip(~isnan(z_plot_clip));
% h2 = cdfplot(z_plot_clip);
% set( h2, 'LineStyle', '-', 'Color', 'b','LineWidth',.5); 
%     % Simulations
% for i = 1:size(pmf_simul_grid.z_realiz,3)
%     z_plot_clip = pmf_simul_grid.z_realiz(:,:,i);
%     z_plot_clip(~in_) = nan;
%     z_plot_clip = z_plot_clip(:);
%     z_plot_clip = z_plot_clip(~isnan(z_plot_clip));
%     h3 = cdfplot(z_plot_clip(:));
%     set( h3, 'LineStyle', '-', 'Color', [0.6 0.6 0.6] ,'LineWidth',0.1);
% end
% h1 = cdfplot(z_cal);
% set( h1, 'LineStyle', '-', 'Color', 'r','LineWidth',2);
% z_plot_clip = pmf_simul_grid.z_mean_realiz(:);
% z_plot_clip(~in_) = nan;
% z_plot_clip = z_plot_clip(:);
% z_plot_clip = z_plot_clip(~isnan(z_plot_clip));
% h2 = cdfplot(z_plot_clip);
% set( h2, 'LineStyle', '-', 'Color', 'b','LineWidth',.5); 
% % ylim([0 1.1]);
% 
% legend({'Calibration set', 'E-type simulation', 'sHER realizations'});
% xlabel('z concentration');
% ylabel('Probability');
% pbaspect([1.5 1 1]);
% 
% % Infogram
% figure;
% hold on;
%     %calibration set
% plot(her.bin_centers_distance_classes, her.H_diff_z_by_class,'Marker', '.', 'Color', 'red', 'MarkerSize', 20); 
%     %E-type
% z_plot_clip = pmf_simul_grid.z_mean_realiz(:);	%simulated PMF mean
%     z_plot_clip(~in_) = nan;
%     z_plot_clip = z_plot_clip(:);
%     z_plot_clip = z_plot_clip(~isnan(z_plot_clip));
% her_ = f_infogram_realizplot(x_plot_clip, y_plot_clip, z_plot_clip, her);
% plot(her_.bin_centers_distance_classes, her_.H_diff_z_by_class,'Marker', '.', 'Color', 'blue', 'MarkerSize', 10); 
% % z_plot_clip = pmf_simul_grid.realiz_mean(:);	%simulated mean of realizations
% %     z_plot_clip(~in_) = nan;
% %     z_plot_clip = z_plot_clip(:);
% %     z_plot_clip = z_plot_clip(~isnan(z_plot_clip));
% % her_ = f_infogram_realizplot(x_plot_clip, y_plot_clip, z_plot_clip, her);
% % plot(her_.bin_centers_distance_classes, her_.H_diff_z_by_class,'Marker', '.', 'Color', 'green', 'MarkerSize', 10); 
% % z_plot_clip = pmf_pred_grid.z_mean_pred(:);     %predicted mean
% %     z_plot_clip(~in_) = nan;
% %     z_plot_clip = z_plot_clip(:);
% %     z_plot_clip = z_plot_clip(~isnan(z_plot_clip));
% % her_ = f_infogram_realizplot(x_plot_clip, y_plot_clip, z_plot_clip, her);
% % plot(her_.bin_centers_distance_classes, her_.H_diff_z_by_class,'Marker', '.', 'Color', 'magenta', 'MarkerSize', 10); 
% 
%     % realizations
% for i = 1:size(pmf_simul_grid.z_realiz,3)
%     z_plot_clip = pmf_simul_grid.z_realiz(:,:,i);
%     z_plot_clip(~in_) = nan;
%     z_plot_clip = z_plot_clip(:);
%     z_plot_clip = z_plot_clip(~isnan(z_plot_clip));
%     her_ = f_infogram_realizplot(x_plot_clip, y_plot_clip, z_plot_clip, her);
%     h2 = plot(her_.bin_centers_distance_classes, her_.H_diff_z_by_class,'Marker', '.', 'MarkerSize', 1);
%     set( h2, 'LineStyle', '-', 'Color', [0.6 0.6 0.6] ,'LineWidth',0.1);
%     i
% end
% %calibration set
% plot(her.bin_centers_distance_classes, her.H_diff_z_by_class,'Marker', '.', 'Color', 'red', 'MarkerSize', 20); 
% z_plot_clip = pmf_simul_grid.z_mean_realiz(:);
%     z_plot_clip(~in_) = nan;
%     z_plot_clip = z_plot_clip(:);
%     z_plot_clip = z_plot_clip(~isnan(z_plot_clip));
% her_ = f_infogram_realizplot(x_plot_clip, y_plot_clip, z_plot_clip, her);
% plot(her_.bin_centers_distance_classes, her_.H_diff_z_by_class,'Marker', '.', 'Color', 'blue', 'MarkerSize', 10); 
% 
% title(escape({strcat('Infogram / ', her.txt);''}));
% pbaspect([1.5 1 1]);
% %     line([0 hmax],[her.H_diff_z her.H_diff_z], 'LineWidth', 2);
% for i = 2 : length(her.edges_distance_classes)
%     line([her.edges_distance_classes(i),her.edges_distance_classes(i)], get(gca, 'ylim'),'Color','black','LineStyle','--');
% end
% legend({'Calibration set', 'E-type simulation','sHER realizations'},'Location','southwest');
% xlabel('Euclidean distance');
% ylabel('Entropy [bit]');
% 
% clear her_ h1 h2 h3 i in_ x_plot_clip y_plot_clip z_plot_clip clims_logPb clims_HlogPb


