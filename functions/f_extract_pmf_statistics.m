function [z_entropy, z_mean, z_mode, varargout] = f_extract_pmf_statistics(pmf_pred_nn, edges_z, bin_centers_edges_z, varargin)
%% function to extract PMF entropy, mean, median, mode, exceed probability

% -------------- Input --------------
% - pmf_pred_nn         {1,T}   predicted z PMF for targets 
% - bin_centers_edges_z [1,n]   bin centers of the z PMF
% - edges_z             [1,n]   edges of the z PMF
% - varargin            h       threshold of z

% -------------- Output --------------
% - z_entropy           [T,1]   entropy of the z_PMFs
% - z_mean              [T,1]   mean of the z_PMFs
% - z_mode              [T,1]   mode of the z_PMFs
% - varargout           [T,1]   probability of z_PMFs exceeding varargin (threshold of z)      

% -------------- Version --------------
% - 2020/03/20 Stephanie Thiesen: intial version
% - 2020/04/10 Stephanie Thiesen: x_target and y_target removed

% -------------- Script --------------
    
    z_entropy = NaN(1,numel(pmf_pred_nn)); 
    z_mean = NaN(1,numel(pmf_pred_nn));
    z_mode = NaN(1,numel(pmf_pred_nn));
    z_prob = NaN(1,numel(pmf_pred_nn));
    
    
    for target_ = 1:numel(pmf_pred_nn)
        %calculate the entropy of the z PMFs for all predicted points 
        z_entropy(1, target_) = f_entropy(cell2mat(pmf_pred_nn(1,target_))); %calculate the entropy
        
        %mean
        z_mean(1,target_) = (sum(bin_centers_edges_z .* (pmf_pred_nn{1,target_}))) / sum(pmf_pred_nn{1,target_}); %interpolated mean (not the bin center)
       
        %mode
        [~,idx_] = max(cell2mat(pmf_pred_nn(1,target_)));
        z_mode(1,target_) = bin_centers_edges_z(idx_); %bin center
    end
    
    if ~isnan(varargin{1})
        if length(varargin) >= 1 %nargin >= 4
            thres = varargin{1};
            for target_ = 1:numel(pmf_pred_nn) 
                %probability of exceeding thres
                bin_thres_ = sum(edges_z <= thres);
                cmf_probab = cumsum(pmf_pred_nn{1,target_});             
                z_prob(1,target_) = 1 - cmf_probab(bin_thres_);
%                 cpmf_prob_ = cumsum(pmf_pred_nn{1,target_});
%                 bin_thres_ = find(bin_centers_edges_z <= 0.8, 1, 'last');
%                 x_prob_ = [bin_centers_edges_z(bin_thres_-1) bin_centers_edges_z(bin_thres_) bin_centers_edges_z(bin_thres_+1)];
%                 y_prob_ = [cpmf_prob_(bin_thres_-1) cpmf_prob_(bin_thres_) cpmf_prob_(bin_thres_+1)];
%                 cpmf_thresh_ = interp1(x_prob_,y_prob_,0.8); %interpolated median (not the bin center)
%                 z_target_probability_pred(1,target_) = 1 - cpmf_thresh_;
            end
        end
    else
        z_prob = NaN;
    end

    if nargout >= 2
        varargout{1} = z_prob;
    end

end