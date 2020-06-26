function [G, ksi_fraction, edges_interv, interv_size] = f_G(pmf_pred_nn, edges_z, z_true)
%% function to calculate the Goodness statistics (Deutsch, 1997, and Goovaerts, 2001)
% (Goovaerts, 2001 & Meirvenne and Goovaerts, 2001)

% -------------- Input --------------
% - pmf_pred_nn         {1,T}   predicted z PMF for targets 
% - edges_z             [1,n]   bin edges of the z PMF

% -------------- Output --------------
% - z_lq                [1,T]   lower quartile of the z_PMFs
% - z_median            [1,T]   median of the z_PMFs   
% - z_up                [1,T]   upper quartile of the z_PMFs

% -------------- Version --------------
% - 2020/06/24 Stephanie Thiesen: intial version

% -------------- Script --------------
    p = 0:0.01:0.5;

    for i = 1:length(p) 
        edges_interv(i,1)  = 0+p(i); %symmetric intervals
        edges_interv(i,2)  = 1-p(i); %symmetric intervals
        interv_size(i) = edges_interv(i,2)-edges_interv(i,1);
        for target_ = 1:numel(z_true)            
            cmf_ = round([0 cumsum(pmf_pred_nn{target_})],2);
            idx_  = find(cmf_ >= edges_interv(i,1), 1); %1st edge of the interval
            if idx_ == 1
                z_1s_ = edges_z(idx_);      
            else
                z_1s_ = edges_z(idx_-1); %left edge
            end
            
            idx_  = find(cmf_ >= edges_interv(i,2), 1); %2nd edge of the interval
%             if idx_ == numel(edges_z)
%                 z_2s_ = edges_z(idx_);      
%             else
%                 z_2s_ = edges_z(idx_); %right edge 
%             end
            z_2s_ = edges_z(idx_); %right edge
            
            if z_true(target_) > z_1s_ && z_true(target_) <= z_2s_
                ksi_(target_) = 1;
            else
                ksi_(target_) = 0;
            end         
        end
        ksi_fraction(i) = sum(ksi_) / numel(z_true);
        if ksi_fraction(i) >= interv_size(i)
            w_(i) = 1;
        else 
            w_(i) = -2;
        end            
    end
    for k = 1:length(w_)
        a_ = sum( w_(k) * abs(ksi_fraction(k) - interv_size(k)) ) / length(k);
    end
        G = 1 - a_;
end

