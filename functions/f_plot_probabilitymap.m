function A = f_plot_probabilitymap(z_prob_plot, z_thresh, txt, x_target_grid, y_target_grid, x, y, z, idx_cal,idx_val, varargin)
%% function to plot the probabilty of z > z_thresh (HER probability map)
% -------------- Input -------------- 
% - z_prob_plot         [T,T]          GRID of probability results to be ploted
% - z_thresh             t             z threshold 
% - txt                 char           dataset name
% - x_target_grid; y_target_grid   [T,T]    x,y coordinates of the GRID
% - x; y                [n,1]          x,y coordinates of the original dataset
% - z                   [n,1]          z true values the predicted locations
% - idx_cal             [1,c]          index of the calibration set
% - varargin (shp_basin)  struc        basin shapefile

% -------------- Version --------------
% - 2020/03/20 Stephanie Thiesen: intial version
% - 2020/03/23 Stephanie Thiesen: removed possibility to plot without GRID
% - 2020/03/23 Stephanie Thiesen: clip to shp

% -------------- Script --------------

    if length(varargin) >= 1
        shp_basin = varargin{1};
    end
sz1 = 60;
sz2 = 35;
    % probabilty map
    figure;
    hold on;
    z_prob_plot_clip = z_prob_plot;
    if exist('shp_basin','var')
        in_ = inpolygon(x_target_grid, y_target_grid, shp_basin.X, shp_basin.Y);
        z_prob_plot_clip(~in_) = nan;
    end
    pcolor(x_target_grid, y_target_grid, z_prob_plot_clip);
    xlabel('x');
    ylabel('y');
    shading flat;
    colormap(flipud(gray(60)));
    % set(gca)
    h=colorbar;  
%     scatter(x(idx_cal), y(idx_cal), 1+70*normalize(z(idx_cal),'range'),'r+','LineWidth',1);
    z_cal = z(idx_cal); x_cal = x(idx_cal); y_cal = y(idx_cal);
    z_val = z(idx_val); x_val = x(idx_val); y_val = y(idx_val);
    idx_ = (z_cal > z_thresh);
    scatter(x_cal(idx_), y_cal(idx_), sz1, 'MarkerEdgeColor',[1 0 0]);
    idx_ = (z_cal <= z_thresh);
    scatter(x_cal(idx_), y_cal(idx_),sz1 ,'MarkerEdgeColor',[0 0 0]);
    idx_ = (z_val > z_thresh);
    scatter(x_val(idx_), y_val(idx_),sz2 ,'s','MarkerEdgeColor',[1 0 0]);
    idx_ = (z_val <= z_thresh);
    scatter(x_val(idx_), y_val(idx_),sz2 ,'s','MarkerEdgeColor',[0 0 0]);
    caxis([0 1]);
    title(escape({strcat('Probability map [Probability Z>',num2str(z_thresh),'] / ', txt);''}));
    pbaspect([1 1 1]);
    limits = [0 15 30 45 60 75 100];
    colormap([repmat([.9 .9 .9],limits(2)-limits(1)+1,1);repmat([0.75 0.75 0.75],limits(3)-limits(2)+1,1); repmat([0.6 0.6 0.6],limits(4)-limits(3)+1,1);repmat([0.4 0.4 0.4],limits(5)-limits(4)+1,1); repmat([0.3 0.3 0.3],limits(6)-limits(5)+1,1); repmat([0.1 0.1 0.1],limits(7)-limits(6)+1,1)]);
    if exist('shp_basin','var')
        symspec_ = makesymbolspec('Polygon', ...
           {'ROCK', 'basin','FaceColor', [1 1 1], 'FaceAlpha',0,'LineWidth',2, 'EdgeColor', [0 0 0]});
        hold on
        mapshow(shp_basin,'SymbolSpec', symspec_);
    end
    legend(strcat('Probability Z >',num2str(z_thresh)), 'cal. set > thresh.','cal. set <= thresh.', 'val. set > thresh.','val. set <= thresh.','Location','northwest');

end