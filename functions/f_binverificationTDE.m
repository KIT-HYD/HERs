%calculate the Cross entropy of the full dataset and binned
% Thiesen, S., Darscheid, P. and Ehret, U.: Identifying rainfall-runoff 
% events in discharge time series: A data-driven method based on Information 
% Theory, Hydrol. Earth Syst. Sci., 23(2), 1015–1034, 
% doi:https://doi.org/10.5194/hess-23-1015-2019, 2019.

% for a sample, create a infogram and for the 1st distance class (lag) take
% the pair points and create a histogram. Validate the binwidth of this
% histogram

function [HPQ_x, sample_sizes] = f_binverificationTDE(zm,xm, ym, binwidth, lg, class, step_)
    % binwidth = delta z binwidth
    %lg = class size (lag width in meters if the input x,y is given in meters)
    %class = analyzed distance class
    num_rep = 500; % number of sampling repetitions 
%     step_ = 5;% step for the sample size 
    % calculate the euclidean distance and the difference between all point pairs
    mat_euc_distance_xy = NaN(length(xm),length(xm)); %matrix for the distance between 2 pair points
    mat_diff_z = NaN(length(zm),length(zm)); %matrix for the diff_z 
    for i = 1:length(xm) %repeat to each east coordinates (x)
        for j = 1:length(xm) %repeat to each north coordinates (y)
            mat_euc_distance_xy(i,j) = f_euclidean_dist(xm(i), ym(i), xm(j), ym(j)); %calculate the euclidean distance
            mat_diff_z(i,j) = f_diff(zm(i), zm(j)); %calculate the diff_z
        end
    end

    %define the distance classes (lag)
    hmax = max(mat_euc_distance_xy(:)); 
    nlag = ceil(hmax/lg); %total number of classes
    distance_classes = [0:lg:lg*nlag]; %edges of the distance classes
    bin_centers_distance_classes = distance_classes(1:end-1) + lg/2; %bin centers of the distance classes

    % compute bin edges
    max_diff_z = max(mat_diff_z(:));
    numbins = 2*floor(max_diff_z/binwidth + 1); %symmetric
    mini_ = -numbins/2*binwidth;
    maxi_ =  numbins/2*binwidth;
    edges_diff_z = linspace(mini_,maxi_,numbins+1); %edges of the diff_z pmf

    % compute delta_z PMF
    countpmf_diff_z_all_obs_plus1 = nan(1,length(edges_diff_z)-1); 
    pmf_diff_z_all_obs = nan(1,length(edges_diff_z)-1);

    %delta_z PMF for the full dataset
    idx_ = mat_euc_distance_xy > 0; %ignore the index of the observation in relation to itself

    % step 6: Calculathe the diff_z PMF by distance class
    %put diff_z in the distance classes
    obs_diff_z_by_class = cell(1, length(distance_classes)-1); %observations of diff_z by lag class (each cell is one lag class) 
    for i = 1 : length(distance_classes)-1 %for each class
        idx_ = mat_euc_distance_xy > distance_classes(i) & mat_euc_distance_xy <= distance_classes(i+1); %find the diff_z within the current lag class
        obs_diff_z_by_class{i} = mat_diff_z(idx_); %save the diff_z 
    end

    % sample_sizes to be tested
    data = obs_diff_z_by_class{class}; %(N,1)
    sample_sizes = [1:step_:size(data)]';%[50;100;150;200;250;300;350;400];
    samplingstrategy = 'continuous'; % sampling strategy

    % evaluate no predictor series  
    edges = cell(1,1); %(1,N)
    edges{1} = edges_diff_z;
    [data_binned_x, data_histcounts_x] = f_histcounts_anyd(data, edges);
    [H_x, DKL_x, HPQ_x, ~, ~, ~] = f_infomeasures_from_samples(data, edges, data_binned_x, data_histcounts_x, sample_sizes, num_rep, samplingstrategy);

end
