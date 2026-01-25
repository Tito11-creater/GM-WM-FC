function dispersion_of_gradient()
% Calculate dispersion metrics for gradient data across networks.

clear all;
clc;

%% Load gradient data and define networks
cd('/data/results/gradients');
gradient_GG = load('gradient_GG_FC_HCP_7T_test_group.txt');
gradient_GW = load('gradient_GW_FC_HCP_7T_test_group.txt');

% Normalize gradients to [0,1] range
gradient_GG = (gradient_GG - min(gradient_GG)) ./ (max(gradient_GG) - min(gradient_GG));
gradient_GW = (gradient_GW - min(gradient_GW)) ./ (max(gradient_GW) - min(gradient_GW));

% Network boundaries: VN, SMN, DAN, VAN, VIN, FPN, DMN
GM = [1,61,62,138,139,184,185,231,232,257,258,309,310,400];
num_networks = 7;

%% Calculate centroids for each network
center_GG = zeros(1, num_networks);
center_GW = zeros(1, num_networks);

for i = 1:num_networks
    start_idx = GM(2*i-1);
    end_idx = GM(2*i);
    center_GG(1,i) = mean(gradient_GG(start_idx:end_idx));
    center_GW(1,i) = mean(gradient_GW(start_idx:end_idx));
end

%% Calculate distances from each ROI to its network centroid
% Initialize storage for distances
GG_distances = cell(num_networks, 1);
GW_distances = cell(num_networks, 1);

network_names = {'VN', 'SMN', 'DAN', 'VAN', 'VIN', 'FPN', 'DMN'};

for net = 1:num_networks
    start_idx = GM(2*net-1);
    end_idx = GM(2*net);
    roi_count = end_idx - start_idx + 1;
    
    GG_distances{net} = zeros(1, roi_count);
    GW_distances{net} = zeros(1, roi_count);
    
    for idx = 1:roi_count
        roi_global_idx = start_idx + idx - 1;
        GG_distances{net}(idx) = abs(gradient_GG(roi_global_idx) - center_GG(net));
        GW_distances{net}(idx) = abs(gradient_GW(roi_global_idx) - center_GW(net));
    end
end

%% paired-sample t-tests between GG and GW distances for each network
t_stats = zeros(num_networks, 1);
p_values = zeros(num_networks, 1);

for net = 1:num_networks
    [~, p, ~, stats] = ttest(GG_distances{net}, GW_distances{net});
    p_values(net) = p;
    t_stats(net) = stats.tstat;
end

% FDR correction
fdr_p = mafdr(p_values, 'BHFDR', true);
writematrix(fdr_p, 'FDR_p_network_dispersion.txt');

%% Save network-specific dispersion distances
cd('/data/results/gradients');

for net = 1:num_networks
    writematrix(GG_distances{net}, ['GG_dispersion_', network_names{net}, '.txt']);
    writematrix(GW_distances{net}, ['GW_dispersion_', network_names{net}, '.txt']);
end

%% Global dispersion (distance to global centroid)
center_global_GG = mean(gradient_GG);
center_global_GW = mean(gradient_GW);

Dis_GG = abs(gradient_GG - center_global_GG);
Dis_GW = abs(gradient_GW - center_global_GW);

[~, p_g, ~, t_g] = ttest(Dis_GG, Dis_GW);
t_g = t_g.tstat;

writematrix(Dis_GG, 'GG_dispersion_global.txt');
writematrix(Dis_GW, 'GW_dispersion_global.txt');

%% Between-network dispersion
Dis_GG_between = zeros(num_networks, num_networks);
Dis_GW_between = zeros(num_networks, num_networks);

for i = 1:num_networks
    for j = 1:num_networks
        Dis_GG_between(i,j) = abs(center_GG(i) - center_GG(j));
        Dis_GW_between(i,j) = abs(center_GW(i) - center_GW(j));
    end
end

% Extract unique between-network distances (upper triangle without diagonal)
Dis_GG_b = [];
Dis_GW_b = [];

for i = 2:num_networks
    for j = 1:i-1
        Dis_GG_b(end+1) = Dis_GG_between(i,j);
        Dis_GW_b(end+1) = Dis_GW_between(i,j);
    end
end

[~, p_b, ~, t_b] = ttest(Dis_GG_b, Dis_GW_b);
t_b = t_b.tstat;

writematrix(Dis_GG_b, 'GG_dispersion_between_networks.txt');
writematrix(Dis_GW_b, 'GW_dispersion_between_networks.txt');

end