function moran_test_for_gradient()

% Test gradient correlation between GG-GW
% Load gradient data
cd('/data/results/gradients');

gradient_GG = load('gradient_GG_FC_HCP_7T_test_group.txt');
gradient_GW = load('gradient_GW_FC_HCP_7T_test_group.txt');

% Calculate Spearman correlation
[r,p] = corr(gradient_GG,gradient_GW,'type','Spearman');

% Plot scatter
figure(1);
plot(gradient_GG,gradient_GW,'*');

% Moran test
% Load W matrix
W = load('distance_ROI_center.txt');
W = W(201:600,201:600);

MEM = compute_mem(W);

% Generate null distribution
surrogate_gradient = moran_randomization(gradient_GG,MEM,5000,'procedure','singleton');
surrogate_gradient = squeeze(surrogate_gradient);

% Calculate null distribution correlations
r_moran_null = corr(surrogate_gradient,gradient_GW);

% One-tailed p-value
P_G = length(find(r_moran_null>r))/5000;

% Plot null distribution histogram
figure(2);
hist(r_moran_null);

% Save null distribution
writematrix(r_moran_null, 'null_distribution_of_gradients.txt');

end