function functional_gradient_GW_FC()

% Calculate gradients for GM-WM connectivity.

% Load FC matrix
cd('/data/results/_FC');
FC = load('FC_HCP_7T_test_group.txt');
FC = FC(201:600, 1:200);  % Extract GM-WM connections

% Calculate gradients
G = GradientMaps('kernel', 'cs', 'approach', 'dm', 'n_components', 1);
G_all = G.fit(FC, 'sparsity', 90);

G1 = G_all.gradients{1}(:,1);

% Save gradients
cd('/data/results/gradients');
writematrix(G1, 'gradient_GW_FC_HCP_7T_test_group.txt', 'Delimiter', 'space');

end