function FC_weighted_degree_for_WM_networks()

% Calculate degree centrality for WM networks from GW-regressed FC matrix
cd('/data/results/FC');

% Load GW-regressed FC matrix
FC_GW = load('FC_GW_regress_HCP_7T_test_group.txt');

% WM network node indices:
% 1-24: WM DMN
% 25-57: WM FPN
% 58-121: WM SMN
% 122-146: WM AN
% 147-173: WM VN
% 174-187: WM BSN
% 188-200: WM CBN

% Initialize degree centrality matrix (400×7)
DC_net = zeros(400, 7);

% Calculate degree centrality for each WM network
DC_net(:, 1) = sum(FC_GW(:, 1:24), 2);   % DMN
DC_net(:, 2) = sum(FC_GW(:, 25:57), 2);  % FPN
DC_net(:, 3) = sum(FC_GW(:, 58:121), 2); % SMN
DC_net(:, 4) = sum(FC_GW(:, 122:146), 2);% AN
DC_net(:, 5) = sum(FC_GW(:, 147:173), 2);% VN
DC_net(:, 6) = sum(FC_GW(:, 174:187), 2);% BSN
DC_net(:, 7) = sum(FC_GW(:, 188:200), 2);% CBN

% Save individual network degree centrality vectors
writematrix(DC_net(:, 1), fullfile(net_dc_path, 'DC_WM_DMN.txt'), 'delimiter', 'space');
writematrix(DC_net(:, 2), fullfile(net_dc_path, 'DC_WM_FPN.txt'), 'delimiter', 'space');
writematrix(DC_net(:, 3), fullfile(net_dc_path, 'DC_WM_SMN.txt'), 'delimiter', 'space');
writematrix(DC_net(:, 4), fullfile(net_dc_path, 'DC_WM_AN.txt'), 'delimiter', 'space');
writematrix(DC_net(:, 5), fullfile(net_dc_path, 'DC_WM_VN.txt'), 'delimiter', 'space');
writematrix(DC_net(:, 6), fullfile(net_dc_path, 'DC_WM_BSN.txt'), 'delimiter', 'space');
writematrix(DC_net(:, 7), fullfile(net_dc_path, 'DC_WM_CBN.txt'), 'delimiter', 'space');

end