function FC_weighted_degree_calculation()

% Load group averaged FC

cd('/data/results/FC/');

FC = load('FC_HCP_7T_test_group.txt');

% Extract GM-WM FC

FC_GW = FC(201:600,1:200);

DC = sum(FC_GW,2);      % Weighted degree calculation

writematrix(DC,'DC_HCP_7T_test_group.txt','delimiter','space');

end