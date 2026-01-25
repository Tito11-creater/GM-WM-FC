function regress_FC()

% Regression analysis between FC and weighted degree

% Load data
cd('/data/results/FC');

FC = load('FC_HCP_7T_test_group.txt');
FC_GW = FC(201:600, 1:200);
DC = load('DC_HCP_7T_test_group');  % Global weighted degree

% Regress out weighted degree from FC
FC_regress = zeros(400, 200);

for i = 1:200
    % Regress from each FC column
    b = regress(FC_GW(:, i), [ones(400, 1), DC]);
    FC_regress(:, i) = FC_GW(:, i) - sum(b' .* [ones(400, 1), DC], 2);
end

% Save regressed FC matrix
writematrix(FC_regress,'FC_GW_regress_HCP_7T_test_group.txt','Delimiter', 'space');

end
