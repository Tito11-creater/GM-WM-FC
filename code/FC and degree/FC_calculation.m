function FC_calculation()

% Calculate individual FC matrices and group-averaged FC matrices for HCP 7T test and retest groups.

%% Load network ordering labels
cd('/data/parcellation');
network_label = load('network label sort.txt');

%% Process HCP 7T test group
fprintf('Processing HCP 7T test group...\n');
filename = dir('/data/HCP_7T_test/');
sub_num = length(filename);

FC_test_group = zeros(600, 600);

for fileth = 1:sub_num
    
    % Load BOLD data
    sub_path = filename(fileth).name;
    cd(['/data/HCP_7T_test/', sub_path]);
    
    BOLD = load('BOLD_600ROI_test.txt');
    
    % Calculate and reorder FC matrix
    FC = corr(BOLD');
    FC = FC(network_label,network_label);
    
    % Save individual FC
    writematrix(FC, 'FC_test_ind.txt', 'delimiter', 'space');
    
    % Fisher Z-transformed FC for group average
    FC_z = 0.5 * log((1 + FC) ./ (1 - FC));
    FC_test_group = FC_test_group + FC_z ./ sub_num;
    
    fprintf('HCP 7T test subject: %d/%d\n', fileth, sub_num);
end

% Convert group average back to correlation
FC_test_group = (exp(2 * FC_test_group) - 1) ./ (exp(2 * FC_test_group) + 1);

for i = 1:1:600
    FC_test_group(i,i) = 0;  % Set diagonal to zero
end

% Save group result
cd('/data/results/FC/');
writematrix(FC_test_group, 'FC_HCP_7T_test_group.txt', 'delimiter', 'space');
fprintf('Processing HCP 7T test group done.\n\n');

%% Process HCP 7T retest group
fprintf('Processing HCP 7T retest group...\n');
filename = dir('/data/HCP_7T_retest/');
sub_num = length(filename);

FC_retest_group = zeros(600, 600);

for fileth = 1:sub_num
    % Load BOLD data
    sub_path = filename(fileth).name;
    cd(['/data/HCP_7T_retest/', sub_path]);
    
    BOLD = load('BOLD_600ROI_retest.txt');
    
    % Calculate and reorder FC matrix
    FC = corr(BOLD');
    FC(isnan(FC)) = 0;
    FC = FC(network_label,network_label);
    
    % Save individual FC
    writematrix(FC, 'FC_retest_ind.txt', 'delimiter', 'space');
    
    % Fisher Z-transformed FC for group average
    FC_z = 0.5 * log((1 + FC) ./ (1 - FC));
    FC_retest_group = FC_retest_group + FC_z ./ sub_num;
    
    fprintf('HCP 7T retest subject: %d/%d\n', fileth, sub_num);
end

% Convert group average back to correlation space
FC_retest_group = (exp(2 * FC_retest_group) - 1) ./ (exp(2 * FC_retest_group) + 1);

for i = 1:1:600
    FC_retest_group(i,i) = 0;  % Set diagonal to zero
end

% Save group result
cd('/data/results/FC/');
writematrix(FC_retest_group, 'FC_HCP_7T_retest_group.txt', 'delimiter', 'space');
fprintf('Processing HCP 7T retest group done.\n\n');

end
