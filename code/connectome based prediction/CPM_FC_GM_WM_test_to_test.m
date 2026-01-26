function CPM_FC_GM_WM_test_to_test()

% This function trains on test session data and validates on retest session data, followed by permutation testing.

%% Load General Intelligence Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('/data/CPM');
Gf = load('Gf_HCP_7T.txt');

%% Load Functional Connectivity Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = dir('/data/HCP_7T_test/');
sub_num = length(filename);

cd('/data/HCP_7T_test/');

% Initialize matrix for all subjects' FC data
FC_all = zeros(102, 80000);

for subj_idx = 1:sub_num
    % Extract and FC data for each subject
    subj_dir = filename(subj_idx).name;
    cd(['/data/HCP_7T_test/', subj_dir]);
    
    FC = load('FC_ind_test.txt');
    
    % Extract specific submatrix and vectorize
    FC_GW = FC(201:600, 1:200);
    FC_GW = reshape(FC_GW, 1, 80000);
    
    FC_all(subj_idx, :) = FC_GW;
end

%% LOOCV Prediction Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gf_pre = zeros(102, 1);     % Store predicted Gf scores

p_thre = 0.01;  % Significance threshold for edge selection

for test_subj = 1:102  % Leave-One-Out Cross-Validation
    % Split data into training and testing
    FC_train = FC_all;
    FC_train(test_subj, :) = [];
    
    Gf_train = Gf;
    Gf_train(test_subj) = [];
    
    FC_test = FC_all(test_subj, :);
    
    % Calculate correlation between Gf and each FC edge
    [r, p] = corr(Gf_train, FC_train);
    
    % Identify significant positive and negative edges
    pos_edges = find(p < p_thre & r > 0);
    neg_edges = find(p < p_thre & r < 0);
    
    % Track edge stability across LOOCV folds
    edge_pos(pos_edges) = edge_pos(pos_edges) + 1;
    edge_neg(neg_edges) = edge_neg(neg_edges) + 1;
    
    % Sum positive and negative edges for each subject
    FC_pos_sum = sum(FC_train(:, pos_edges), 2);
    FC_neg_sum = sum(FC_train(:, neg_edges), 2);
    
    % Fit linear regression model
    b = regress(Gf_train, [FC_pos_sum, FC_neg_sum, ones(101, 1)]);
    
    % Predict held-out subject
    FC_pos_test = sum(FC_test(:, pos_edges), 2);
    FC_neg_test = sum(FC_test(:, neg_edges), 2);
    
    Gf_pre(test_subj) = b(1) * FC_pos_test + b(2) * FC_neg_test + b(3);
end

% Visualize prediction accuracy
figure(1);
plot(Gf, Gf_pre, '*');
xlabel('Observed Gf');
ylabel('Predicted Gf');

[r_raw, p_raw] = corr(Gf, Gf_pre);

%% Permutation Testing for Significance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_iter = 5000;
p_thre = 0.01;

r_perm = zeros(1, n_iter);

for iter = 1:n_iter
    % Permute Gf scores to create null distribution
    Gf_perm = Gf(randperm(102));
    Gf_perm_pre = zeros(102, 1);
    
    for test_subj = 1:102  % LOOCV with permuted labels
        % Split data
        FC_train = FC_all;
        FC_train(test_subj, :) = [];
        
        Gf_train_perm = Gf_perm;
        Gf_train_perm(test_subj) = [];
        
        FC_test = FC_all(test_subj, :);
        
        % Calculate correlations
        [r, p] = corr(Gf_train_perm, FC_train);
        
        % Identify significant edges
        pos_edges = find(p < p_thre & r > 0);
        neg_edges = find(p < p_thre & r < 0);
        
        % Sum edges
        FC_pos_sum = sum(FC_train(:, pos_edges), 2);
        FC_neg_sum = sum(FC_train(:, neg_edges), 2);
        
        % Fit model
        b = regress(Gf_train_perm, [FC_pos_sum, FC_neg_sum, ones(101, 1)]);
        
        % Predict
        FC_pos_test = sum(FC_test(:, pos_edges), 2);
        FC_neg_test = sum(FC_test(:, neg_edges), 2);
        
        Gf_perm_pre(test_subj) = b(1) * FC_pos_test + b(2) * FC_neg_test + b(3);
    end
    
    r_perm(iter) = corr(Gf_perm, Gf_perm_pre);
    
    % Display progress
    disp(['Permutation test iteration: ', num2str(iter)]);
end

%% Calculate P-value from Null Distribution %%%%%%%%%%%%%%%%%%%%%%%%%
p_perm = length(find(r_perm > r_raw)) / n_iter;

%% Save Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

writematrix([Gf, Gf_pre], 'CPM_GM-WM_test_to_test.txt');
writematrix(r_perm, 'r_perm_test_to_test.txt');

end