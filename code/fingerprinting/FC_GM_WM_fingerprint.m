function FC_GM_WM_fingerprint()

% This function calculates similarity between test and retest FC matrices,performs identity matching, permutation test

%% Calculate Similarity Between Test and Retest FC Matrices %%%%%%%%%%%
Sim_FC = zeros(102, 102);

% Initialize filename storage variable
filename = dir('/data/HCP_7T_test/');

sub_num = length(filename);     % Number of subjects

cd('/data/HCP_7T_test/');

%% Main Loop: Process Each Subject in Test Group %%%%%%%%%%%%%%%%%%%%%
for test_subj = 1:sub_num
    
    % Extract and Process Test FC Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    test_dir = filename(test_subj).name;
    cd(['/data/HCP_7T_test/', test_dir]);
    
    FC = load('FC_test_ind.txt');
    
    % Extract GM-WM matrix
    FC_test = FC(201:600, 1:200);
    FC_test = reshape(FC_test, 80000, 1);    % Vectorize for correlation
    
    % Inner Loop: Compare with Each Retest Subject %%%%%%%%%%%%%%%%%%%%
    for retest_subj = 1:sub_num
        
        retest_dir = filename(retest_subj).name;
        cd(['/data/HCP_7T_retest/', retest_dir]);
        
        FC = load('FC_retest_based_whole_parcellation.txt');
        
        % Extract GM-WM matrix
        FC_retest = FC(201:600, 1:200);
        FC_retest = reshape(FC_retest, 80000, 1);
        
        % Calculate similarity (Pearson correlation)
        Sim_FC(test_subj, retest_subj) = corr(FC_test, FC_retest);
    end
    
    % Display progress
    disp(['FC fingerprint for test subject: ', num2str(test_subj)]);
end

%% Identity Matching %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iden_test = zeros(102, 1);      % Best match for each test subject
iden_retest = zeros(102, 1);    % Best match for each retest subject

for i = 1:102
    % Find highest similarity in each row (test) and column (retest)
    [~, iden_test(i)] = max(Sim_FC(i, :));
    [~, iden_retest(i)] = max(Sim_FC(:, i));
end

% Calculate accuracy
acc_test = 0;
acc_retest = 0;

for i = 1:102
    if i == iden_test(i)
        acc_test = acc_test + 1;
    end
    if i == iden_retest(i)
        acc_retest = acc_retest + 1;
    end
end

% Overall accuracy (average of test and retest directions)
acc = (acc_test + acc_retest) / (2 * 102);

%% Permutation Testing for Significance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
acc_perm = zeros(5000, 1);

for iter = 1:5000
    % Generate random permutation of subject labels
    sub_rand = randperm(102);
    
    acc_t = 0;
    acc_re = 0;
    
    % Calculate accuracy with permuted labels
    for i = 1:102
        if sub_rand(i) == iden_test(i)
            acc_t = acc_t + 1;
        end
        if sub_rand(i) == iden_retest(i)
            acc_re = acc_re + 1;
        end
    end
    
    acc_perm(iter) = (acc_t + acc_re) / (2 * 102);
end

% Save Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('/data/results/fingerprint');

writematrix([acc_test, acc_retest, acc], 'accuracy_raw.txt');
writematrix(acc_perm, 'accuracy_permutation.txt');

end