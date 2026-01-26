function Regression_model_neurotrans()

% Performs linear regression to investigate neurotransmitter basis of weighted degree pattern

%%% Load weighted deegree data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('/data/results/FC');

DC = ('DC_HCP_7T_test_group.txt');
DC = zscore(DC);

%%% Load neurotransmitter receptor data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('/data/results/neurotransmitter');

dirpath_tran = '/data/results/neurotransmitter';
filelist_tran = dir(dirpath_tran);
n_tran = length(filelist_tran);  % Total of 8 receptor systems

% Initialize receptor data matrix
neurotrans = zeros(400, n_tran);

% Load each receptor data file
for i = 1:n_tran
    neurotrans(:, i) = zscore(load(filelist_tran(i).name));
end

%%% Fit regression model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform linear regression
mdl_DC = fitglm(neurotrans, DC, 'linear', 'Distribution', 'normal');

% Extract model coefficients
GLM = (mdl_DC.Coefficients.Estimate)';

% Calculate predicted DC values
Pre_DC = sum(GLM .* [ones(400, 1), neurotrans], 2);

% Calculate adjusted R-squared
R_squa = mdl_DC.Rsquared.Adjusted;

% Save observed and predicted DC values
writematrix([DC, Pre_DC], 'Weighted_degree_from_regression_model.txt');

% Save model fit statistics
writematrix(R_squa, 'Adj_R_squa.txt');

end