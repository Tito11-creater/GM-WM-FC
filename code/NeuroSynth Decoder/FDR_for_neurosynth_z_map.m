function FDR_for_neurosynth_z_map()

% Loads z-value NeuroSynth data, converts to p-values and applies FDR correction

% Get list of z-value files
cd('/data/results/DC_WM_networks/neurosynth');
filename = dir('/data/results/DC_WM_networks/neurosynth/*z.txt');

% Initialize matrix for z-values
z = zeros(24, 7);

% Load z-values from all files
for i = 1:7
    z(:, i) = load(filename(i).name);
end

% Convert z-values to p-values
p = 2 * (1 - normcdf(abs(z)));

% Reshape p-values for FDR correction
p_vec = reshape(p, [24 * 7, 1]);

% Apply FDR correction
p_fdr = mafdr(p_vec, 'BHFDR', true);

% Reshape back to original dimensions
p_fdr = reshape(p_fdr, [24, 7]);

% Save p-values and FDR-corrected p-values
for i = 1:7
    % Save raw p-values
    writematrix(p(:, i), [filename(i).name(1:end-5), 'p.txt']);
    
    % Save FDR-corrected p-values
    writematrix(p_fdr(:, i), [filename(i).name(1:end-5), 'p_fdr.txt']);
end

end