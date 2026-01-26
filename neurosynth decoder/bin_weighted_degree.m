function bin_weighted_degree()
% Create weighted degree bin mask files

% Load network labels
cd('/data/parcellation');
network_label = load('network label sort.txt');

% Read brain parcellation data
[parcellation, ~, ~, H] = y_ReadAll('brain_parcellation_600.nii');

% Initialize index label array
voxel_count = length(find(parcellation));
index_label = zeros(voxel_count, 1);

% Extract non-zero parcellation values
count = 1;
for i = 1:61
    for j = 1:73
        for k = 1:61
            if parcellation(i, j, k) > 0
                index_label(count) = parcellation(i, j, k);
                count = count + 1;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load weighted degree of each WM network
cd('/data/results/DC_WM_networks');
dirpath = '/data/results/DC_WM_networks';
filelist = dir(dirpath);

n_DC = length(filelist);

for i = 1:n_DC
    
    % Load weighted degree data
    DC = load(filelist(i).name);
    
    % Sort weighted degree values in ascending order
    [~, I] = sort(DC, 'ascend');
    
    % Create 20 masks per file
    for j = 1:20
        DC_mask = zeros(400, 1);
        
        % Select indices for current mask
        mask_indices = I((j-1)*20 + 1 : j*20);
        DC_mask(mask_indices) = 1;
        
        % Initialize 3D mask
        DC_mask_3D = zeros(61, 73, 61);
        
        % Map 1D mask to 3D space using network labels
        for count = 201:600
            label = network_label(count);
            
            for ii = 1:61
                for jj = 1:73
                    for kk = 1:61
                        if parcellation(ii, jj, kk) == label
                            DC_mask_3D(ii, jj, kk) = DC_mask(count - 200);
                        end
                    end
                end
            end
        end
        
        % Save 3D mask as NIfTI file
        cd('/data/results/DC_WM_networks/neurosynth');
        output_name = [filelist(i).name(1:end-4), '_mask', num2str(j), '.nii'];
        y_Write(DC_mask_3D, H, output_name);
    end
end
end