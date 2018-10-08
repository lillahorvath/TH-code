function lumin = th_luminance(lumin)

% This function computes the luminance values of the visual stimuli
% presented to the subjects on the treasure hunt task using the formula 
%           luminance(image) := 0.299*R + 0.587*G + 0.114*B
%
% Inputs
%       lumin:       structure with fields
%        .stim_dir:  directory containing the presented visual stimuli 
%        .grid_ims:  image names: grid cell position
%        .treas_im:  image name: treasure
%        .obs_im:    image names: observation bars
%        .arrow_ims: image names: decision prompt arrows
%
% Outputs
%       lumin:       input structure with additional fields
%        .grid_lum:  luminance values of the grid cell position images
%        .treas_lum: luminance value of the treasure image
%        .obs_lum:   luminance values of the observation bar images
%        .arrow_lum: luminance values of the decision prompt arrow images
%
% Copyright (C) Lilla Horvath, Dirk Ostwald
% -------------------------------------------------------------------------

% initialize
% ------------------------------------------------------------------------

% unpack input structure
stim_dir    = lumin.stim_dir;
grid_fname  = lumin.grid_ims;
treas_fname = lumin.treas_im;
obs_fname   = lumin.obs_ims;
arrow_fname = lumin.arrow_ims;

% compute luminance
% -------------------------------------------------------------------------

% ---------------------- load stimuli into workspace ----------------------

% grid position 
grid_pict = cell(1,numel(grid_fname));
for i = 1:numel(grid_fname)
    grid_pict{i} = double(imread(fullfile(stim_dir,grid_fname{i})));
end

% treasure
treas_pict = double(imread(fullfile(stim_dir,treas_fname)));
                                               
% observation bars
obs_pict = cell(1,numel(obs_fname));
for i = 1:numel(obs_fname)
    obs_pict{i} = double(imread(fullfile(stim_dir,obs_fname{i})));
end

% decision prompt arrows
arrow_pict = cell(1,numel(obs_fname));
for i = 1:numel(obs_fname)
    arrow_pict{i} = double(imread(fullfile(stim_dir,arrow_fname{i})));
end

% ---------------------- compute relative luminance -----------------------

% grid position
grid_lum = NaN(numel(grid_fname),1);
for i = 1:numel(grid_fname)
    grid_lum(i) = mean(mean(grid_pict{i}(:,:,1)*0.299+grid_pict{i}(:,:,2)*0.587+grid_pict{i}(:,:,3)*0.114));
end

% treasure
treas_lum = mean(mean(treas_pict(:,:,1)*0.299+treas_pict(:,:,2)*0.587+treas_pict(:,:,3)*0.114));

% observation bars
obs_lum = NaN(numel(obs_fname),1);
for i = 1:numel(obs_fname)
    obs_lum(i) = mean(mean(obs_pict{i}(:,:,1)*0.299+obs_pict{i}(:,:,2)*0.587+obs_pict{i}(:,:,3)*0.114));
end

% decision prompt arrows
arrow_lum = NaN(numel(arrow_fname),1);
for i = 1:numel(arrow_fname)
    arrow_lum(i) = mean(mean(arrow_pict{i}(:,:,1)*0.299+arrow_pict{i}(:,:,2)*0.587+arrow_pict{i}(:,:,3)*0.114));
end

% -------------------------- set output structure -------------------------

lumin.grid_lum  = grid_lum;
lumin.treas_lum = treas_lum;
lumin.obs_lum   = obs_lum;
lumin.arrow_lum = arrow_lum;

end

