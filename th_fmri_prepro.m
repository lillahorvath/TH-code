function th_fmri_prepro(fmri_pre_inps)

% This function preprocesses the fMRI data acquired from 19 participants 
% completing the treasure hunt paradigm in a Siemens TrioTim 3T scanner. 
%
% Inputs
%       fmri_pre_inps: structure with fields
%        .data_dir:    group level data directory
%        .func_dir:    functional data subdirectory: participant-level
%        .pre_dir:     fmri preprocessing subdirectory: participant-level
%        .subs:        BIDS style participant IDs
%        .runs:        participant run numbers
%
% Outputs
%       none
%
% Copyright (C) Dirk Ostwald, Lilla Horvath
% -------------------------------------------------------------------------

% initialization
% -------------------------------------------------------------------------

% invoke SPM12
addpath(genpath([pwd '\spm12']))

% unpack input structure
data_dir = fmri_pre_inps.data_dir;
func_dir = fmri_pre_inps.func_dir;
pre_dir  = fmri_pre_inps.pre_dir;
subs     = fmri_pre_inps.subs;
runs     = fmri_pre_inps.runs;

% cycle over participants and preprocess the data
% -------------------------------------------------------------------------
for i = 1:numel(subs)
   
    spm_fmri_th_preprocess(data_dir, subs{i}, func_dir, pre_dir, runs{i})
    
end

end

function spm_fmri_th_preprocess(data_dir, sub, func_dir, pre_dir, run)

% This function implements a series of preprocessing steps for a single
% participant fMRI dataset:
% 
% (0) extract compressed nifti datafiles (nii.gz )
% (1) realignment (estimate and reslice) 
% (2) normalization
% (3) smoothing
% (4) visualization
% 
% Inputs
%       data_dir: group-level data directory
%       sub:      subject ID
%       func_dn:  subject functional data directory
%       pre_dn:   subject fmri preprocessing directory 
%       run_num:  subject specific run numbers
%
% Outputs
%       none, writes proprocessed data to disk
%
% Copyright (C) Dirk Ostwald, Lilla Horvath
% -------------------------------------------------------------------------

% initialization
% -------------------------------------------------------------------------

% set SPM defaults
spm('defaults','fmri')
spm_jobman('initcfg');

% get SPM graphics window 
fg = spm_figure('GetWin','Graphics');

% full path: subject functional data directory
sub_func_dir = [data_dir, filesep, sub, filesep, func_dir, filesep]; 

% full path: preprocessing target directory that will contain the created .mat and image files 
tgt_dir_pre = [data_dir, filesep, sub, filesep, pre_dir, filesep];

% delete target directory and its content if pre-existent and create a new one 
if exist(tgt_dir_pre, 'dir')
    rmdir(tgt_dir_pre,'s')
    mkdir(tgt_dir_pre)
else
    % create target directory if non-existent
    mkdir(tgt_dir_pre)
end

% data preprocessing
% -------------------------------------------------------------------------

% ---------------- (0) data extraction and (1) realignment ----------------

% cycle over runs
for r = 1:length(run)
    
    % if data is in compressed form: unzip
    if ~exist([sub_func_dir, [sub, '_task-th', '_run-0' num2str(r) '_bold.nii']], 'file')
        gunzip([sub_func_dir, [sub, '_task-th', '_run-0' num2str(r) '_bold.nii.gz']]);
    end
     
    % get number of frames 
    N      = [];
    fr_num = [];
    N      = nifti([sub_func_dir, [sub, '_task-th', '_run-0' num2str(r) '_bold.nii']]);
    fr_num = N.dat.dim(4);
     
    % create SPM style filenames 
    f = [];
    f = [sub_func_dir, [sub, '_task-th', '_run-0' num2str(r) '_bold.nii']];
    
    % get list of EPI volumes for realigment 
    fs = cell(fr_num, 1);
    for fr = 1:fr_num
       fs(fr) = cellstr([f ',' num2str(fr)]);
    end
    
    % assign to job structure: realign to the first volume of the first run
    job{1}.spm.spatial.realign.estwrite.data{r} = fs;

end

% job structure specification
job{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;                 % transformation parameter estimation ("est") options - SPM12 default values quality parameter
job{1}.spm.spatial.realign.estwrite.eoptions.sep     = 4;                   % separation 
job{1}.spm.spatial.realign.estwrite.eoptions.fwhm    = 5;                   % smoothing (FWHM)
job{1}.spm.spatial.realign.estwrite.eoptions.rtm     = 0;                   % register to first: images are registered to the first image in the series
job{1}.spm.spatial.realign.estwrite.eoptions.interp  = 2;                   % interpolation polynomial order
job{1}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];             % no warping
job{1}.spm.spatial.realign.estwrite.eoptions.weight  = '';                  % no weighting
job{1}.spm.spatial.realign.estwrite.roptions.which   = [2 1];               % resampling ("write") options - SPM8 default values
job{1}.spm.spatial.realign.estwrite.roptions.interp  = 4;                   % 4th degree B-Spline interpolation
job{1}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];             % no warping
job{1}.spm.spatial.realign.estwrite.roptions.mask    = 1;                   % Masking
job{1}.spm.spatial.realign.estwrite.roptions.prefix  = 'r';                 % create new images with prefix 'r' (realignment)

% save the job variable to disc
save([tgt_dir_pre, 'job_realignment.mat'], 'job');

% inform user
fprintf(['Realignment ', sub, '\n']);

% run job
spm_jobman('run', job);

% save SPM window 
saveas(fg, [tgt_dir_pre, [sub, '_rp.pdf']]);

% delete created .ps file (window-linked)
delete *.ps

% clear job
clear job;

% --------------------------- (2) normalization ---------------------------

% get tissue probability map for inhomegeneity correction
tpm = {[pwd '\spm12\tpm\TPM.nii']};

% initialize filename array
fs_all = [];

% cycle over runs
for r = 1:length(run)

    % get number of frames 
    N      = [];
    fr_num = [];
    N      = nifti([sub_func_dir, [sub, '_task-th', '_run-0' num2str(r) '_bold.nii']]);
    fr_num = N.dat.dim(4);
     
    % create SPM style filenames 
    f = [];
    f = [sub_func_dir, ['r', sub, '_task-th', '_run-0' num2str(r) '_bold.nii']];
    
    % get list of realigned EPI volumes for normalization
    fs = cell(fr_num, 1);
    for fr = 1:fr_num
       fs(fr) = cellstr([f ',' num2str(fr)]);
    end
    
    % concatenate over runs
    fs_all = [fs_all; fs];
    
end

% job structure specification
job{1}.spm.spatial.normalise.estwrite.subj.vol          = {fs_all{1}};      % source image - this image is normalized
job{1}.spm.spatial.normalise.estwrite.subj.resample     = fs_all;           % images to normalize based on parameters estimated for source image
job{1}.spm.spatial.normalise.estwrite.eoptions.biasreg  = 0.0001;           % regulariation
job{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;               % smmohting
job{1}.spm.spatial.normalise.estwrite.eoptions.tpm      = tpm;              % tissue probability map
job{1}.spm.spatial.normalise.estwrite.eoptions.affreg   = 'mni';            % normalize to MNI space
job{1}.spm.spatial.normalise.estwrite.eoptions.reg      = [0 0.001 0.5 0.05 0.2];
job{1}.spm.spatial.normalise.estwrite.eoptions.fwhm     = 0;                % writing smoothing
job{1}.spm.spatial.normalise.estwrite.eoptions.samp     = 3;                % resampling
job{1}.spm.spatial.normalise.estwrite.woptions.bb       = [-78 -112 -70
                                                            78 76 85];      % bounding box
job{1}.spm.spatial.normalise.estwrite.woptions.vox      = [2 2 2];          % voxel size for re-interpoloation
job{1}.spm.spatial.normalise.estwrite.woptions.interp   = 4;                % interpolation constant
job{1}.spm.spatial.normalise.estwrite.woptions.prefix   = 'n';              % file prefix

% save the job variable to disc
save([tgt_dir_pre, 'job_normalization.mat'], 'job');

% inform user
fprintf(['Normalization ', sub, '\n'])

% run job
spm_jobman('run', job)

% clear job
clear job

% ----------------------------- (3) smoothing -----------------------------

% initialize filename array
fs_all = [];

% cycle over runs
for r = 1:length(run)
    
    % get number of frames 
    N      = [];
    fr_num = [];
    N      = nifti([sub_func_dir, [sub, '_task-th', '_run-0' num2str(r) '_bold.nii']]);
    fr_num = N.dat.dim(4);
     
    % create SPM style filenames
    f = [];
    f = [sub_func_dir, ['nr', sub, '_task-th', '_run-0' num2str(r) '_bold.nii']];
    
    % get list of realigned and normalized EPI volumes for smoothing
    fs = cell(fr_num, 1);
    for fr = 1:fr_num
       fs(fr) = cellstr([f ',' num2str(fr)]);
    end
    
    % concatenate over runs
    fs_all = [fs_all; fs];
 
end

% job structure specification
job{1}.spm.spatial.smooth.data   = fs_all;                                  % filenames 
job{1}.spm.spatial.smooth.fwhm   = [8 8 8];                                 % Gaussian smoothing kernel FWHM
job{1}.spm.spatial.smooth.dtype  = 0;                                       % file type
job{1}.spm.spatial.smooth.im     = 0;                                       % no implicit masking
job{1}.spm.spatial.smooth.prefix = 's';                                     % create new images with prefix 's' (smoothing)

% save the job variable to disc
save([tgt_dir_pre, 'job_smoothing.mat'], 'job');

% inform user
fprintf(['Smoothing ', sub, '\n'])

% run job
spm_jobman('run', job)

% clear job
clear job

% --------------------------- (4) visualization ---------------------------

% cycle over runs
for r = 1:length(run)
     
    % get the raw, realigned, normalized and smoothed run specific .nii
    % files
    fs    = cell(4,1);
    fs(1) = cellstr([sub_func_dir, [sub, '_task-th', '_run-0' num2str(r) '_bold.nii'] ',' num2str(1)]);
    fs(2) = cellstr([sub_func_dir, ['r', sub, '_task-th', '_run-0' num2str(r) '_bold.nii'] ',' num2str(1)]);
    fs(3) = cellstr([sub_func_dir, ['nr', sub, '_task-th', '_run-0' num2str(r) '_bold.nii'] ',' num2str(1)]);
    fs(4) = cellstr([sub_func_dir, ['snr', sub, '_task-th', '_run-0' num2str(r) '_bold.nii'] ',' num2str(1)]);
    
    % assign to job structure
    job{1}.spm.util.checkreg.data = fs;
    
    % save the job variable to disc
    save([tgt_dir_pre, 'job_vis_run-0' num2str(r) '.mat'], 'job');

    % run job
    spm_jobman('run', job)
    
    % save SPM window
    saveas(fg, [tgt_dir_pre, ['pre_', sub, '_task-th_run-0' num2str(r) '.pdf']]);
    
    % delete created .ps file (window-linked)
    delete *.ps
    
    % clear job
    clear job
    
end

end