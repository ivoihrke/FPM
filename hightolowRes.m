%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is to convert the high resolution image
% into low resolution images with the knowledge of their
% corresponding dirac peak positions.
%
% by John Meshreki, john.meshreki@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% helper functions
%Fourier and inverse Transforms
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));

% operator to crop region of O
downsamp = @(x,cen,Np) x(cen(1)-floor(Np(1)/2):cen(1)-floor(Np(1)/2)+Np(1)-1,...
    cen(2)-floor(Np(2)/2):cen(2)-floor(Np(2)/2)+Np(2)-1);

%% Specify dependencies and Date/time tags
% use today's date to create new output directory
todaysdatetime = string(datetime('now','Format','yyyy_MM_dd_HH_mm_ss'));
todaysdate = string(datetime('today','Format','yyyy_MM_dd'));

% add path to dependencies
addpath('../dependencies/natsortfiles');
addpath('../dependencies/export_fig');
addpath('../dependencies/labelpoints');
addpath('../dependencies/min_max_elements_index_n_values');

% add path for functions files
addpath('FP_Func/');

%% specify which sample to use as input
sample_name= 'USAF'; %'stained', 'USAF', 'hela'

%% map container of output directory path
out_dir_name = containers.Map({'stained'; 'USAF'; 'hela'},...
{strcat('../out_dir/Tian14_StainedHistologySlide/',todaysdate,'/',todaysdatetime,'/');...
    strcat('../out_dir/Tian14_ResTarget/',todaysdate,'/',todaysdatetime,'/');...
    strcat('../out_dir/Out_Tian15_inVitroHeLa/',todaysdate,'/',todaysdatetime,'/')});

out_dir = out_dir_name(sample_name);
mkdir(out_dir);

% keep a log
diary(strcat(out_dir,'/','log_',todaysdatetime,'.txt'));

%% read in all images into the memory first
fprintf(['loading the high resolution image...\n']);
tic;

high_res_input_dir_name = containers.Map({'stained'; 'USAF'; 'hela'},...
{'../data/Tian14/1LED/tif/';...
   %'/home/ads/jm095624/microscopy3d/fourier_ptychography/out_dir/Tian14_ResTarget/28_05_2022/28_05_2022_19_13_42/';...
   '/home/ads/jm095624/microscopy3d/fourier_ptychography/out_dir/Tian14_ResTarget/2022_06_29/2022_06_29_17_56_18/';...
   '../data/Tian15_inVitroHeLa/data/'});

high_res_filedir = high_res_input_dir_name(sample_name);

% load the high res file
load([high_res_filedir, 'output_high_res.mat']);


% high res complex field
hig_res_O = O; 

toc;

%% saving the high res image for testing it is loaded correctly
imwrite(uint16(hig_res_O), strcat(out_dir,'high_res_O_image.tif'));
f1 = figure('visible','off');imshow(hig_res_O,[]);
title('(high res O)');
export_fig(f1,strcat(out_dir,'high_res_O_figure.tif'),'-m4');


%% go to Fourier space
hig_res_O_fourier = F(hig_res_O);


%% define ROI, read pupil, and no. of images
Np = [2160, 2560];

P = P;

% No. images
Nimg=293;

% read led indices
idx_led = idx_led;

%read dirac centers
dirac_cen = dirac_cen;

%% Estimate low res from high res

% all intensities for estimated low res
I_all_est_low_res = zeros(Np(1),Np(2),Nimg);

for m = 1:Nimg

    % find Dirac peak position for image m (1 LED per image)
    index_m = find(idx_led==m);
    dirac_cen_pos = dirac_cen(index_m, :);
    

    % crope and apply the circular filter(P) at the dirac peak position 
    O_cropped = downsamp(hig_res_O_fourier,dirac_cen_pos,Np).*P;

    % go to real space
    O_est_low_res = Ft(O_cropped);

    % get estimated low res intensity from Object
    I_est_low_res = abs(O_est_low_res).^2;
    
    % save all estimated low res intensities in one array
    I_all_est_low_res(:,:,m) = I_est_low_res;
    
    %% saving variables to matlab files (save every 10 images together!) 
    if rem(m,10)==0 | m==293
        fprintf('\nfinished m = %2d images\n', m);
        low_res_est_I_matfile = fullfile(strcat(out_dir,'low_res_est_I_till_m_',num2str(m),'.mat'));
        save(low_res_est_I_matfile, 'I_all_est_low_res', '-v7.3');
    end
end    

