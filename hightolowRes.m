%% helper functions
%Fourier and inverse Transforms
F = @(x) fftshift(fft2(x));
Ft = @(x) fftshift(ifft2(ifftshift(x)));

% operator to crop region of O
downsamp = @(x,cen,Np) x(cen(1)-floor(Np(1)/2):cen(1)-floor(Np(1)/2)+Np(1)-1,...
    cen(2)-floor(Np(2)/2):cen(2)-floor(Np(2)/2)+Np(2)-1);

%% Specify dependencies and Date/time tags
% use today's date to create new output directory
todaysdatetime = string(datetime('now','Format','dd_MM_yyyy_HH_mm_ss'));
todaysdate = string(datetime('today','Format','dd_MM_yyyy'));

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
   '/home/ads/jm095624/microscopy3d/fourier_ptychography/out_dir/Tian14_ResTarget/28_05_2022/28_05_2022_19_13_42/';...
   '../data/Tian15_inVitroHeLa/data/'});

high_res_filedir = high_res_input_dir_name(sample_name);

% load the big file to get w_NA
load([high_res_filedir, 'output_high_res.mat']);


% high res complex field
hig_res_O = O; 

%% remove the white pixels from the unprocessed O
nrows = size(O, 1);
ncols = size(O, 2);

for i = 1:nrows
    for j = 1:ncols
      if abs(hig_res_O(i,j))>25
        hig_res_O(i,j) = complex(25, imag(hig_res_O(i,j)));
      end
    end
end

toc;

%% saving the high res image for testing everything is loaded correctly
imwrite(uint16(hig_res_O), strcat(out_dir,'high_res_O_image.tif'));
f1 = figure('visible','off');imshow(hig_res_O,[]);
title('(high res O)');
export_fig(f1,strcat(out_dir,'high_res_O_figure.tif'),'-m4');


%% go to Fourier space
hig_res_O_fourier = Ft(hig_res_O);


%% define ROI, read pupil, and no. of images
Np = [2160, 2560];
%W_Na = w_NA;

pupil = P;

% No. images
Nimg=293;

%% Estimate low res from high res

% all intensities for estimated low res
I_all_est_low_res = zeros(Np(1),Np(2),Nimg);

% all Object complex fields for estimated low res
O_all_est_low_res = zeros(Np(1),Np(2),Nimg);

% all intensities for measured low res
I_all_meas_low_res = zeros(Np(1),Np(2),Nimg);

for m = 1:Nimg
    % map container of input directory path to
    % multiplex reading of low res image
    low_res_input_dir_name = containers.Map({'stained'; 'USAF'; 'hela'},...
    {'../data/Tian14/1LED/tif/';...
       '../data/Tian14_ResTarget/1LED/';...
       '../data/Tian15_inVitroHeLa/data/'});

    low_res_filedir = low_res_input_dir_name(sample_name);

    % Generate the image list, in 'tif' image format (depending on your image format)
    low_res_imglist = dir([low_res_filedir,'Iled_',int2str(m),'.tif']);

    fn = [low_res_filedir,low_res_imglist.name];

    % read low res image
    I_meas_low_res = double(imread(fn));
    
    % save all measured low res intensities in one array
    I_all_meas_low_res(:,:,m) = I_meas_low_res;

    % find Dirac peak position for image m (1 LED per image)
    index_m = find(idx_led==m);
    dirac_cen_pos = dirac_cen(index_m, :);
    

    % crope and apply the circular filter(P) at the dirac peak position 
    O_cropped = downsamp(hig_res_O_fourier,dirac_cen_pos,Np).*P;

    % go to real space
    O_est_low_res = F(O_cropped);
    
    % save all estimated low res Object complex fields in one array
    O_all_est_low_res(:,:,m) = O_est_low_res;

    % get estimated low res intensity from Object
    I_est_low_res = abs(O_est_low_res).^2;
    
    % save all estimated low res intensities in one array
    I_all_est_low_res(:,:,m) = I_est_low_res;

    % saving the estimated low res intensities
    imwrite(uint16(I_est_low_res), strcat(out_dir,['I_est_low_res_',int2str(m),'_image.tif']));

    % saving the measured low res intensities
    imwrite(uint16(I_meas_low_res), strcat(out_dir,'I_meas_low_res_',int2str(m),'_image.tif'));
end    

%% saving variables to matlab files (careful: huge output file, 
% best to put an if cond when adding to the arrays to select a few images!)
%low_res_both_meas_est_I_matfile = fullfile(out_dir, ['low_res_both_meas_est_I','.mat']);
%save(low_res_both_meas_est_I_matfile, 'I_all_est_low_res', 'I_all_meas_low_res', '-v7.3');