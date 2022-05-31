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


%% map container of output directory path to
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


% all image data
hig_res_O = O; 

% remove the white pixels from the unprocessed O
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

% saving the high res image for testing everything is loaded correctly
imwrite(uint16(hig_res_O), strcat(out_dir,'high_res_O_image.tif'));
f1 = figure('visible','off');imshow(hig_res_O,[]);
title('(high res O)');
export_fig(f1,strcat(out_dir,'high_res_O_figure.tif'),'-m4');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply a circle filter on the high res
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fourier and inverse Transforms
F = @(x) fftshift(fft2(x));
Ft = @(x) fftshift(ifft2(ifftshift(x)));

% go to Fourier space
hig_res_O_fourier = Ft(hig_res_O);


% saving the high res image in Fourier space
%imwrite(hig_res_O_fourier, strcat(out_dir,'high_res_O_fourier_image.png'));
%f1 = figure('visible','off');imshow(hig_res_O_fourier);
%title('(high res fourier O)');
%export_fig(f1,strcat(out_dir,'high_res_O_fourier_figure.png'),'-m4');

%% define ROI, dirac positions and read pupil
Np = [2160, 2560];
%W_Na = w_NA;
pupil = P;

%% Get dirac peaks positions
% No. images
Nimg=293;

% all intensities for estimated low res
I_all_est_low_res = zeros(Np(1),Np(2),Nimg);
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
    disp(fn);

    % read low res image
    I_meas_low_res = double(imread(fn));
    
    % save all measured low res intensities in one array
    I_all_meas_low_res(:,:,m) = I_meas_low_res;

    index_m = find(idx_led==m);
    dirac_cen_pos = dirac_cen(index_m, :);

    %% operator to crop region of O
    downsamp = @(x,cen) x(dirac_cen_pos(1)-floor(Np(1)/2):dirac_cen_pos(1)-floor(Np(1)/2)+Np(1)-1,...
        dirac_cen_pos(2)-floor(Np(2)/2):dirac_cen_pos(2)-floor(Np(2)/2)+Np(2)-1);

    % cropy and apply the circular filter(W_Na) at the dirac peak position 
    O_cropped = downsamp(hig_res_O_fourier,dirac_cen_pos).*P;

    % go to real space
    O_est_low_res = F(O_cropped);

    % get estimated low res intensity from Object
    I_est_low_res = abs(O_est_low_res).^2;
    
    % save all estimated low res intensities in one array
    I_all_est_low_res(:,:,m) = I_est_low_res;

    % saving the estimated low res
    imwrite(uint16(I_est_low_res), strcat(out_dir,['I_est_low_res_',int2str(m),'_image.tif']));
    %f_I_est_low_res_147 = figure('visible','off');imshow(I_est_low_res_147, []);
    %title('(estimated low res I 147)');
    %export_fig(f_I_est_low_res_147,strcat(out_dir,'I_est_low_res_147_figure.tif'),'-m4');

    % saving the measured low res
    imwrite(uint16(I_meas_low_res), strcat(out_dir,'I_meas_low_res_',int2str(m),'image.tif'));
    %f_I_meas_low_res = figure('visible','off');imshow(I_meas_low_res, []);
    %title('(measured low res I 147)');
    %export_fig(f_I_meas_low_res,strcat(out_dir,'I_meas_low_res_figure.tif'),'-m4');
end    


% saving variables to matlab files
%low_res_both_meas_est_I_matfile = fullfile(out_dir, ['low_res_both_meas_est_I','.mat']);
%save(low_res_both_meas_est_I_matfile, 'I_all_est_low_res', 'I_all_meas_low_res', '-v7.3');