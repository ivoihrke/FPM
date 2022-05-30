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

% specify which sample to use as input
sample_name= 'USAF'; %'stained', 'USAF', 'hela'

% map container of input directory path to
% multiplex reading of images
input_dir_name = containers.Map({'stained'; 'USAF'; 'hela'},...
{'../data/Tian14/1LED/tif/';...
   '/home/ads/jm095624/microscopy3d/fourier_ptychography/out_dir/Tian14_ResTarget/28_05_2022/28_05_2022_19_13_42/';...
   '../data/Tian15_inVitroHeLa/data/'});

filedir = input_dir_name(sample_name);

% Generate the image list, in 'tif' image format (depending on your image format)
imglist = dir([filedir,'O_Np_2160_2560_minIter_1_maxIter_3_image.png']);

% map container of output directory path to
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

% load the big file to get w_NA
load([filedir, 'wna.mat']);


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


W_Na = w_NA;
dirac_cen = [3241,3841];
toc;

% saving the high res image for testing everything is loaded correctly
imwrite(hig_res_O, strcat(out_dir,'high_res_O_image.png'));
f1 = figure('visible','off');imshow(hig_res_O);
title('(high res O)');
export_fig(f1,strcat(out_dir,'high_res_O_figure.png'),'-m4');


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

%% define processing ROI
Np = [2160, 2560];

%% operator to crop region of O
downsamp = @(x,cen) x(dirac_cen(1)-floor(Np(1)/2):dirac_cen(1)-floor(Np(1)/2)+Np(1)-1,...
    dirac_cen(2)-floor(Np(2)/2):dirac_cen(2)-floor(Np(2)/2)+Np(2)-1);

% cropy and apply the circular filter(W_Na) at the dirac peak position 
O_cropped = downsamp(hig_res_O_fourier,dirac_cen).*W_Na;

% go to real space
I_low_res_147 = real(F(O_cropped));


% saving the high res image in Fourier space
imwrite(I_low_res_147, strcat(out_dir,'I_low_res_147_image.png'));
f1 = figure('visible','off');imshow(I_low_res_147);
title('(low res O 147)');
export_fig(f1,strcat(out_dir,'I_low_res_147_figure.png'),'-m4');

