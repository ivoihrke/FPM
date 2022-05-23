%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main file to implement Fourier Ptychography reconstruction algorithm
% ref
% Lei Tian, et.al, Biomedical Optics Express 5, 2376-2389 (2014).
%
% last modified on 10/07/2015
% by Lei Tian, lei_tian@alum.mit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To do list for the user: (marked by 'TODO#')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) specify where the data located in 'filedir'
% 2) specify where you want to store the results in 'out_dir'
% 3) Find coordinates for estimating background levels and input into 'bk1'
% and 'bk2'.
% 4) specify a threshold value, above which the estimated background value
% is signal rather than noise.
% 5) make sure the LED index (used for taking the images) are properly defined
% in 'lit'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use today's date to create new output directory
todaysdatetime = string(datetime('now','Format','dd_MM_yyyy_HH_mm_ss'));
todaysdate = string(datetime('today','Format','dd_MM_yyyy'));



%% Reconstruction library locates here
%clear all;
addpath('../dependencies/natsortfiles');
addpath('../dependencies/export_fig');
addpath('../dependencies/labelpoints');


% add path for functions files
addpath('FP_Func/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO 1: specify the file directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify which sample to use as input
sample_name= 'hela'; %'stained', 'USAF', 'hela'

% map container of input directory path to
% multiplex reading of images
input_dir_name = containers.Map({'stained'; 'USAF'; 'hela'},...
{'../data/Tian14/1LED/tif/';...
   '../data/Tian14_ResTarget/1LED/';...
   '../data/Tian15_inVitroHeLa/data/'});

filedir = input_dir_name(sample_name);

% Generate the image list, in 'tif' image format (depending on your image format)
imglist = dir([filedir,'*.tif']);
nstart = [100, 100];
N = natsortfiles({imglist.name});%sorting the images in 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO 2: specify output folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% map container of output directory path to
out_dir_name = containers.Map({'stained'; 'USAF'; 'hela'},...
{strcat('../out_dir/Tian14_StainedHistologySlide/',todaysdate,'/',todaysdatetime,'/');...
    strcat('../out_dir/Tian14_ResTarget/',todaysdate,'/',todaysdatetime,'/');...
    strcat('../out_dir/Out_Tian15_inVitroHeLa/',todaysdate,'/',todaysdatetime,'/')});

out_dir = out_dir_name(sample_name);
mkdir(out_dir);

% keep a log
diary(strcat(out_dir,'/','log_',todaysdatetime,'.txt'));


%% define # of LEDs used to capture each image
numlit = 1;
% raw image size
n1 = 2160; n2 = 2560;


%% read in all images into the memory first
fprintf(['loading the images...\n']);
tic;
Nimg = length(imglist);
Iall = zeros(n1,n2,Nimg,'uint16'); %3d array; each pixel in the 2d image for all images
Ibk = zeros(Nimg,1);
for m = 1:Nimg

    %     fn = [filedir,imglist(m).name]; %removed for sorting purpose
    fn = [filedir,N{m}];  %added for sorting purpose
    disp(fn);
    % all image data
    I = double(imread(fn));
    Iall(:,:,m) = I  ; %filling 3d array by replacing zeros in 1st&2nd elements by the image values and 3rd element is index of image
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TODO 3: specify background region coordinates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % backgr4ound noise esimtation 
    bk1 = mean2(double(Iall(1:100,1:100,m)));
    bk2 = mean2(double(Iall(492:600,2380:2520,m)));

    Ibk(m) = mean([bk1,bk2]);
    % brightfield background processing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO 4: if Ibk is larger than some threshold, it is not noise! (chnage 
    % this value correpondingly)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if Ibk(m)>300*(2^8-1)/(2^16-1)&& m>1
%         Ibk(m) = Ibk(m-1);
%     end
    if Ibk(m)>300 
        Ibk(m) = Ibk(m-1);
    end
end
fprintf(['\nfinish loading images\n']);
toc;


%f_test = figure('visible','off');imshow(Iall(:,:,128));
%title('img of Iall, index==128');
%export_fig(f_test,strcat(out_dir,'Iall_index_128.png'),'-m4');

%% define processing ROI
%Np = 344;
%Np = 1000;
%Np = 2048;
%Np = 2160;
%Np = [344,344];
Np = [2160, 2560]


%nstart = [1,1];

%% read system parameters
USAF_Parameter();
%if sample_name == 'USAF'
%    USAF_Parameter();
%else
%    test_setup();
%end

% overwrite nstart in USAF_Parameter.m
nstart = [1,1];

%% load in data: read in the patch from the memory
Imea = double(Iall(nstart(1):nstart(1)+Np(1)-1,nstart(2):nstart(2)+Np(2)-1,:)); % why arrray 344x344x293??
%Imea = double(Iall(1:n1-1,1:n2-1,:)); % why arrray 344x344x293??

%f_test = figure('visible','off');imshow(uint16(Imea(:,:,128)));
%title('img of Imea, index==128; patched into 2048x2048');
%export_fig(f_test,strcat(out_dir,'Imea_index_128.png'),'-m4');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO 5: Define the LED index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in this example, the LEDs are scanned sequentially within illumination NA
% defined in the systemsetup();
% In other cases, make the correponding changes. In the end, all we need is
% the matrix 'lit' that specifies the LED index corresponding image.

ledidx = 1:Nled;
ledidx = reshape(ledidx,numlit,Nimg);
lit = Litidx(ledidx); %index of the sequence of LEDs that were used in experiment corresponding to the 293 LEDs in the 32x32 grid!
lit = reshape(lit,numlit,Nimg);
% reorder LED indices based on illumination NA
[dis_lit2,idx_led] = sort(reshape(illumination_na_used,1,Nled)); %reshape NA of the 293 LEDs to be a vector from 1:293 and sort them so that
                                                                 %they store SEPARATELY the index and values of the NA after sorting them!
                                                                 %i.e., idx_led is the SORTED NA-indicies
                                                                 % while dis_lit2 is
                                                                 % the SORTED NA-values
                                                                 
Nsh_lit = zeros(numlit,Nimg); % a vector (in case of numlit==1) of size Nimg==293
Nsv_lit = zeros(numlit,Nimg);

for m = 1:Nimg
    % corresponding index of spatial freq for the LEDs are lit
    lit0 = lit(:,m); %read index of sequence of LEDs that were used in experiment for each image
                     %in case of 1 LED sequence; lit0 just stores the
                     %index of that LED that was on for the corresponding
                     %image!
                     %Therefore, for that image, which LED was on!!
                     
    % this is a translation from the LED plane represented by lit0 or lit
    % to the pupil plane represented by idx_u and idx_v
    Nsh_lit(:,m) = idx_u(lit0); % for the LED that was on, get dirac peak position in u-coordinate (of the horizontal) relative to the center of image
    Nsv_lit(:,m) = idx_v(lit0); % for the LED that was on, get dirac peak position in v-coordinate (of the vertical) relative to the center of image
end

% reorder the LED indices and intensity measurements according the previous
% dis_lit or rather their indices: idx_led
Ns = []; % 3D array where for each LED, the corresponding index of k_u and index of k_v are to be stored
Ns(:,:,1) = Nsv_lit;
Ns(:,:,2) = Nsh_lit;

% Images indexed by LEDs indicies at illumination plane
%f_image_illum_plane = figure('visible','off'); 
%pcolor(LitCoord);
%labelpoints(LitCoord, string(illumination_na_used), 'FontSize', 5);
%title(['Images pattern indexed by LEDs indicies at illumination plane for all 32x32 LEDs/n with their NAs']);
%export_fig(f_image_illum_plane,strcat(out_dir,'image_indexed_pattern_illum_plane.png'),'-m4');


%Dirac centers of the 32x32 LEDs
f_dirac_32 = figure('visible','off'); scatter(idx_u, idx_v);
title(['Dirac peaks centers at O plane for all 32x32 LEDs']);
export_fig(f_dirac_32,strcat(out_dir,'dirac_centers_32_32.png'),'-m4');

%Dirac centers of 293 LEDs sorted according to their NA
f_dirac_293= figure('visible','off'); 
scatter(Nsh_lit, Nsv_lit);
labelpoints(Nsh_lit, Nsv_lit, string(idx_led), 'FontSize', 5);
title(['Dirac peaks centers at O plane for all 293 LEDs after NA sorting']);
export_fig(f_dirac_293,strcat(out_dir,'dirac_centers_293_NA_sorted.png'),'-m4');

Imea_reorder = Imea(:,:,idx_led);
Ibk_reorder = Ibk(idx_led);

%f_test = figure('visible','off');imshow(uint16(Imea_reorder(:,:,5)));
%title('img of Imea_reordered, old_index==128, reordered_index==5; patched into 344x344');
%export_fig(f_test,strcat(out_dir,'Imea_reodered_index_128.png'),'-m4');

% gain = a(idx_led);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processing the data to DENOISING is IMPORTANT
% background subtraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ithresh_reorder = Imea_reorder;
for m = 1:Nimg
    Itmp = Ithresh_reorder(:,:,m);
    Itmp = Itmp-Ibk_reorder(m);
%     Itmp = awgn(Itmp,0,'measured');
    Itmp(Itmp<0) = 0;
    Ithresh_reorder(:,:,m) = Itmp;
    
end

%f_test = figure('visible','off');imshow(uint16(Ithresh_reorder(:,:,5)));
%title('img of Ithresh_reorder, old_index==128, reordered_index==5; patched into 344x344');
%export_fig(f_test,strcat(out_dir,'Ithresh_reorder_index_128.png'),'-m4');

%load('C:\Users\Muneeb\Desktop\ajmal fpm\laura\multiplexed fpm\Ns_cal289.mat\');
% Ns_reorder = Ns(:,idx_led,:);
% Ithresh_reorder = Ithresh_reorder(:,:,1:89);
Ns_reorder = Ns(:,idx_led,:); %why?
                              % 1st element is empty!
                              % 3rd element has the indicies of the k_u and k_v
                              % now store in 2nd element the NA-indicies
                              % idx_led is the NA-indices, e.g.
                              % for USAF samples, the following:
                              % idx_led(1) = 147,
                              % illumination_na_used(147)=0.0012,
                              % dis_lit2(1)=0.0012
                              % Therefore, for each NA-index, store indicies of the corresponding k_u and k_v
clear Imea;
%% reconstruction algorithm
% select the index of images that will be used in the processing
Nused = 293;
idx_used = 1:Nused;
I = Ithresh_reorder(:,:,idx_used);
Ns2 = Ns_reorder(:,idx_used,:);


f_test = figure('visible','off');imshow(uint16(Ithresh_reorder(:,:,5)));
title('img of I_input_to_AlterMin, old_index==128, reordered_index==5; patched into 344x344');
export_fig(f_test,strcat(out_dir,'I_input_to_AlterMin_index_128.png'),'-m4');

%% reconstruction algorithm options: opts
%   tol: maximum change of error allowed in two consecutive iterations
    %   maxIter: maximum iterations 
    %   minIter: minimum iterations
    %   monotone (1, default): if monotone, error has to monotonically dropping
    %   when iters>minIter
%   display: display results (0: no (default) 1: yes)
    %   saveIterResult: save results at each step as images (0: no (default) 1: yes)
    %   mode: display in 'real' space or 'fourier' space.
    %   out_dir: saving directory
%   O0, P0: initial guesses for O and P
    %   OP_alpha: regularization parameter for O
    %   OP_beta: regularization parameter for P
%   scale: LED brightness map
%   H0: known portion of the aberration function, 
        % e.g. sample with a known defocus induce a quardratic aberration
        % function can be defined here
%   poscalibrate: flag for LED position correction using
    % '0': no correction
    % 'sa': simulated annealing method
        % calbratetol: parameter in controlling error tolence in sa
    % 'ga': genetic algorithm
    % caution: takes consierably much longer time to compute a single iteration
%   F, Ft: operators of Fourier transform and inverse
opts.tol = 1;
%opts.maxIter = 1;
opts.maxIter = 3; 
%opts.maxIter = 10; 
%opts.maxIter = 20;
opts.minIter = 1;
%opts.minIter = 2;
%opts.minIter = 20;
%opts.minIter = 10;
opts.monotone = 1;
% 'full', display every subroutin,
% 'iter', display only results from outer loop
% 0, no display
opts.display = 'full';%0;%'iter';
opts.out_dir = out_dir;

upsamp = @(x) padarray(x,[(N_obj(1)-Np(1))/2,(N_obj(2)-Np(2))/2]);
opts.O0 = F(sqrt(I(:,:,1)));
opts.O0 = upsamp(opts.O0);
opts.P0 = w_NA;
opts.Ps = w_NA;
opts.iters = 1;
opts.mode = 'fourier';
%opts.mode = 'real';
opts.scale = ones(Nused,1);
opts.OP_alpha = 1;
opts.OP_beta = 1e3 ;
opts.poscalibrate =0; % using correction or not
opts.calbratetol = 1e-1;
opts.F = F;
opts.Ft = Ft;
opts.StepSize = 0.1;
opts.saveIterResult = 1;
disp(opts);
disp(['nstart: ', num2str(nstart)]);
disp(['Np: ', num2str(Np)]);
f88 = [];
%% algorithm starts
[O,P,err_pc,c,Ns_cal] = AlterMin(I,[N_obj(1),N_obj(2)],round(Ns2),opts);

%% save results
%fn = ['RandLit-',num2str(numlit),'-',num2str(Nused)];
%save([out_dir,'\',fn],'O','P','err_pc','c','Ns_cal');

%f1 = figure; imagesc(-angle(O),[-.6,1]); axis image; colormap gray; axis off 

% add tag for Np + maxIter
Np_Iter = ['Np_',num2str(Np(1)),'_',num2str(Np(2)),'_minIter_',num2str(opts.minIter),'_maxIter_',num2str(opts.maxIter)];

%I = mat2gray(real(O));
real_O = mat2gray(real(O));
f1 = figure('visible','off');imshow(real_O);
title(['(', opts.mode, ')',' (O)']);
export_fig(f1,strcat(out_dir,'O_',Np_Iter,'_',opts.mode,'.png'),'-m4');

angle_O = angle(O);
f2 = figure('visible','off');imshow(angle_O,[]);
title(['(', opts.mode, ')',' angle (O)']);
export_fig(f2,strcat(out_dir,'angle_O_',Np_Iter,'_',opts.mode,'.png'),'-m4');

abs_O = abs(O);
f3 = figure('visible','off');imshow(abs_O,[]);
title(['(', opts.mode, ')',' abs (O)']);
export_fig(f3,strcat(out_dir,'abs_O_',Np_Iter,'_',opts.mode,'.png'),'-m4');

f4 = figure('visible','off'); 
plot(1:length(err_pc), (err_pc))
title(['(', opts.mode, ')', ' err/iter']);
xlabel('err');
ylabel('iter');
export_fig(f4,strcat(out_dir,'err_',Np_Iter,'_',opts.mode,'.png'),'-m4');

f_err_log = figure('visible','off');
plot(1:length(err_pc), (log(err_pc)))
title(['(', opts.mode, ')', ' err/iter']);
xlabel('log err');
ylabel('iter');
export_fig(f_err_log,strcat(out_dir,'err_log_',Np_Iter,'_',opts.mode,'.png'),'-m4');

% figure(3);imagesc(-angle(O));colormap gray;

% saving variables to matlab files
hig_res_O_matfile = fullfile(out_dir, ['hig_res_O_',Np_Iter,'_',opts.mode,'.mat']);
save(hig_res_O_matfile, 'O', 'real_O');
all_vars_matfile = fullfile(out_dir, ['all_vars_',Np_Iter,'_',opts.mode,'.mat']);
save(all_vars_matfile);
% save all but these large variables
%save(all_vars_matfile, '-regexp','^(?!(I|Iall|Imea_reorder|Ithresh_reorder)$).');



%post-processing of the white pixels w/ high intensity
proc_abs_O = abs(O);
%proc_abs_O(abs(O)>25) = 25;
proc_abs_O(abs(O)>18) = 18;
f5 = figure('visible','off');imshow(proc_abs_O,[]);
title(['(', opts.mode, ')',' processed abs (O)']);
export_fig(f5,strcat(out_dir,'proc_abs_O_',Np_Iter,'_',opts.mode,'.png'),'-m4');


fprintf('processing completes\n');
