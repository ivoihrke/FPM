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

%% Reconstruction library locates here
%clear all;
addpath('../dependencies/natsortfiles');
addpath('../dependencies/export_fig');

% add path for functions files
addpath('FP_Func/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO 1: specify the file directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Multiplex image directory
filedir = ['../data/Tian14/1LED/tif/'];
%filedir = ['../data/Tian14_ResTarget/1LED/'];
%filedir = ['../data/Tian15_inVitroHeLa/data/'];


% Generate the image list, in 'tif' image format (depending on your image format)
imglist = dir([filedir,'*.tif']);
nstart = [100, 100];
N = natsortfiles({imglist.name});%sorting the images in 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO 2: specify output folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_dir = ['../out_dir/Tian14_StainedHistologySlide/'];
%out_dir = ['../out_dir/Tian14_ResTarget/'];
%out_dir = ['../out_dir/Out_Tian15_inVitroHeLa/'];
%mkdir(out_dir);

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


f_test = figure('visible','off');imshow(Iall(:,:,128));
title('img of Iall, index==128')
temp=[out_dir,filesep,'Iall_index_128.png'];
saveas(gca,temp);


%% define processing ROI
Np = 344;
%Np = 2048;

%nstart = [1,1];


%% read system parameters
% USAF_Parameter();
%U2OS_Parameter();
test_setup();
%test_setup_tian15_invitrohela();
%% load in data: read in the patch from the memory
Imea = double(Iall(nstart(1):nstart(1)+Np-1,nstart(2):nstart(2)+Np-1,:)); % why arrray 344x344x293??
%Imea = double(Iall(1:n1-1,1:n2-1,:)); % why arrray 344x344x293??

I_uint16 = uint16(Imea(:,:,128))
f_test = figure('visible','off');imshow(I_uint16);
title('img of Imea, index==128; patched into 2048x2048')
temp=[out_dir,filesep,'Imea_index_128.png'];
saveas(gca,temp);


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
[dis_lit2,idx_led] = sort(reshape(illumination_na_used,1,Nled)); %reshape NA of the 293 LEDs to be a vector from 1:293 and sort them while saving their NA-indices to idx_led!
                                                                 %i.e., idx_led is the NA-indicies
                                                                 % while dis_lit2 is
                                                                 % the NA-values
                                                                 
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
    % to the image/object plane represented by idx_u and idx_v
    Nsh_lit(:,m) = idx_u(lit0); % for the LED that was on, get dirac peak position in u-coordinate (of the horizontal) relative to the center of image
    Nsv_lit(:,m) = idx_v(lit0); % for the LED that was on, get dirac peak position in v-coordinate (of the vertical) relative to the center of image
end

% reorder the LED indices and intensity measurements according the previous
% dis_lit
% THIS IS WRONG!!!
% THERE IS NO USE OF dis_lit,
Ns = []; % 3D array where for each LED, the corresponding index of k_u and index of k_v are to be stored
Ns(:,:,1) = Nsv_lit;
Ns(:,:,2) = Nsh_lit;

Imea_reorder = Imea(:,:,idx_led);
Ibk_reorder = Ibk(idx_led);

I_reordered_uint16 = uint16(Imea_reorder(:,:,5))
f_test = figure('visible','off');imshow(I_reordered_uint16);
title('img of Imea_reordered, old_index==128, reordered_index==5; patched into 344x344')
temp=[out_dir,filesep,'Imea_reodered_index_128.png'];
saveas(gca,temp);

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

Ithresh_reorder_uint16 = uint16(Ithresh_reorder(:,:,5))
f_test = figure('visible','off');imshow(Ithresh_reorder_uint16);
title('img of Ithresh_reorder, old_index==128, reordered_index==5; patched into 344x344')
temp=[out_dir,filesep,'Ithresh_reorder_index_128.png'];
saveas(gca,temp);


%load('C:\Users\Muneeb\Desktop\ajmal fpm\laura\multiplexed fpm\Ns_cal289.mat\');
% Ns_reorder = Ns(:,idx_led,:);
% Ithresh_reorder = Ithresh_reorder(:,:,1:89);
Ns_reorder = Ns(:,idx_led,:); %why?
                              % 1st element is empty!
                              % 3rd element has the indicies of the k_u and k_v
                              % now store in 2nd element the NA-indicies
                              % idx_led is the NA-indices, e.g.
                              % idx_led(1) = 147,
                              % illumination_na_used(147)=0.035,
                              % dis_lit2(1)=0.035
                              % Therefore, for each NA-index, store indicies of the corresponding k_u and k_v
clear Imea
%% reconstruction algorithm
% select the index of images that will be used in the processing
Nused = 293;
idx_used = 1:Nused;
I = Ithresh_reorder(:,:,idx_used);
Ns2 = Ns_reorder(:,idx_used,:);


I_input_to_AlterMin_uint16 = uint16(Ithresh_reorder(:,:,5))
f_test = figure('visible','off');imshow(I_input_to_AlterMin_uint16);
title('img of I_input_to_AlterMin, old_index==128, reordered_index==5; patched into 344x344')
temp=[out_dir,filesep,'I_input_to_AlterMin_index_128.png'];
saveas(gca,temp);



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
opts.maxIter = 10; 
opts.minIter = 2;
opts.monotone = 1;
% 'full', display every subroutin,
% 'iter', display only results from outer loop
% 0, no display
opts.display = 'full';%0;%'iter';
opts.out_dir = out_dir

upsamp = @(x) padarray(x,[(N_obj-Np)/2,(N_obj-Np)/2]);
opts.O0 = F(sqrt(I(:,:,1)));
opts.O0 = upsamp(opts.O0);
opts.P0 = w_NA;
opts.Ps = w_NA;
opts.iters = 1;
%opts.mode = 'fourier';
opts.mode = 'real';
opts.scale = ones(Nused,1);
opts.OP_alpha = 1;
opts.OP_beta = 1e3 ;
opts.poscalibrate =0;
opts.calbratetol = 1e-1;
opts.F = F;
opts.Ft = Ft;
opts.StepSize = 0.1;
opts.saveIterResult = 1;
f88 = [];
%% algorithm starts
[O,P,err_pc,c,Ns_cal] = AlterMin(I,[N_obj,N_obj],round(Ns2),opts);

%% save results
fn = ['RandLit-',num2str(numlit),'-',num2str(Nused)];
%save([out_dir,'\',fn],'O','P','err_pc','c','Ns_cal');

%f1 = figure; imagesc(-angle(O),[-.6,1]); axis image; colormap gray; axis off 

I = mat2gray(real(O));
f1 = figure('visible','off');imshow(I);
title('Spatial (O)')
temp=[out_dir,filesep,'spatial_O.png'];
saveas(gca,temp);

f2 = figure('visible','off');imshow(angle(O),[]);
title('Spatial angle(O)')
temp=[out_dir,filesep,'spatial_angle_O.png'];
saveas(gca,temp);

f3 = figure('visible','off');imshow(abs(O),[])
title('Spatial abs(O)')
temp=[out_dir,filesep,'spatial_abs_O.png'];
saveas(gca,temp);

f4 = figure('visible','off'); h=barh(err_pc)
title('err/iter')
xlabel('err')
ylabel('iter')
temp=[out_dir,filesep,'err.png'];
saveas(gca,temp);

% figure(3);imagesc(-angle(O));colormap gray;
fprintf('processing completes\n');
