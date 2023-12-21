function [alldoflag]= ReconHiLo(waveL,raw_file, wf_file,hilo_file,stacknum_file,isstackwrite)
alldoflag=0;disp('1!');
N = 2;
input_parameters=[waveL(2),waveL(3)];
filename = raw_file;
%raw_data = bfopen(raw_file);
if isequal(filename,0)
    disp('Cancel data reading！');
    return;
end
raw_data = readMTiff(raw_file);
nz = size(raw_data,3)/N;
%nz = size(raw_data{1},1)/6;

[Np1,Np2] = size(raw_data(:,:,1)); % 只备注1024
NPixel = max(Np1,Np2);
calibrationFlag = false;
if calibrationFlag == 1
    % calibration map is saved as ratio*10000 in pixel
    calib10 = double(imread(calib1path));
    calib20 = double(imread(calib2path));
    calib1 = zeros(NPixel,NPixel);
    calib2 = zeros(NPixel,NPixel);
    if Np1 == Np2 || Np2 > Np1
       calib1(1:Np1,:)=calib10(1:Np1,:);
    else
       calib2(:,1:Np2)=calib20(:,1:Np2);
    end 
else
    calib1 = ones(NPixel,NPixel);
    calib2 = ones(NPixel,NPixel);
end
calib(:,:,1) = calib1; calib(:,:,2) = calib2;

%[M,N]=size(raw_data[1][1,1]);
%calib1 = ones(M,N); calib2 = ones(M,N);
%fov = size(raw_data{1}{1,1},1);
%calib{1} = imresize(calib1,[fov,fov]);
%calib{2} = imresize(calib2,[fov,fov]);

% input parameters
% NA = 1.49;                      % objective NA
                   % wavelength(nm) of the emission light
% pixel_size = 65;                % pixel size(nm) corresponding to the sample
M_pattern = 0.6;                % sin pattern contrast, for calculating scaling factor
Lo_scaling = 1;                 
% Lateral_resolution = 15.385;     % Lateral resolution in px/um
If_scaling_own = 0;             % whether to select the scaling factor manually
Lo_scaling_own = 0.4;           % manually selected scaling factor


%% iter
for kk =  1 : nz
    raw_sim = zeros(NPixel,NPixel,2);

    for ll = 1 : 2
        raw_sim(1:Np1,1:Np2,ll) = raw_data(:,:,(kk-1)*2+ll);
%         sum = sum +  raw_sim(:,:,ll);
    end
    [M,N]=size(raw_sim);
    uniform_file = zeros(M,N/2,3);
    for ii = 1
        uniform_file(:,:,ii) = (raw_sim(:,:,2*ii-1) + raw_sim(:,:,2*ii))./2;
        img = uniform_file(:,:,ii);
       %% Define parameters 
       % res = Lateral_resolution*1e6*0.0254;    % Lateral resolution
       % uniform_file = (raw_data{1}{1,1} + raw_data{1}{2,1} + raw_data{1}{3,1} + raw_data{1}{4,1}+ raw_data{1}{5,1} + raw_data{1}{6,1})/6;
       [h, w] = size(img);
%        res = 0.61 * waveL / (NA);                   % resolution
%        k_m = w / (res / pixel_size);           % objective cut-off frequency ???
       k_m=w/(10*input_parameters(1));
       kc = nearest(k_m * 0.18);                % cut-off frequency between hp and lp filter
       k_side_HP = kc/2;                     % cutoff frequency of high-pass filter for differential image
       sigmaLP = kc*2/2.355;                   % Finding sigma value for low pass
       lambda = nearest(w/(2*kc));             % Side length of contrast
       % evalutation mask by rounding to
       % nearest odd integer
       if mod(lambda,2) == 0                   % lambda must be odd
           lambda = lambda+1;
       else
       end
       If_correction=0;
       if If_correction == 0                   % Check if correction file is selected
           gauss_corr = 1;
       else
           gauss_corr = single(imread(dir_corr));	% Create gaussian field illumination correction matrix if file was selected
           gauss_corr = gauss_corr+0.2*max(gauss_corr(:));
           gauss_corr = gauss_corr/max(gauss_corr(:));
       end
       nh = ones(lambda);                      % create kernel for local contrast
       lp_conv2 = ones(lambda*2);
       h = h+2*lambda;                         % increase image size by lambda for
       w = w+2*lambda;                         % padding
       Nk = sum(nh(:));                        % evalutation
       %% Creating filters, pre-allocating and reading image stacks
       
       % Create band pass, high and low pass filters
       lp = lpgauss(h,w,sigmaLP);
       hp = hpgauss(h,w,sigmaLP);
       hp_1 = hpgauss(h,w,sigmaLP * 2);
       hp_side = hp_OneSide(h,w,k_side_HP);
       u = single(uniform_file(:,:,ii));
       uni = padarray(u./gauss_corr,[lambda lambda],'symmetric');
       s = single(raw_sim(:,:,2*ii-1));
       rat = padarray(s./gauss_corr,[lambda lambda],'symmetric');
%        rat = rat./uni;
       %% One-side highpass filtering ratio image stack and Contrast evaluation
       rat_hp_s_f = fft2(rat).*hp_side;
       rat_hp_s_f = rat_hp_s_f .* hp_1;
       rat_hp_s = ifft2(rat_hp_s_f);
       weight = (sqrt(rat_hp_s .* conj(rat_hp_s)));
       intermediate_weight = weight(lambda+1:end-lambda,lambda+1:end-lambda,:);
       weight = weight/max(intermediate_weight(:));
       Lo = weight .* uni;
       %% Filtering
       Hi = real(ifft2(fft2(uni).*hp));
       Hi = Hi(lambda+1:end-lambda,lambda+1:end-lambda,:);
       Lo = real(ifft2(fft2(Lo).*lp));
       Lo = Lo(lambda+1:end-lambda,lambda+1:end-lambda,:);
       Lo_raw = real(ifft2(fft2(uni).*lp));          % for debugging
       Lo_raw = Lo_raw(lambda+1:end-lambda,lambda+1:end-lambda,:);
       %% Scaling
       PSFe = fspecial('gaussian',30,3.0);
       Hi=edgetaper(Hi,PSFe);
        Hif=fft2(Hi);Lof=fft2(Lo);
        nabla =abs((sum(sum(Hif)))./(sum(sum(Lof))));
       LoWeight = nabla*input_parameters(2);
       % % method two
       % nabla = abs(fft2(Hi(1,kc+1,:)))./(abs(fft2(Lo(1,kc+1,:))));
       % nabla_mu = mean(nabla);
       % nabla_med = median(nabla);
       % nabla_min = min(nabla);
       % LoWeight = Lo_scaling * nabla_med;
       
       % set manually
       if If_scaling_own
           LoWeight = Lo_scaling_own;
       end
       Lo = Lo*LoWeight;
       
       %% Reconstruction and 16bit conversion
       HiLo = Hi+Lo;
       HiLo(HiLo<0)=0;
       HiLo = HiLo/max(HiLo(:));
       HiLo_dir(:,:,ii) = HiLo;
%        HiLo_dir(:,:,ii) = im2uint16(HiLo);
    end
%     HiLo01 = HiLo_dir(:,:,1);
%     imwrite(HiLo01, [dataname '6张图重建HiLo\HiLo01_', num2str(kk, '%.3d'), '.tif'])
%     HiLo02 = HiLo_dir(:,:,2);
%     imwrite(HiLo02, [dataname '6张图重建HiLo\HiLo02_', num2str(kk, '%.3d'), '.tif'])
%     HiLo03 = HiLo_dir(:,:,3);
%     imwrite(HiLo03, [dataname '6张图重建HiLo\HiLo03_', num2str(kk, '%.3d'), '.tif'])
%     HiLo_dir(:,:,1) =  HiLo_dir(:,:,1)./max(max(HiLo_dir(:,:,1)));
%     HiLo_dir(:,:,2) =  HiLo_dir(:,:,2)./max(max(HiLo_dir(:,:,2)));
%     HiLo_dir(:,:,3) =  HiLo_dir(:,:,3)./max(max(HiLo_dir(:,:,3)));
    HiLo_sum = HiLo_dir(:,:,1);
    HiLo_sum = HiLo_sum/ max(HiLo_sum(:));
    HiLo = im2uint16(HiLo_sum);
       
   
    HiLo_dir_uint = im2uint16(HiLo_dir);
    
    WF = sum(raw_sim,3)/2;
   %% saving
  WF1=WF;
  HiLo1=HiLo;
  clear WF HiLo
  Hilo=zeros(Np1,Np2); WF=zeros(Np1,Np2);
  HiLo(:,:)= HiLo1(1:Np1,1:Np2); WF(:,:,:)= WF1(1:Np1,1:Np2);
    if kk==1
    imwrite(HiLo, hilo_file);
    imwrite(im2uint16(mat2gray(WF)), wf_file);
    
    elseif isstackwrite
%     imwrite(im2uint16(pSIM), pSIM_file,'WriteMode','append')

    imwrite(HiLo, hilo_file,'WriteMode','append');
    imwrite(im2uint16(mat2gray(WF)), wf_file,'WriteMode','append');
   
    fp1 = fopen(stacknum_file, 'w');
    fprintf(fp1,'%s',num2str(kk));
    fclose(fp1);
    else 
    imwrite(HiLo, hilo_file);
    imwrite(im2uint16(mat2gray(WF)), wf_file);
  
    end
end
clear;close all;clc;
alldoflag=1;
end

function [ out ] = hpgauss(H,W,SIGMA)
%   Creates a 2D Gaussian filter for a Fourier space image of height H and
%   width W. SIGMA is the standard deviation of the Gaussian.
out=1-lpgauss(H,W,SIGMA);
end

function [ out ] = lpgauss(H,W,SIGMA)
%   Creates a 2D Gaussian filter for a Fourier space image
%   W is the number of columns of the source image and H is the number of
%   rows. SIGMA is the standard deviation of the Gaussian.
H = single(H);
W = single(W);
kcx = (SIGMA);
kcy = ((H/W)*SIGMA);
temp0 = -floor(W/2);
[x,y] = meshgrid(-floor(W/2):floor((W-1)/2), -floor(H/2):floor((H-1)/2));
temp = -(x.^2/(kcx^2)+y.^2/(kcy^2));
out = ifftshift(exp(temp));
% out = ifftshift(exp(-(x.^2/(kcx^2)+y.^2/(kcy^2))));
end

function [ out ] = hp_OneSide(H,W,SIGMA)
%   Creates a 2D Gaussian filter for a Fourier space image of height H and
%   width W. SIGMA is the standard deviation of the Gaussian.

H = single(H);
W = single(W);
kcx = (SIGMA);
kcy = ((H/W)*SIGMA);
[x,y] = meshgrid(-floor(W/2):floor((W-1)/2), -floor(H/2):floor((H-1)/2));
lp_out = exp(-(x.^2/(kcx^2)+y.^2/(kcy^2)));

out=1-lp_out;

H_c = floor(H/2);
W_c = floor(W/2);

out(:, 1:W_c) = 0;
out(1: (H_c+1), W_c+1) = 0;
out = ifftshift(out);
end

function ld = psim_sim2d_recon1(raw_sim2d, theta, calib)
%% Calculate frequency components on orientational dimension
% obtain polarized wide field images
% [M,N,t]=size(raw_sim2d);
% ld=zeros(M,N,3);
% ld(:,:,1) = raw_sim2d(:,:,1)*1.5; 
% ld(:,:,2)  = raw_sim2d(:,:,2)*1.5./calib{1}; 
% ld(:,:,3) = raw_sim2d(:,:,3)*1.5./calib{2};

[M,N,t]=size(raw_sim2d);
ld=zeros(M,N,3);
ld(:,:,1) = raw_sim2d(:,:,1); 
ld(:,:,2)  = raw_sim2d(:,:,2); 
ld(:,:,3) = raw_sim2d(:,:,3);


% % cal fft
% d1_f = fft2(ld(:,:,1)); 
% d2_f = fft2(ld(:,:,2)); 
% d3_f = fft2(ld(:,:,3));
% 
% % calulate matrixes
% pol = mod(theta + pi/2,pi);
% mat_ld = 0.5*[1 0.5*exp(2i*pol(1)) 0.5*exp(-2i*pol(1));
%     1 0.5*exp(2i*pol(2)) 0.5*exp(-2i*pol(2));
%     1 0.5*exp(2i*pol(3)) 0.5*exp(-2i*pol(3))];
% mat_ld_inv = inv(mat_ld);
% % cal frequency components
% ld_fo = mat_ld_inv(1,1)*d1_f + mat_ld_inv(1,2)*d2_f + mat_ld_inv(1,3)*d3_f;%DC
% ld_fp = mat_ld_inv(2,1)*d1_f + mat_ld_inv(2,2)*d2_f + mat_ld_inv(2,3)*d3_f;%AC
% ld_fm = mat_ld_inv(3,1)*d1_f + mat_ld_inv(3,2)*d2_f + mat_ld_inv(3,3)*d3_f;
% % cal ld image
% ld_f = ld_fo;
% ld_f(:,:,2) = ld_fp; 
% ld_f(:,:,4) = ld_fm;
% 
% ld = abs(ifft(ifft(ifft(ld_f,[],1),[],2),[],3));

end

function [ouf, om, wf, cm, ld_ouf,h1] = recon_pm(pm, pol,cmin, cmax, calib)   %pm为空域ld；pol为theta角，cmin为ld最小~~~，wf为宽场图，om为pwf
%% test module
% pol = [0, pi/4, pi/4*2, pi/4*3];
% pm = cos(2*(pol-pi/4));

% pol = mod(pol + pi/2,pi);

% calibration
if (nargin < 6)
    calib = ones([size(pm, 1), size(pm, 2), 2]);
end

calib_1 = double(calib(:,:,1)) / 10000;
calib_2 = double(calib(:,:,2)) / 10000;
calib_1(calib_1 < 0.1) = 1;
calib_1(calib_1 > 10) = 1;
calib_2(calib_2 < 0.1) = 1;
calib_2(calib_2 > 10) = 1;

pm(:,:,2) = pm(:,:,2) ./ calib_1;
pm(:,:,3) = pm(:,:,3) ./ calib_2;

%%
if length(size(pm))<3
    pm = reshape(pm, 1, 1, length(pm));
end
%% cal matrix
mat_pm = zeros(length(pol),3);
for kk = 1 : size(mat_pm,1)
    mat_pm(kk,:) = [1 cos(2*pol(kk)) sin(2*pol(kk))];
end
mat_pm_inv = pinv(mat_pm);
% cal pm factors
pm_f = zeros(size(pm,1), size(pm,2), 3);
for kk = 1 : size(mat_pm_inv,1)
    for ll = 1 : size(mat_pm_inv,2)
        pm_f(:,:, kk)= pm_f(:,:,kk) + mat_pm_inv(kk,ll)*pm(:,:,ll);
    end
end
% obtain dc, ac
dc = pm_f(:,:,1);
ac = sqrt(pm_f(:,:,2).^2+pm_f(:,:,3).^2);

% obtain alpha
alpha = zeros(size(dc));
alpha1 = mod(atan(pm_f(:,:,3)./pm_f(:,:,2))/2, pi);
alpha2 = mod((atan(pm_f(:,:,3)./pm_f(:,:,2))+pi)/2, pi);
alpha(pm_f(:,:,2)>=0) = alpha1(pm_f(:,:,2)>=0);
alpha(pm_f(:,:,2)<0) = alpha2(pm_f(:,:,2)<0);

alpha = mod(alpha,pi);
% cal ouf
ouf = 2*ac./(ac+dc);
%% generate om image
%
h = alpha/pi;       % 将alpha转变为0~1
s = 0.6*ones(size(h));
% s = ouf;
%
wf = dc;
wf = max(wf-cmin, 0); 
wf = min(wf/cmax, 1);
v = wf;
%
om_hsv = cat(3, h, s, v);
om = hsv2rgb(om_hsv);
%% generate cm image
%
cm = zeros(100, 100);
s = 0.6*ones(size(cm));
%
xx = 0:size(cm,1)-1; yy = xx; c0 = size(cm,1)/2-0.5;
[xx,yy] = meshgrid(xx,yy);
radius = sqrt((xx-c0).^2+(yy-c0).^2);
mask = (radius <= 50) .* (radius>30);  
v = ones(size(cm)).*mask;
%
phy = atan((c0-yy)./(xx-c0));
phy = mod(phy, pi);
h = phy/pi;
%
cm_hsv = cat(3, h, s, v);
cm = hsv2rgb(cm_hsv);


%% 绘制OUF的伪彩图
h_ouf = ouf;
s_ouf = 0.6*ones(size(h_ouf));
v_ouf = wf;
ld_om_hsv = cat(3, h_ouf, s_ouf, v_ouf);
ld_ouf = hsv2rgb(ld_om_hsv);
% figure; imshow(ld_ouf,[])

 colorb = sort(h_ouf,'descend');
 [~,c2] = find(colorb == max(colorb(:)));
 maxb = colorb(:,c2);
 sizebar=size(maxb);
a1=1;a2=sizebar(:,1);a3=floor(a2/2);
 width = 50;
 h_col = repmat(maxb,1,width);
 s_col = 0.6*ones(size(h_col));
  v_col = ones(size(h_col));
ld_col_hsv = cat(3, h_col, s_col, v_col);
ld_col = hsv2rgb(ld_col_hsv);
h1 = figure('Visible', 'off'); 
set (gca,'position',[0.05,0.05,0.85,0.85] )
imshow(ld_col,[])
la = min(maxb); lb = max(maxb);
text(width+10,length(maxb),num2str(la))
text(width+10,1,num2str(lb))
yticks([]);
yticklabels({});
xticks([]);
xticklabels({});
frame = getframe(gcf);
h1 = frame2im(frame);
end
function raw_data = readMTiff(raw_file)
  info = imfinfo(raw_file);
  frames = numel(info);
  
  raw_data = zeros(info(1).Height, info(1).Width, frames, 'uint16');
  for k = 1:frames
      raw_data(:,:,k) = im2uint16(imread(raw_file, k));
  end
end