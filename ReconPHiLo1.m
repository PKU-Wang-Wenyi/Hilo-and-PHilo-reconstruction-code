function [alldoflag]=ReconPHiLo1(waveL,theta,raw_file,wf_file,PWF_file, hilo_file, philo_file, philo_ouf_file, philo_ouf_pseudo_file, philo_ouf_pseudo_bar_file, cm_file,stacknum_file,isstackwrite,If_correction,dir_corr)
N = 6;alldoflag=0;disp('1!');
input_parameters=[waveL(2),waveL(3)];
% ratio1 = 0.25;
% ratio2 =0.55;
% hue1=0.9;
% hue2=0.05;
theta
ratio1 = waveL(4);
ratio2 =waveL(5);
hue1=waveL(6);
hue2=waveL(7);
waveL=waveL(1);
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
calibrationFlag = If_correction;
if calibrationFlag == 1
    calib=double(readMTiff(dir_corr));
    calib10 = calib(:,:,1);
    calib20 = calib(:,:,2);
    calib1 = zeros(NPixel,NPixel);
    calib2 = zeros(NPixel,NPixel);
    if Np1 == Np2 || Np2 > Np1
        calib1(1:Np1,:)=calib10(1:Np1,:);
        calib2(1:Np1,:)=calib20(1:Np1,:);
    else
        calib1(:,1:Np1)=calib10(:,1:Np1);
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
NA = 1.49;                      % objective NA
% waveL = 620;                    % wavelength(nm) of the emission light
pixel_size = 65;                % pixel size(nm) corresponding to the sample
M_pattern = 0.6;                % sin pattern contrast, for calculating scaling factor
Lo_scaling = 1;
% Lateral_resolution = 15.385;     % Lateral resolution in px/um
If_scaling_own = 0;             % whether to select the scaling factor manually
Lo_scaling_own = 0.4;           % manually selected scaling factor
Lo_scaling_own = 0.2;           % manually selected scaling factor

%% iter
for kk =  1 : nz
    raw_sim = zeros(NPixel,NPixel,6);

    for ll = 1 : 6
        %         raw_sim(:,:,ll) =raw_data(:,:,(kk-1)*6+ll);
        raw_sim(1:Np1,1:Np2,ll) = raw_data(:,:,(kk-1)*6+ll);
        %         sum = sum +  raw_sim(:,:,ll);
    end
    [M,N]=size(raw_sim);
    uniform_file = zeros(M,N/6,3);
    for ii = 1:3
        uniform_file(:,:,ii) = (raw_sim(:,:,2*ii-1) + raw_sim(:,:,2*ii))./2;
        img = uniform_file(:,:,ii);
        %% Define parameters
        % res = Lateral_resolution*1e6*0.0254;    % Lateral resolution
        % uniform_file = (raw_data{1}{1,1} + raw_data{1}{2,1} + raw_data{1}{3,1} + raw_data{1}{4,1}+ raw_data{1}{5,1} + raw_data{1}{6,1})/6;
        [h, w] = size(img);
        res = 0.61 * waveL / (NA);                   % resolution
        k_m=w/(10*input_parameters(1));
        kc = nearest(k_m * 0.2);                % cut-off frequency between hp and lp filter
        k_side_HP = kc/2;                     % cutoff frequency of high-pass filter for differential image
        sigmaLP = kc*2/2.355;                   % Finding sigma value for low pass
        lambda = nearest(w/(2*kc));             % Side length of contrast
        % evalutation mask by rounding to
        % nearest odd integer
        if mod(lambda,2) == 0                   % lambda must be odd
            lambda = lambda+1;
        else
        end
        if If_correction == 0                   % Check if correction file is selected
            gauss_corr = 1;
        else
            gauss_corr = double(imread(dir_corr));	% Create gaussian field illumination correction matrix if file was selected
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
        u = double(uniform_file(:,:,ii));
        uni = padarray(u./gauss_corr,[lambda lambda],'symmetric');
        s = double(raw_sim(:,:,2*ii-1));
        rat = padarray(s./gauss_corr,[lambda lambda],'symmetric');
        %         rat = rat./uni;
        %% One-side highpass filtering ratio image stack and Contrast evaluation
        rat_hp_s_f = fft2(rat).*hp_side;
        rat_hp_s_f = rat_hp_s_f .* hp_1;
        rat_hp_s = ifft2(rat_hp_s_f);
        weight = sqrt(rat_hp_s .* conj(rat_hp_s));
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
        % % method one
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
    HiLo_sum = HiLo_dir(:,:,1) + HiLo_dir(:,:,2) + HiLo_dir(:,:,3);
    HiLo_sum = HiLo_sum/ max(HiLo_sum(:));
    HiLo = (HiLo_sum);


    HiLo_dir_uint = im2uint16(HiLo_dir);
    ld = psim_sim2d_recon1(HiLo_dir_uint, theta, calib);
    ld1 = psim_sim2d_reconWF(raw_sim, theta, calib);
    %[ouf_pHiLo, pHiLo, wf, cm, ld_ouf, h1] = recon_pm(ld, theta, min(ld(:)), max(ld(:)), calib);
    [ouf_pHiLo, pHiLo, ~, cm, ld_ouf, h1] = recon_pm(ld,ld1, theta, min(ld(:)), max(ld(:)), calib,ratio1,ratio2,hue1,hue2);

    % WF = wf_display1(raw_sim, calib);
    WF=mean(uniform_file,3);
    aa=double(max(WF(:)));
    HiLo=HiLo.*aa;
    ld = psim_sim2d_reconWF(raw_sim, theta, calib);
    [~, pWF, ~, ~, ~, ~] = recon_pm(ld, ld1,theta,min(ld(:)), max(ld(:)), calib,ratio1,ratio2,hue1,hue2);
    %% saving
    WF1=WF;
    HiLo1=HiLo;pWF1=pWF;ouf_pHiLo1=ouf_pHiLo;ld_ouf1=ld_ouf;pHiLo1=pHiLo;
    clear WF HiLo pWF ouf_pHiLo ld_ouf pHiLo
    HiLo=zeros(Np1,Np2); WF=zeros(Np1,Np2);pWF=zeros(Np1,Np2,3);ouf_pHiLo=zeros(Np1,Np2);ld_ouf=zeros(Np1,Np2,3);pHiLo=zeros(Np1,Np2,3);
    HiLo(:,:)= HiLo1(1:Np1,1:Np2); WF(:,:,:)= WF1(1:Np1,1:Np2);pWF(:,:,:)= pWF1(1:Np1,1:Np2,:);ouf_pHiLo(:,:,:)= ouf_pHiLo1(1:Np1,1:Np2);ld_ouf(:,:,:)= ld_ouf1(1:Np1,1:Np2,:);pHiLo(:,:,:)= pHiLo1(1:Np1,1:Np2,:);
    if kk==1
        imwrite(uint16(HiLo), hilo_file);
        imwrite(im2uint16((WF/65535)), wf_file);
        imwrite(uint8(pHiLo*255), philo_file);
        imwrite(uint8(pWF*255), PWF_file);
        imwrite(uint8(ouf_pHiLo*255), philo_ouf_file);
        imwrite(uint8(cm*255), cm_file);
        imwrite(uint8(ld_ouf*255), philo_ouf_pseudo_file);
        imwrite(h1, philo_ouf_pseudo_bar_file);
        fp1 = fopen(stacknum_file, 'w');
        fprintf(fp1,'Author:WenyiWang\n');
        fprintf(fp1,'processed num %s data',num2str(kk));
        fclose(fp1);
    else
        %     imwrite(im2uint16(pSIM), pSIM_file,'WriteMode','append')

        imwrite(uint16(HiLo), hilo_file,'WriteMode','append');
        imwrite(im2uint16(mat2gray(WF)), wf_file,'WriteMode','append');
        imwrite(uint8(pHiLo*255), philo_file,'WriteMode','append');
        imwrite(uint8(pWF*255), PWF_file,'WriteMode','append');
        imwrite(uint8(ouf_pHiLo*255), philo_ouf_file,'WriteMode','append');
        imwrite(uint8(cm*255), cm_file,'WriteMode','append');
        imwrite(uint8(ld_ouf*255), philo_ouf_pseudo_file,'WriteMode','append');
        imwrite(h1, philo_ouf_pseudo_bar_file,'WriteMode','append');
        fp1 = fopen(stacknum_file, 'w');
        fprintf(fp1,'Author:WenyiWang\n');
        fprintf(fp1,'processed num %s data',num2str(kk));
        fclose(fp1);
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
function ld = psim_sim2d_recon1(raw_sim2d, theta, calib)
%% Calculate frequency components on orientational dimension
% obtain polarized wide field images
% [M,N,t]=size(raw_sim2d);
% ld=zeros(M,N,3);
% ld(:,:,1) = mean(raw_sim2d(:,:,1:2),3)*1.5;
% ld(:,:,2)  = mean(raw_sim2d(:,:,3:4),3)*1.5./calib{1};
% ld(:,:,3) = mean(raw_sim2d(:,:,5:6),3)*1.5./calib{2};


[M,N,t]=size(raw_sim2d);
ld=zeros(M,N,3);
ld(:,:,1) = mean(raw_sim2d(:,:,1),3);
ld(:,:,2)  = mean(raw_sim2d(:,:,2),3);
ld(:,:,3) = mean(raw_sim2d(:,:,3),3);


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
function WF = wf_display1(raw_sim2d, calib)
%% Calculate frequency components on orientational dimension
% obtain polarized wide field images
[M,N,t]=size(raw_sim2d);
ld_wf=zeros(M,N,3);
ld_wf(:,:,1) = mean(raw_sim2d(:,:,1),3)*1.5;
ld_wf(:,:,2)  = mean(raw_sim2d(:,:,2),3)*1.5./calib(:,:,1);
ld_wf(:,:,3) = mean(raw_sim2d(:,:,3),3)*1.5./calib(:,:,2);
WF = mean(ld_wf,3);
WF = uint16(WF);
end
function [ out ] = lpgauss(H,W,SIGMA)
%   Creates a 2D Gaussian filter for a Fourier space image
%   W is the number of columns of the source image and H is the number of
%   rows. SIGMA is the standard deviation of the Gaussian.
H = double(H);
W = double(W);
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

H = double(H);
W = double(W);
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

function ld = psim_sim2d_reconWF(raw_sim2d, theta, calib)
%% Calculate frequency components on orientational dimension
% obtain polarized wide field images
% [M,N,t]=size(raw_sim2d);
% ld=zeros(M,N,3);
% ld(:,:,1) = raw_sim2d(:,:,1)*1.5;
% ld(:,:,2)  = raw_sim2d(:,:,2)*1.5./calib{1};
% ld(:,:,3) = raw_sim2d(:,:,3)*1.5./calib{2};

[M,N,t]=size(raw_sim2d);
ld=zeros(M,N,3);
ld(:,:,1) = mean(raw_sim2d(:,:,1:2),3);
ld(:,:,2)  = mean(raw_sim2d(:,:,3:4),3);
ld(:,:,3) = mean(raw_sim2d(:,:,5:6),3);


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

function [ouf, om, wf, cm, ld_ouf,h1] = recon_pm(pm,pm1, pol,cmin, cmax, calib,ratio1,ratio2,hue1,hue2)   %pm为空域ld；pol为theta角，cmin为ld最小~~~，wf为宽场图，om为pwf
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

% pm(:,:,2) = pm(:,:,2) ./ calib_1;
% pm(:,:,3) = pm(:,:,3) ./ calib_2;
pm1(:,:,2) = pm1(:,:,2) ./ calib_1;
pm1(:,:,3) = pm1(:,:,3) ./ calib_2;
%%
if length(size(pm))<3
    pm = reshape(pm, 1, 1, length(pm));
end
if length(size(pm1))<3
    pm1 = reshape(pm1, 1, 1, length(pm1));
end
%% cal matrix
mat_pm = zeros(length(pol),3);
for kk = 1 : size(mat_pm,1)
    mat_pm(kk,:) = [1 cos(2*pol(kk)) sin(2*pol(kk))];
end
mat_pm_inv = pinv(mat_pm);
% cal pm factors
pm_f = zeros(size(pm1,1), size(pm1,2), 3);
for kk = 1 : size(mat_pm_inv,1)
    for ll = 1 : size(mat_pm_inv,2)
        pm_f(:,:, kk)= pm_f(:,:,kk) + mat_pm_inv(kk,ll)*pm1(:,:,ll);
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
wf = mean(pm,3);
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
ouf(isnan(ouf))=0;
ouf(ouf>ratio2) = ratio2;
% OUFP245A(OUFP245A>ratio2) = 0;
ouf(ouf<ratio1) = ratio1;
ouf= interp1([ratio1, ratio2], [hue1, hue2],ouf);
ouf(isnan(ouf))=0;
h_ouf = ouf;
s_ouf = 0.6*ones(size(h_ouf));
v_ouf = wf;
ld_om_hsv = cat(3, h_ouf, s_ouf, v_ouf);
ld_ouf = hsv2rgb(ld_om_hsv);
% figure; imshow(ld_ouf,[])

colorb = sort(h_ouf,'descend');
[~,c2] = find(colorb == max(colorb(:)));
maxb = colorb(:,c2(1));
maxb=linspace(0,1,size(h_ouf,1)).';
sizebar=size(maxb);
a1=1;a2=sizebar(:,1);a3=floor(a2/2);
width = 50;
% h_col = repmat(maxb,1,50);
h_col = repmat(maxb,1,width);
h_col = interp1([0, 1], [hue1, hue2], h_col);
s_col = 0.6*ones(size(h_col));
v_col = ones(size(h_col));
ld_col_hsv = cat(3, h_col, s_col, v_col);
ld_col = hsv2rgb(ld_col_hsv);
h1 = figure('Visible', 'off');
set (gca,'position',[0.05,0.05,0.85,0.85] )
imshow(ld_col,[])
la = ratio1; lb = ratio2;
text(width+10,length(maxb),num2str(la))
text(width+10,1,num2str(lb))
yticks([]);
yticklabels({});
xticks([]);
xticklabels({});
frame = getframe(gcf);
h1 = frame2im(frame);
end
function [om] = recon_WF1(pm, pol,cmin, cmax, calib)   %pm为空域ld；pol为theta角，cmin为ld最小~~~，wf为宽场图，om为pwf
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

end
function raw_data = readMTiff(raw_file)
info = imfinfo(raw_file);
frames = numel(info);

raw_data = zeros(info(1).Height, info(1).Width, frames, 'uint16');
for k = 1:frames
    raw_data(:,:,k) = im2uint16(imread(raw_file, k));
end
end