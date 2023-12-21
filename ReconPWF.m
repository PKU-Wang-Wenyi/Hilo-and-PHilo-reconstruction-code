function ReconPWF(waveL,theta,raw_file,WF_file, pWF_file, pWF_ouf_file, pWF_ouf_pseudo_file, pWF_ouf_pseudo_bar_file, cm_file,stacknum_file,isstackwrite,If_correction,dir_corr)
%     function ReconPHiLo(waveL,theta,raw_file,wf_file,PWF_file, hilo_file, philo_file, philo_ouf_file, philo_ouf_pseudo_file, philo_ouf_pseudo_bar_file, cm_file,stacknum_file,isstackwrite,If_correction,dir_corr)
N = 6;
ratio1 = waveL(2);
ratio2 =waveL(3);
hue1=waveL(4);
hue2=waveL(5);
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


%% iter
for kk =  1 : nz
    raw_sim = zeros(NPixel,NPixel,N);

    for ll = 1 : N
        raw_sim(:,:,ll) = double((raw_data(:,:,(kk-1)*N+ll)));
        %         sum = sum +  raw_sim(:,:,ll);
    end


    WF = wf_display1(raw_sim, calib);
    ld = psim_sim2d_reconWF(raw_sim, theta, calib);
    [ouf_pHiLo, pWF, ~, cm,ld_ouf,h1] = recon_pm(ld, theta,min(ld(:)), max(ld(:)), calib,ratio1,ratio2,hue1,hue2);
    %% saving

    if kk==1

        imwrite(im2uint16(mat2gray(WF)),WF_file);

        imwrite(uint16(pWF*65535), pWF_file);
        imwrite(uint16(ouf_pHiLo*65535), pWF_ouf_file);
        imwrite(uint16(cm*65535), cm_file);
        imwrite(uint16(ld_ouf*65535), pWF_ouf_pseudo_file);
        imwrite(h1, pWF_ouf_pseudo_bar_file);
    elseif isstackwrite
        %     imwrite(im2uint16(pSIM), pSIM_file,'WriteMode','append')


        imwrite(im2uint16(mat2gray(WF)),WF_file,'WriteMode','append');

        imwrite(uint16(pWF*65535), pWF_file,'WriteMode','append');
        imwrite(uint16(ouf_pHiLo*65535), pWF_ouf_file,'WriteMode','append');
        imwrite(uint16(cm*65535), cm_file,'WriteMode','append');
        imwrite(uint16(ld_ouf*65535), pWF_ouf_pseudo_file,'WriteMode','append');
        imwrite(h1, pWF_ouf_pseudo_bar_file,'WriteMode','append');
        fp1 = fopen(stacknum_file, 'w');
        fprintf(fp1,'%s',num2str(kk));
        fclose(fp1);
    else

        imwrite(im2uint8(mat2gray(WF)),WF_file);

        imwrite(uint8(pWF*255), pWF_file);
        imwrite(uint8(ouf_pHiLo*255), pWF_ouf_file);
        imwrite(uint8(cm*255), cm_file);
        imwrite(uint8(ld_ouf*255), pWF_ouf_pseudo_file);
        imwrite(h1, pWF_ouf_pseudo_bar_file);
    end
end
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

function [ouf1, om, wf, cm, ld_ouf,h1] = recon_pm(pm, pol,cmin, cmax, calib,ratio1,ratio2,hue1,hue2)   %pm为空域ld；pol为theta角，cmin为ld最小~~~，wf为宽场图，om为pwf
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

ouf1=ouf;
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