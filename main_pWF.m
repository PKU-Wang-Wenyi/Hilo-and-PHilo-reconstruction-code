clear all; close all; clc

%% Input
raw_file = 'pGC_RAW_C2_M01_Z001_T001.tif';
waveL = [525,0.00,5,0.9,0.05];

WF_file = '.\pWF\WF_001.tif';
pWF_file='.\pWF\pWF_001.tif';
pWF_ouf_file = '.\pWF\pWF_outf_001.tif';
pWF_ouf_pseudo_file = '.\pWF\pWF_outf_pseudocolor_001.tif';
pWF_ouf_pseudo_bar_file = '.\pWF\pWF_outf_pseudocolor_bar_001.tif';
cm_file = '.\pWF\cm_001.tif';
mkdir('.\pWF');
theta = deg2rad([-87.031690 , -27.004184 , 32.979669 ]);%rad
stacknum_file = '.\pWF\Stacknum.txt';
isstackwrite=1;
If_correction=0;
dir_corr= 'C:\program code\cpSIM\data\optical_section\20190122\correction.tif';
% 输出WF，pWF，OUF，OUF_color，colorbar各一张
tic
ReconPWF(waveL,theta,raw_file,WF_file, pWF_file, pWF_ouf_file, pWF_ouf_pseudo_file, pWF_ouf_pseudo_bar_file, cm_file,stacknum_file,isstackwrite,If_correction,dir_corr);
toc