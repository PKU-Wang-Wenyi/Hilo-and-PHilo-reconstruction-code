clear all; close all; clc

%% Input
raw_file = '12.tif';
waveL = [525,0.5,0.1,0.25,0.55,0.9,0.05];
hilo_file = '.\pHiLo\HiLo_001.tif';
wf_file = '.\pHiLo\WF_001.tif';
PWF_file='.\pHiLo\PWF_001.tif';
philo_file = '.\pHiLo\pHiLo_001.tif';
philo_ouf_file = '.\pHiLo\pHiLo_outf_001.tif';
philo_ouf_pseudo_file = '.\pHiLo\pHiLo_outf_pseudocolor_001.tif';
philo_ouf_pseudo_bar_file = '.\pHiLo\pHiLo_outf_pseudocolor_bar_001.tif';
cm_file = '.\pHiLo\cm_001.tif';
mkdir('.\pHiLo');
theta = deg2rad([-87.031690 , -27.004184 , 32.979669 ]);%rad
% theta =[-1.3900 0.7040 -0.3430]
stacknum_file = '.\pHiLo\Stacknum.txt';
isstackwrite=1;
If_correction=1;
dir_corr= '.\calib_sim.tif';
% 输出HiLo，pHiLo，OUF，OUF_color，colorbar各一张
% ratio1 = 0.25; ratio2 =0.55;
% hue1=0.9;
% hue2=0.05;
tic
[alldoflag]=ReconPHiLo(waveL,theta,raw_file, wf_file,PWF_file,hilo_file, philo_file, philo_ouf_file, philo_ouf_pseudo_file, philo_ouf_pseudo_bar_file, cm_file,stacknum_file,isstackwrite,If_correction,dir_corr);
toc