clear all; close all; clc

%% Input
raw_file = '12.tif';
% raw_file = 'RAW.tif';
waveL = [525,0.35,0.5]; 
waveL = [525,0.95,0.5]; 
waveL = [525,0.35,0.15];
hilo_file = '.\HiLo\HiLo_001.tif';
wf_file = '.\HiLo\WF_001.tif';
mkdir('.\HiLo');
stacknum_file = '.\HiLo\Stacknum.txt';
isstackwrite=1;
% 输出HiLo，pHiLo，OUF，OUF_color，colorbar各一张
input_parameters=[0.35,0.5];
tic
[alldoflag]=ReconHiLo(waveL,raw_file, wf_file,hilo_file,stacknum_file,isstackwrite);
toc