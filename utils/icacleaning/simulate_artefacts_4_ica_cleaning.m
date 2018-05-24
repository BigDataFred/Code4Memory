%%
restoredefaultpath;
addpath('/bcbl/home/home_a-f/froux/fieldtrip-20140527/');
ft_defaults;
%%
path2file1 = '/bcbl/home/home_a-f/froux/MEG/pilot04_vslmeg_adapted2/';
file1 = 'pilot4_preproc_pre_ICA.mat';

path2file2 = '/bcbl/home/home_a-f/froux/MEG/pilot22_vslmeg_adapted2/';
file2 = 'Pilot22_preproc_pre_ICA.mat';
%%
dat1 = load([path2file1,file1]);
dat2 = load([path2file2,file2]);