%%
addpath('~rouxf/MATLAB/utils/');
%%
p2d = '~rouxf/Data/EM/Neuralynx_data/other/2016-05-16_18-05-58/sort/5/';
fn = dir([p2d,'A46_sorted_new.mat']);

spdat = load([p2d,fn.name]);
%%
visualize_sorting(spdat);