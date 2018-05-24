%%
addpath('~rouxf/MATLAB/toolboxes/fieldtrip-20160309/');
ft_defaults;
%%
p2d ='/home/rouxf/Data/EM/Neuralynx_data/other/2016-05-17_11-59-29/A/preproc/';
files = dir([p2d,'spike_triggeredAVG_cluster_*_part*.mat']);
%%
ID = cell(length(files),1);
for it = 1:length(files)

	ix(1) = regexp(files(it).name,'cluster_')+8;
	ix(2) = max(regexp(files(it).name,'_'))-1;

	ID{it} = files(it).name(ix(1):ix(2));
end;

ID = unique(ID);
 
%%
for it = 1:length(ID)
	fn = dir([p2d,'spike_triggeredAVG_cluster_',ID{it},'_part*.mat']);

	tmp = cell(1,length(fn));
	for jt = 1:length(fn)

		load([p2d,fn(jt).name]);
		k = fn(jt).name(regexp(fn(jt).name,'part')+4);
		tmp{jt} = eval(['dsCSC',k]);

	end;
	
	dat = ft_appenddata([],tmp{:});
	clear tmp;

	sn = [fn(jt).name(1:regexp(fn(jt).name,'part')-1),'concatenated.mat'];

	save([p2d,sn],'dat');
	clear dat;

end;
