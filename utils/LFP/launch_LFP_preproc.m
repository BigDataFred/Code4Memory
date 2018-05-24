function [CSC_preproc] = launch_LFP_preproc(dataset,trl)

CSC_preproc = cell(size(trl,1),1);

x = preproc_LFP(dataset,trl(1,:));
x = whos('x');
gb = round((x.bytes/10e5)/1000);

% Create a pool using default settings and disable the use of SPMD.
parpool('SpmdEnabled', false);
if gb < 2    
    parfor it = 1:size(trl,1)
        [CSC_preproc{it}] = preproc_LFP(dataset,trl(it,:));
    end;
else
    fprintf('Warning: datasize exceeds parallel-job capacity - running in serial mode\n');
    for it = 1:size(trl,1)
        [CSC_preproc{it}] = preproc_LFP(dataset,trl(it,:));
    end;
end;

[CSC_preproc] = ft_appenddata([],CSC_preproc{:});
delete(gcp);