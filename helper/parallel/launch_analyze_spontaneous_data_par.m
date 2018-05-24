function launch_analyze_spontaneous_data_par( pID, CSCidx, makeSpikes, makeLFP , makeMUA)
%%
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/custom/'));% custom code
addpath(genpath('/media/rouxf/rds-share/Common/releaseDec2015/'));
addpath(genpath('/media/rouxf/rds-share/Common/wave_clus/'));%needed for spike detection & sorting
addpath(genpath('/media/rouxf/rds-share/Common/osort-v3-rel/'));%needed for spike detection & sorting
fn = dir('/media/rouxf/rds-share/Common/fieldtrip-*');
addpath(['/media/rouxf/rds-share/Common/',fn.name]);% spectral analysis, signal processing, spike detection
ft_defaults;
%%
if nargin == 0
    pID = 'P06';
    CSCidx = 1:12;
    makeSpikes = 1;
    makeLFP = 1;
    makeMUA = 0;
end;

%%
savePath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/Spontaneous/']

chck = dir(savePath);
if isempty(chck)
    mkdir(savePath);
end;

%%
p2NLX = ['/media/rouxf/rds-share/Archive/MICRO/',pID,'/Spontaneous/'];

[NLXfiles] = dir([p2NLX,'2017-03-*']);

%%
for seshSel = length( NLXfiles )
    
    NLXsesh = NLXfiles(seshSel).name;
    
    %%
    chck = dir([savePath,NLXsesh, filesep, 'spike_dat']);
    if isempty(chck)
        mkdir([savePath,NLXsesh, filesep, 'spike_dat']);
    end;
    
    chck = dir([savePath,NLXsesh, filesep, 'lfp_dat']);
    if isempty(chck)
        mkdir([savePath,NLXsesh, filesep, 'lfp_dat']);
    end;
    
    chck = dir([savePath,NLXsesh, filesep, 'mua_dat']);
    if isempty(chck)
        mkdir([savePath,NLXsesh, filesep, 'mua_dat']);
    end;
    
    %% open pool of workers
    if isempty(gcp('nocreate'))
        parpool(length(CSCidx),'SpmdEnabled',false);
    end;
    
    %%
    parfor jt = 1:length(CSCidx)
        
        [data,readme,chan] = analyze_spontaneous_data(p2NLX,NLXsesh,CSCidx(jt),makeSpikes,makeLFP,makeMUA);
        
        if makeSpikes ==1
            chck = regexp(readme,'Spike');
            ix = [];
            for nt = 1:length(chck)
                ix(nt) = ~isempty(chck{nt});
            end;
            ix = find(ix);
            spike_data = data(ix);
            
            [outname] = [ pID , '_' , 'TS' , num2str(seshSel) , '_spike_data_', chan , '_' , NLXsesh , '.mat' ];
            par_save( [ savePath , NLXsesh,  filesep, 'spike_dat' , filesep, outname ] , [] , {spike_data});
        end;
        
        if makeLFP ==1
            chck = regexp(readme,'LFPdata');
            ix = [];
            for nt = 1:length(chck)
                ix(nt) = ~isempty(chck{nt});
            end;
            ix = find(ix);
            LFP_data = data(ix);
            
            [outname] = [ pID , '_' , 'TS' , num2str(seshSel) , '_lfp_data_', chan , '_' , NLXsesh , '.mat' ];
            par_save( [ savePath , NLXsesh,  filesep, 'lfp_dat' , filesep, outname ] , [] , {LFP_data});
        end;
        
        if makeMUA ==1
            chck = regexp(readme,'MUAdata');
            ix = [];
            for nt = 1:length(chck)
                ix(nt) = ~isempty(chck{nt});
            end;
            ix = find(ix);
            MUA_data = data(ix);
            
            [outname] = [ pID , '_' , 'TS' , num2str(seshSel) , '_mua_data_', chan , '_' , NLXsesh , '.mat' ];
            par_save( [ savePath , NLXsesh,  filesep, 'mua_dat' , filesep, outname ] , [] , {MUA_data});
        end;
        
    end;
    
    %%
    delete('dummy_CSC_data_tmpCSC_*.dg_01.lab');
    delete('dummy_CSC_data_tmpCSC_*.dg_01');
    delete('spc_log.txt');
    
end;
%%

delete(gcp);
exit;