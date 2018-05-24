fn = dir('/media/rouxf/rds-share/Common/fieldtrip-*');
addpath(['/media/rouxf/rds-share/Common/',fn.name]);% spectral analysis, signal processing, spike detection
ft_defaults;

%%
if isempty(gcp('nocreate'))
    parpool(36,'SpmdEnabled',false);
end;

%%
basepath2 = '/media/rouxf/rds-share/iEEG_DATA/MICRO/';
pIDs = dir(basepath2);
pIDs(1:2) = [];

%%
for yt = 2%5:length( pIDs )
    
    basepath2 = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pIDs(yt).name,filesep];
    nd = dir(basepath2);
    nd(1:2) = [];
    chck = regexp({nd.name},'EM','match');
    sel = [];for it = 1:length(chck); sel(it) = isempty(chck{it});end;
    basepath2 = [basepath2,nd(find(sel==0)).name,filesep];
    nd = dir(basepath2);
    nd(1:2) = [];
    nd = nd(find([nd.isdir]==1));
    
    %%
    tlck = cell(1,length( nd ));
    itpc = cell(1,length( nd ));
    
    pow1 = cell(1,length( nd ));
    pow2 = cell(1,length( nd ));
    pow3 = cell(1,length( nd ));
    pow4 = cell(1,length( nd ));
    
    tfr1 = cell(1,length( nd ));
    tfr2 = cell(1,length( nd ));
    tfr3 = cell(1,length( nd ));
    tfr4 = cell(1,length( nd ));
    
    %%
    for zt = 1:length( nd )
        
        p2d = [basepath2,nd(zt).name,filesep,'lfp_dat',filesep];
        files = dir([p2d,'*.mat']);
        
        %     %%
        %     chck = regexp({files(:).name},'CSC_L');
        %     nL = [];
        %     ixL = [];
        %     k = 0;
        %     for it = 1:length(chck)
        %
        %         if ~isempty(chck{it})
        %             k = k+1;
        %             ix2 = regexp(files(it).name(chck{it}:end),'_');
        %             ix2 = ix2(2);
        %             ixL(k) = it;
        %             nL(k) = str2double(files(it).name(chck{it}+5:chck{it}+ix2-2));
        %         end;
        %     end;
        %
        %     chck = regexp({files(:).name},'CSC_R');
        %     nR = [];
        %     ixR = [];
        %     k = 0;
        %     for it = 1:length(chck)
        %
        %         if ~isempty(chck{it})
        %             k = k+1;
        %             ix2 = regexp(files(it).name(chck{it}:end),'_');
        %             ix2 = ix2(2);
        %             ixR(k) = it;
        %             nR(k) = str2double(files(it).name(chck{it}+5:chck{it}+ix2-2));
        %         end;
        %     end;
        %
        %     [~,s_idx] = sort(nL);
        %     ixL = ixL(s_idx);
        %
        %     [~,s_idx] = sort(nR);
        %     ixR = ixR(s_idx);
        %
        %     files = files([ixL ixR]);
        
        %%
        if ~isempty( files )
            dat = cell( 1 , length( files ) );
            parfor it = 1:length(dat)
                
                dum = load([p2d,files(it).name]);
                dat{it} = dum.save_data{1}{1}{1};
                dat{it}
            end;
            
            %%
            cfg                     = [];
            
            [dat] = ft_appenddata( cfg , dat{:} );
            
            %%
            cfg                 = [];
            cfg.lpfilter        = 'yes';
            cfg.lpfreq          = 20;
            
            [dum]               = ft_preprocessing( cfg , dat );
            
            %%
            cfg                 = [];
            cfg.lpfilter        = 'yes';
            cfg.lpfreq          = 20;
            cfg.hilbert         = 'complex';
            
            [dum2] = ft_preprocessing( cfg , dat );
            
            %%
            cfg                 = [];
            cfg.latency         = [-1 3];
            
            [dum]               = ft_selectdata( cfg , dum );
            
            cfg                 = [];
            cfg.latency         = [-1 3];
            
            [dum2]               = ft_selectdata( cfg , dum2 );
            
            %%
            cfg                     = [];
            cfg.keeptrials          = 'no';
            
            [tlck{zt}] = ft_timelockanalysis( cfg , dum );
            clear dum;
            
            %%
            cfg                     = [];
            cfg.keeptrials          = 'yes';
            
            [tlck2] = ft_timelockanalysis( cfg , dum2 );
            clear dum2;
            
            %%
            N = size(tlck2.trial,1);
            
            itpc{zt} = tlck2.trial./abs( tlck2.trial );
            itpc{zt} = sum( itpc{zt} , 1 );
            itpc{zt} = abs( itpc{zt} )./N;
            itpc{zt} = squeeze( itpc{zt} );
            clear tlck2;
            
            %%
            pf = 0.05;
            pf = (1/pf-1);
            
            if pf ==0
                pf=1;
            end;
            
            cfg                     = [];
            cfg.latency             = [-1 -1e-3];
            
            [dum1] = ft_selectdata( cfg , dat );
            
            cfg                     = [];
            cfg.latency             = [2 3];
            
            [dum2] = ft_selectdata( cfg , dat );
            
            cfg                     = [];
            cfg.detrend             = 'yes';
            cfg.demean              = 'yes';
            
            [dum1] = ft_preprocessing( cfg , dum1 );
            
            cfg                     = [];
            cfg.detrend             = 'yes';
            cfg.demean              = 'yes';
            
            [dum2] = ft_preprocessing( cfg , dum2 );
                                    
            zp = zeros(length(dum1.label),length(dum1.time{1})+length(dum1.time{1})*ceil(pf/2));
            
            for it = 1:length( dum1.trial )
                dum1.trial{it} = [zp dum1.trial{it} zp];
                dum1.time{it} = 0:1/dum1.fsample:(length(dum1.trial{it})-1)/dum1.fsample;
            end;
            
            zp = zeros(length(dum2.label),length(dum2.time{1})+length(dum2.time{1})*ceil(pf/2));
            
            for it = 1:length( dum2.trial )
                dum2.trial{it} = [zp dum2.trial{it} zp];
                dum2.time{it} = 0:1/dum2.fsample:(length(dum2.trial{it})-1)/dum2.fsample;
            end;
            
            w = length(dum1.time{1})/dum1.fsample;
            
            cfg                     = [];
            cfg.method              = 'mtmfft';
            cfg.pad                 = 'maxperlen';
            cfg.taper               = 'dpss';
            cfg.tapsmofrq           = 1/w;
            cfg.output              = 'fourier';
            
            [freq1] = ft_freqanalysis( cfg , dum1 );
            
            [freq2] = ft_freqanalysis( cfg , dum2 );
            clear dum*;
            
            cfg                     = [];
            cfg.frequency           = [1 20];
            
            [freq1] = ft_selectdata( cfg ,  freq1 );
            [freq2] = ft_selectdata( cfg , freq2 );
            
            dim = size(freq1.fourierspctrm);
            
            cfg                     = [];
            cfg.keeptrials          = 'yes';
            cfg.jackknife           = 'yes';
            cfg.variance            = 'yes';
            
            [pow1{zt}] = ft_freqdescriptives( cfg , freq1 );
            [pow2{zt}] = ft_freqdescriptives( cfg , freq2 );                                   
            
            %%
            cfg                     = [];
            cfg.latency             = [-2 5];
            
            [dum1] = ft_selectdata( cfg , dat );
            
            cfg                     = [];
            cfg.detrend             = 'yes';
            cfg.demean              = 'yes';
            
            [dum1] = ft_preprocessing( cfg , dum1 );
            
            cfg                     = [];
            cfg.method              = 'mtmconvol';
            cfg.pad                 = 'maxperlen';
            cfg.taper               = 'dpss';
            cfg.tapsmofrq           = 1;
            cfg.output              = 'fourier';
            cfg.foi                 = 1:20;
            cfg.t_ftimwin           = ones(1,length(cfg.foi));
            cfg.toi                 = dum1.time{1}(1):0.01:dum1.time{1}(end);
            cfg.keeptrials          = 'yes';
            cfg.polyremoval         = -1;
            
            [freq1] = ft_freqanalysis( cfg , dum1 );
            clear dum*;
                        
            cfg                     = [];
            cfg.keeptrials          = 'yes';
            cfg.jacknife            = 'yes';
            cfg.variance            = 'yes';
            
            [tfr1{zt}] = ft_freqdescriptives( cfg , freq1 );
            clear freq1;
            
            cfg                     = [];
            cfg.latency             = [-1 4];
            cfg.frequency           = [1 20];
            
            [tfr1{zt}] = ft_selectdata( cfg , tfr1{zt} );               
            
            %%
            pf = 1;
            pf = (1/pf-1);
            
            if pf ==0
                pf=1;
            end;
            
            cfg                     = [];
            cfg.latency             = [-1 -1e-3];
            
            [dum1] = ft_selectdata( cfg , dat );
            
            cfg                     = [];
            cfg.latency             = [2 3];
            
            [dum2] = ft_selectdata( cfg , dat );
            
            cfg                     = [];
            cfg.detrend             = 'yes';
            cfg.demean              = 'yes';
            
            [dum1] = ft_preprocessing( cfg , dum1 );
            
            cfg                     = [];
            cfg.detrend             = 'yes';
            cfg.demean              = 'yes';
            
            [dum2] = ft_preprocessing( cfg , dum2 );
                                    
            zp = zeros(length(dum1.label),length(dum1.time{1})+length(dum1.time{1})*ceil(pf/2));
            
            for it = 1:length( dum1.trial )
                dum1.trial{it} = [zp dum1.trial{it} zp];
                dum1.time{it} = 0:1/dum1.fsample:(length(dum1.trial{it})-1)/dum1.fsample;
            end;
            
            zp = zeros(length(dum2.label),length(dum2.time{1})+length(dum2.time{1})*ceil(pf/2));
            
            for it = 1:length( dum2.trial )
                dum2.trial{it} = [zp dum2.trial{it} zp];
                dum2.time{it} = 0:1/dum2.fsample:(length(dum2.trial{it})-1)/dum2.fsample;
            end;
            
            w = length(dum1.time{1})/dum1.fsample;
            
            cfg                     = [];
            cfg.method              = 'mtmfft';
            cfg.pad                 = 'maxperlen';
            cfg.taper               = 'dpss';
            cfg.tapsmofrq           = 10;
            cfg.output              = 'fourier';
            
            [freq1] = ft_freqanalysis( cfg , dum1 );
            
            [freq2] = ft_freqanalysis( cfg , dum2 );
            clear dum*;
            
            cfg                     = [];
            cfg.frequency           = [20 200];
            
            [freq1] = ft_selectdata( cfg ,  freq1 );
            [freq2] = ft_selectdata( cfg , freq2 );
                        
            cfg                     = [];
            cfg.keeptrials          = 'yes';
            cfg.jackknife           = 'yes';
            cfg.variance            = 'yes';
                        
            [pow3{zt}] = ft_freqdescriptives( cfg , freq1 );
            [pow4{zt}] = ft_freqdescriptives( cfg , freq2 );            
            
            %%
            cfg                     = [];
            cfg.latency             = [-2 5];
            
            [dum1] = ft_selectdata( cfg , dat );
            
            cfg                     = [];
            cfg.detrend             = 'yes';
            cfg.demean              = 'yes';            
            
            [dum1] = ft_preprocessing( cfg , dum1 );
            
            cfg                     = [];
            cfg.method              = 'mtmconvol';
            cfg.pad                 = 'maxperlen';
            cfg.taper               = 'dpss';
            cfg.tapsmofrq           = 10;
            cfg.output              = 'fourier';
            cfg.foi                 = 20:200;
            cfg.t_ftimwin           = 0.25*ones(1,length(cfg.foi));
            cfg.toi                 = dum1.time{1}(1):0.01:dum1.time{1}(end);
            cfg.keeptrials          = 'yes';
            cfg.polyremoval         = -1;
            
            [freq1] = ft_freqanalysis( cfg , dum1 );
            clear dum*;
                        
            cfg                     = [];
            cfg.keeptrials          = 'yes';
            cfg.jacknife            = 'yes';
            cfg.variance            = 'yes';
            
            [tfr2{zt}] = ft_freqdescriptives( cfg , freq1 );
            clear freq1;
            
            cfg                     = [];
            cfg.latency             = [-1 4];
            cfg.frequency           = [20 200];
            
            [tfr2{zt}] = ft_selectdata( cfg , tfr2{zt} );                          
            
        end;
        
    end;
    
    %%
    readme = {'baseline','poststim'};
    
    save([basepath2,'event_related_data_',pIDs(yt).name,'.mat'],'tlck','itpc');
    save([basepath2,'low_freq_spectrum_',pIDs(yt).name,'.mat'],'pow1','pow2','readme');
    save([basepath2,'low_freq_tfr_',pIDs(yt).name,'.mat'],'tfr1');
    save([basepath2,'high_freq_spectrum_',pIDs(yt).name,'.mat'],'pow3','pow4','readme');
    save([basepath2,'high_freq_tfr_',pIDs(yt).name,'.mat'],'tfr2');
    
end;

%%
exit;
delete(gcp);
