function par=set_parameters_Bham(sr,varargin)

par = [];
par.fname_in = 'tmp';
par.fnamespc = ['CSC_data_tmp'];
% SPC PARAMETERS 
par.mintemp = 0.00;                  % minimum temperature for SPC 
par.maxtemp = 0.201;                 % maximum temperature for SPC 
par.tempstep = 0.01;                 % temperature steps 
par.num_temp = floor(...
    (par.maxtemp - par.mintemp)/par.tempstep);
par.stab = 0.8;
par.SWCycles = 100;                  % SPC iterations for each temperature (default 100) 
par.KNearNeighb = 11;                % number of nearest neighbors for SPC 
par.randomseed = 0;                  % if 0, random seed is taken as the clock value (default 0) 
%par.randomseed = 147;               % If not 0, random seed 
par.min_clus_abs = 60;               % minimum size of a cluster (default 60) 
%par.max_clus = 33;                  % maximum number of clusters allowed (default 13) 
par.min_clus_rel = 0.005;            % minimum cluster size, relative to the total nr. of spikes (only for batch scripts). 
par.foce_auto  = 'y';             % automatically force membership if temp >3
%par.temp_plot = 'lin';              % temperature plot in linear scale 
%par.temp_plot = 'log';              % temperature plot in log scale 

% DETECTION PARAMETERS
par.tmax= 'all';                     % maximum time to load (default)
%par.tmax= 180;                       % maximum time to load (in sec)
par.tmin= 0;                         % starting time for loading (in sec)
par.w_pre=20;                        % number of pre-event data points stored
par.w_post=44; 60;%                      % number of post-event data points stored
%par.detection = 'pos';               % type of threshold (default 'pos')
%par.detection = 'neg';
par.detection = 'both';
if isempty(varargin)
    par.stdmin = 5;                      % minimum threshold for detection
else
    par.stdmin = varargin{1};
end;
par.stdmax = 50;                     % maximum threshold for detection
par.alignment_window = 10;           % number of points around the sample expected to be the maximum 
par.ref_ms = 1.5;                      % detector dead time, minimum refractory period (in ms) 

% INTERPOLATION PARAMETERS
par.interpolation = 'y';             % interpolation with cubic splines (default)
%par.interpolation = 'n';
par.int_factor = 4;                  % interpolation factor
%par.int_factor = 2;                  % interpolation factor

par.detect_fmin = 500;               %high pass filter for detection (default 300)
par.detect_fmax = 1500;              %low pass filter for detection (default 1000)
par.detect_order = 2;
par.sort_order = 2;
par.sort_fmin = 500;                 %high pass filter for sorting (default 300)
par.sort_fmax = 5000;                %low pass filter for sorting (default 3000)

par.segments_length = 6;             %length of segments in which the data is cutted (default 5min).
par.sr=sr;                           % sampling rate

ref = 1.5;                           % detector dead time (in ms)
par.ref = floor(ref *sr/1000);       % conversion to datapoints
par.segments = 6;                    %number of segments in which data is cutted

% FEATURES PARAMETERS 
par.inputs = 10;                       % number of inputs to the clustering 
par.scales = 4;                        % number of scales for the wavelet decomposition 
par.features = 'wav';                % type of feature ('wav' or 'pca') 
if strcmp(par.features,'pca')
    par.inputs = 3;
end;
 
% FORCE MEMBERSHIP PARAMETERS 
par.template_sdnum = 3;             % max radius of cluster in std devs. 
par.template_k = 10;                % # of nearest neighbors 
par.template_k_min = 10;            % min # of nn for vote 
%par.template_type = 'mahal';       % nn, center, ml, mahal 
par.template_type = 'center';       % nn, center, ml, mahal 
par.force_feature = 'spk';          % feature use for forcing (whole spike shape) 
%par.force_feature = 'wav';         % feature use for forcing (wavelet coefficients). 
par.force_auto = true;              %automatically force membership (only for batch scripts). 
 
% TEMPLATE MATCHING 
par.match = 'y';                    % for template matching 
%par.match = 'n';                   % for no template matching 
par.max_spk = 20000;                % max. # of spikes before starting templ. match. 
par.permut = 'y';                   % for selection of random 'par.max_spk' spikes before starting templ. match. 
% par.permut = 'n';                 % for selection of the first 'par.max_spk' spikes before starting templ. match. 
 
% HISTOGRAM PARAMETERS 
par.nbins = 100;                    % # of bins for the ISI histograms 
par.bin_step = 1;                   % percentage number of bins to plot 
