function [params] = set_up_stimuli(params,backgrnd)
%% 
if nargin == 0
    params.p2f = '/home/rouxf/Exp_EM/test/';
    params.concept_neurons = {'gstq','bb'};
    backgrnd = 'n';
end;
%% search for and get image file-names
[Impath] = [params.p2f];
[Imf] = dir([Impath,filesep,'*.*']);
Imf(1:2) = [];
%if isempty(Imf)
%    error('no image files detected');
%end;
%% get the stimuli
[ID,~] = stimIDmat(Imf);
%% search for a background image
if strcmp(backgrnd,'y')
    [Impath] = [params.basepath,filesep,'image_data',filesep,'background',filesep];
    [Imfb] = dir([Impath,'background*.*']);
    [Imfb2] = dir([Impath,'break*.*']);
%     if isempty(Imfb)
%         error('no background-image file detected');
%     end;
end;
%% make up indices to simulate concept neurons
%The logic here is that each image has a code that links it to an object.
%E.g. gstq links to all images of Quenn Elizabeth. So if a neuron responds
% to one particular object, to find all the images corresponding to that
% object the code of that object has to be entered.

cnx = cell(length(params.concept_neurons),1);
for kt = 1:length(params.concept_neurons)
    [ix,~] = find(strcmp(ID.id,params.concept_neurons(kt)));
    ix = unique(ix);
    
    cnx{kt} = ix;
end; clear ix;


if isempty(cnx{1}) && isempty(cnx{2})
    
    for kt = 1:length(params.concept_neurons)
        [ix] = regexp(ID.id,params.concept_neurons{kt});
        sel = []; k= 0;
        for lt = 1:length(ix)
            if ~isempty(ix{lt})
                k = k+1;
                sel(k) = lt;
            end;
        end;
        
        cnx{kt} = sel;
    end;
end;
%%
%[stim_mat] = organize_stim_mat_concept_neurons(cnx,ID);
% changed 08/07/2016
% %% organize triplet sequences
%[stim_mat] = organize_stim_mat3(cnx,ID); clear cnx;
%if (length(unique(stim_mat.seq(1,:))) == size(stim_mat.seq,2)) ==0
%    error('indices must be unique');
%end;
%%
[stim_mat] = build_unique_cn_pairs(ID,params);

%% randomize sequence order
rand('state',sum(100*clock));
ridx = randperm(size(stim_mat.seq,2));
stim_mat.seq = stim_mat.seq(:,ridx);
ridx = randperm(size(stim_mat.seq,2));
stim_mat.seq = stim_mat.seq(:,ridx);
ridx = randperm(size(stim_mat.seq,2));
stim_mat.seq = stim_mat.seq(:,ridx);

%% generate inter-trial intervals
ntrl = size(stim_mat.seq,2);
[iti1] = generate_iti([1000 2250],ntrl);
[iti2] = generate_iti([250 400],ntrl);
%% save data
params.Impath = Impath;
params.Imf = Imf;
params.Imfb = Imfb;
params.Imfb2 = Imfb2;
params.ID = ID;
params.stim_mat = stim_mat;
params.iti1 = iti1;
params.iti2 = iti2;
return;