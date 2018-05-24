%%
c = 0;
C = [];
for it = 1:5
    for jt = 1:5
        for kt = 1:5
            c = c+1;
            C(c,:) = [it/5 jt/5 kt/5];
        end;
    end;
end;
C(find(diff(diff(C,[],2),[],2)==0),:) = [];
%%
nchans = 24;%length(spike_data{2}.label);
figure;

for jt = 1:nchans
    
    cID = unique(spike_data{2}.unit{jt});
    
    X = squeeze(spike_data{2}.waveform{jt});
    X = X';
    [COEFF, SCORE] = pca(X);
    
    subplot(nchans/8,8,jt);
    hold on;
    plot(SCORE(:,1),SCORE(:,2),'.','Color',[.5 .5 .5]);
    
    x = 1:size(C,1);
    
    ptc = [];
    for it = 1:length(cID)
        ix = find(spike_data{2}.unit{jt} == cID(it));
        
        isi = [];
        trl = spike_data{2}.trial{jt}(ix);
        for kt = 1:length(trl)
            ts = spike_data{2}.time{jt}(find(spike_data{2}.trial{jt} == trl(kt) & spike_data{2}.unit{jt} == cID(it)));
            isi = [isi diff(ts)];
        end;
        pct(it) = length(find(isi < 0.003))/length(isi);
        
        if pct(it) < 0.03
            if ~isempty(x)
                sel = randperm(length(x));
                a = x(sel(1));
                x(sel(1)) = [];
                sel = a;
                
                plot(SCORE(ix,1),SCORE(ix,2),'.','Color',C(sel,:));
                
            end;
        end;
    end;
    
    axis tight;
    title(spike_data{2}.label(jt));
    
    
end;