%%
addpath('/media/rouxf/rds-share/Fred/code/mcode/custom/helper/logfile_readers/');

%%
p2f = '/media/rouxf/rds-share/Archive/MICRO/P08/log/seqPRED/';
fn = dir([p2f,'*.txt'])

%%
D1 = {};%cell( 1, length( fn ) );
D2 = {};%cell( 1, length( fn ) );
D3 = {};%cell( 1, length( fn ) );

b1 = {};%cell( 1, length( fn ) );
b2 = {};%cell( 1, length( fn ) );
b3 = {};%cell( 1, length( fn ) );

%%
cnt = 0;
for lt = 1:length( fn )
    
    fid = fopen([p2f,fn(lt).name],'r')
    
    %%
    [dat] = textscan(fid,'%s\t');
    dat = dat{1};
    
    hdr = dat(1:11);
    hdr(1:3) = [];
    
    begIx = find(strcmp(dat,'Learning-Beg'));
    endIx = find(strcmp(dat,'Learning-End'));
    
    dat = dat(begIx+1:endIx-1);
    
    %%
    LogDat = [];
    idx = 1:length(hdr);
    ntrl = length(dat)/length(hdr);
    for it = 1:ntrl
        LogDat = [LogDat;dat(idx)'];
        idx = idx+length(hdr);
    end;
    
    %%
    if ~isempty( LogDat )
        
        cnt = cnt+1;
        
        ix1 = [];
        ix2 = [];
        c1 = 0;
        c2 = 0;
        
        for it = 1:size(LogDat,1)
            
            if strcmp(LogDat(it,2),'cn_1.jpg') && strcmp(LogDat(it+1,2),'cn_2.jpg')
                c1 =c1+1;
                ix1(c1,1) = it;
            end;
            
            if strcmp(LogDat(it,2),'cn_2.jpg') && strcmp(LogDat(it+1,2),'cn_1.jpg')
                c2 =c2+1;
                ix2(c2,1) = it;
            end;
            
        end;
        
        cn1Ix = ix2+1;
        cn2Ix = ix1+1;
        
        %%
        c1RT = str2double(  LogDat(cn1Ix-1,end) ).*1e3;
        c2RT = str2double(  LogDat(cn2Ix-1,end) ).*1e3;
        
        c1RT2 = str2double(  LogDat(cn1Ix,end) ).*1e3;
        c2RT2 = str2double(  LogDat(cn2Ix,end) ).*1e3;
        
        D1{cnt} = c1RT2-c1RT;
        D2{cnt} = c2RT2-c2RT;
        
        r_ix = setdiff(1:size(LogDat,1),[cn1Ix;cn2Ix]);
        r_ix(r_ix==1) = [];
        
        nRT =  str2double( LogDat(r_ix-1,end) );
        nRT2 = str2double( LogDat(r_ix,end) );
        
        D3{cnt} =  nRT2-nRT;
        
        %%
        X = 1:length(D1{cnt});
        %b1{lt} = robustfit(X',D1{lt});
        b1{cnt} = regress(D1{cnt},[X' ones(length(X),1)]);
        
        X = 1:length(D2{cnt});
        %b2{lt} = robustfit(X',D2{lt});
        b2{cnt} = regress(D2{cnt},[X' ones(length(X),1)]);
        
        X = 1:length(D3{cnt});
        %b3{lt} = robustfit(X',D3{lt});
        b3{cnt} = regress(D3{cnt},[X' ones(length(X),1)]);
    end;
end;

%%
for lt = 1:length(D1)
    
    figure;
    a = [];
    for it = 1:4
        subplot(2,3,it);
        a = [a gca];
        hold(a(it),'on');
    end;
    
    
    %yp = b1{lt}(2)+b1{lt}(1)*[1:length(D1{lt})]';
    yp = [[1:length(D1{lt})]' ones(length(D1{lt}),1)]*b1{lt};
    plot(a(1),D1{lt},'bs-');
    plot(a(1),1:length(yp),yp,'r-');
    title(a(1),'Sequence #1');
    
    %yp = b2{lt}(2)+b2{lt}(1)*[1:length(D2{lt})]';
    yp = [[1:length(D2{lt})]' ones(length(D2{lt}),1)]*b2{lt};
    plot(a(2),D2{lt},'bs-');
    plot(a(2),1:length(yp),yp,'r-');
    title(a(2),'Sequence #2');
    
    %yp = b3{lt}(2)+b3{lt}(1)*[1:length(D3{lt})]';
    yp = [[1:length(D3{lt})]' ones(length(D3{lt}),1)]*b3{lt};
    plot(a(3),D3{lt},'bs-');
    plot(a(3),1:length(yp),yp,'r-');
    title(a(3),'Random order');
    
    bar(a(4),1:3,[b1{lt}(1) b2{lt}(1) b3{lt}(1)]);
    
    axis(a,'tight');
    set(a(4),'XLim',[0 4]);
    set(a(4),'XTick',1:3);
    set(a(4),'XTickLabel',{'S1' 'S2' 'R'});
    ylabel(a(4),'Slope [ms]');
    
    for it = 1:3
        xlabel(a(it),'Repetition');
        ylabel(a(it),'Reaction time [ms]');
    end;
end;

%%
figure;plot(str2double(LogDat(sort([cn1Ix;cn2Ix]),end)),'bs-')
figure;plot(str2double(LogDat(sort([cn1Ix]),end)),'bs-')
figure;plot(str2double(LogDat(sort([cn2Ix]),end)),'bs-')

%%
figure;
subplot(121);
hold on;
h = bar([1 2],[nanmean(fRT) nanmean(pRT)]);
plot([1 1],[nanmean(fRT) nanmean(fRT)+nanstd(fRT)/sqrt(length(fRT)-1)],'k','LineWidth',3);
plot([1 1],[nanmean(fRT)-nanstd(fRT)/sqrt(length(fRT)-1) nanmean(fRT)+nanstd(fRT)/sqrt(length(fRT)-1)],'k','LineWidth',3);
plot([.9 1.1],[nanmean(fRT)+nanstd(fRT)/sqrt(length(fRT)-1) nanmean(fRT)+nanstd(fRT)/sqrt(length(fRT)-1)],'k','LineWidth',3);
plot([.9 1.1],[nanmean(fRT)-nanstd(fRT)/sqrt(length(fRT)-1) nanmean(fRT)-nanstd(fRT)/sqrt(length(fRT)-1)],'k','LineWidth',3);

plot([2 2],[nanmean(pRT)-nanstd(pRT)/sqrt(length(pRT)-1) nanmean(pRT)+nanstd(pRT)/sqrt(length(pRT)-1)],'k','LineWidth',3);
plot([1.9 2.1],[nanmean(pRT)+nanstd(pRT)/sqrt(length(pRT)-1) nanmean(pRT)+nanstd(pRT)/sqrt(length(pRT)-1)],'k','LineWidth',3);
plot([1.9 2.1],[nanmean(pRT)-nanstd(pRT)/sqrt(length(pRT)-1) nanmean(pRT)-nanstd(pRT)/sqrt(length(pRT)-1)],'k','LineWidth',3);

xlim([0 3]);ylim([.45 .65]);
set(h,'FaceColor',[.75 .75 .75]);

subplot(122);
hold on;
h = bar([1 2],[nanmean(c1RT) nanmean(c2RT)]);

plot([1 1],[nanmean(c1RT) nanmean(c1RT)+nanstd(c1RT)/sqrt(length(c1RT)-1)],'k','LineWidth',3);
plot([1 1],[nanmean(c1RT)-nanstd(c1RT)/sqrt(length(c1RT)-1) nanmean(c1RT)+nanstd(c1RT)/sqrt(length(c1RT)-1)],'k','LineWidth',3);
plot([.9 1.1],[nanmean(c1RT)+nanstd(c1RT)/sqrt(length(c1RT)-1) nanmean(c1RT)+nanstd(c1RT)/sqrt(length(c1RT)-1)],'k','LineWidth',3);
plot([.9 1.1],[nanmean(c1RT)-nanstd(c1RT)/sqrt(length(c1RT)-1) nanmean(c1RT)-nanstd(c1RT)/sqrt(length(c1RT)-1)],'k','LineWidth',3);

plot([2 2],[nanmean(c2RT)-nanstd(c2RT)/sqrt(length(c2RT)-1) nanmean(c2RT)+nanstd(c2RT)/sqrt(length(c2RT)-1)],'k','LineWidth',3);
plot([1.9 2.1],[nanmean(c2RT)+nanstd(c2RT)/sqrt(length(c2RT)-1) nanmean(c2RT)+nanstd(c2RT)/sqrt(length(c2RT)-1)],'k','LineWidth',3);
plot([1.9 2.1],[nanmean(c2RT)-nanstd(c2RT)/sqrt(length(c2RT)-1) nanmean(c2RT)-nanstd(c2RT)/sqrt(length(c2RT)-1)],'k','LineWidth',3);

xlim([0 3]);ylim([.45 .65]);
set(h,'FaceColor',[.75 .75 .75]);

%%
sRT1 = [];
sRT1(:,1) = str2double(  LogDat(fIx,end) );
sRT1(:,2) = str2double(  LogDat(fIx+1,end) );

sRT2 = [];
sRT2(:,1) = str2double(  LogDat(pIx,end) );
sRT2(:,2) = str2double(  LogDat(pIx+1,end) );
