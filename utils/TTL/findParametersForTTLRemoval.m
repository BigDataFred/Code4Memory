function [ix] = findParametersForTTLRemoval(ntrl,tw,dataSamples,TTLix,selIxOn,selIxOff,Fs)

ix = [];

muRes1 = zeros(length(ntrl),length(tw));
muRes2 = zeros(length(ntrl),length(tw));
for gt = 1:length(ntrl)
    dum1 = zeros(length(tw),1);
    dum2 = zeros(length(tw),1);
    v1 = ntrl(gt);
    for mt = 1:length(tw)
        [~,res] = cleanARTIFACTfromLFP(dataSamples, TTLix(selIxOn),[32*0.5 32*tw(mt)],v1,2,Fs);
        dum1(mt) = mean(res.^2);
        [~,res] = cleanARTIFACTfromLFP(dataSamples, TTLix(selIxOff),[32*0.5 32*tw(mt)],v1,2,Fs);
        dum2(mt) = mean(res.^2);
    end;
    muRes1(gt,:) = dum1;
    muRes2(gt,:) = dum2;
end;

[ix1,ix2] = find(muRes1 == min(min(muRes1)));
[ix3,ix4] = find(muRes2 == min(min(muRes2)));

ix = [min(ix1) min(ix2) min(ix3) min(ix4)];
