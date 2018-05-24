%%
xb = [0 1 .5];
xb = repmat(xb,[1 100]);
xb = xb + .1*randn(1,length(xb));

y = [1 2 3];
y = repmat(y,[1 100]);

orig = xb';
y = y';

ix1 = 1:length(orig);
ix2 = zeros(length(orig),length(orig)-1);
for it = 1:length(orig)
    ix2(it,:) = setdiff(1:length(orig),ix1(it));
end;

acc = zeros(1,length(ix1));
for it = 1:length(ix1)   
    [c] = classify(orig(ix1(it))',orig(ix2(it,:)),y(ix2(it,:)),'linear');
    acc(it) = (c==y(ix1(it)));
end;

sum(acc)/length(ix1)