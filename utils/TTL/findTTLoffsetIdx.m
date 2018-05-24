function [selIx0] = findTTLoffsetIdx(ttls,selIx,offSet)


dum1 = selIx(ttls(selIx+1)==0)+offSet(1);
dum2 = find(ttls(selIx+1) ~=0);
dum2 = selIx(dum2(ttls(selIx(dum2)+offSet(2)) ==0))+offSet(2);

[selIx0] = sort([dum1 dum2]);