load MioTimesUpdatedExplicit.txt;
tq=MioTimesUpdatedExplicit(1,3);
sp = tq./MioTimesUpdatedExplicit(:,3);
plot(MioTimesUpdatedExplicit(:,2), sp, 'r*-')