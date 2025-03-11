function [spikes,artifacts]=find_spikes_basic(z,f,time,s_range,prominence,art_range)


[~,S1]=min(abs(f-s_range(1)));
[~,S2]=min(abs(f-s_range(2)));

test=mean(z(S2:S1,:),1);
test=test/mean(test);
[~,locs]=findpeaks(test,'MinPeakProminence',prominence);
spikes=time(locs);


[~,S1]=min(abs(f-art_range(1)));
[~,S2]=min(abs(f-art_range(2)));
test=mean(z(S2:S1,:),1);
temp=test/mean(test);
[~,locs]=findpeaks(temp,'MinPeakProminence',prominence);
I=0;
artifacts=zeros(1e6,1);
for i=1:length(locs)
    if test(locs(i))>mean(z(:,locs(i)))
        I=I+1;
        artifacts(I)=time(locs(i));
    end
end
artifacts(I+1:end)=[];