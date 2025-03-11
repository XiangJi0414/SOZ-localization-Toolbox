% Area weighted average for blobs.
function mu=mean_amplitude(temp_x,temp_y,temp_z)
utemp_x=unique(temp_x);
utemp_y=unique(temp_y);
[X,Y]=meshgrid(utemp_x,utemp_y);
I=0;
Z=zeros(size(X));
ZA=zeros(size(X));
for i=1:length(temp_x)
    ii=find(temp_x(i)==X(1,:));
    test=Y(:,1);
    jj=find(test==temp_y(i));
    if ~isempty(jj)
        I=I+1;
        Z(jj,ii)=temp_z(i);
        ZA(jj,ii)=1;
    end
end
integ_1=trapz(X(1,:),Z,2);
integ_2=trapz(Y(:,1),integ_1,1);
integ_1=trapz(X(1,:),ZA,2);
integ_3=trapz(Y(:,1),integ_1,1);
mu=integ_2/integ_3;