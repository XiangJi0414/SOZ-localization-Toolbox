% This function is called by cut_results if it is needed to create unimodal
% blobs. 
function [unimodal,temp_x,temp_y,temp_z,cut_CC]=cut_once(temp_x,temp_y,temp_z)
mu_z=mean(temp_z);
[~,ind]=sort(temp_x);
x=temp_x(ind);
z_x=temp_z(ind);
[X,ia] = unique(x);
Z_X=zeros(length(X),1);
for j=1:length(ia)-1
    Z_X(j)=max(z_x(ia(j):ia(j+1)-1));
end
Z_X(j+1)=max(z_x(ia(j+1):end));
[TF_x,P_x] = islocalmin(Z_X/mu_z,'MinProminence',0.1);
P_x=P_x.*TF_x;
[~,ind]=sort(temp_y);
y=temp_y(ind);
z_y=temp_z(ind);
[Y,ia] = unique(y);
Z_Y=zeros(length(Y),1);
for j=1:length(ia)-1
    Z_Y(j)=max(z_y(ia(j):ia(j+1)-1));
end
Z_Y(j+1)=max(z_y(ia(j+1):end));
[TF_y,P_y] = islocalmin(Z_Y/mu_z,'MinProminence',0.1);
P_y=P_y.*TF_y;
[max_P_x,ind_x]=max(P_x);
[max_P_y,ind_y]=max(P_y);
unimodal=false;
if max_P_x>max_P_y
    cuting_time=X(ind_x);
    inds=temp_x<cuting_time;
    cut_CC.x=temp_x(inds);
    cut_CC.y=temp_y(inds);
    cut_CC.z=temp_z(inds);
    temp_x(inds)=[];
    temp_y(inds)=[];
    temp_z(inds)=[];
else
    if max_P_y>0
        cuting_freq=Y(ind_y);
        inds=temp_y<cuting_freq;
        cut_CC.x=temp_x(inds);
        cut_CC.y=temp_y(inds);
        cut_CC.z=temp_z(inds);
        temp_x(inds)=[];
        temp_y(inds)=[];
        temp_z(inds)=[];
    else
        unimodal=true;
        cut_CC.x=temp_x;
        cut_CC.y=temp_y;
        cut_CC.z=temp_z;
    end
end