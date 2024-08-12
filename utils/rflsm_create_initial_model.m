function [vel,vel_s,x,z] = rflsm_create_initial_model(param)
% This function creates the velocity model from Nowack et al. (2010)
% 
% input:      param -- struct contains the parameters for radon transform
% output:     vel -- P velocity model
%             vel_s -- S velocity model
%             x -- grid lateral position
%             z -- grid vertical poistion
%             
% August, 2024, Yunfeng Chen, write the function
dx = param.dx;
dz = param.dz;
xmax=param.xmax;
ymax=param.ymax;
xpad=param.xpad;
x=0-xpad:dx:xmax+xpad;
z=0:dz:ymax;
nx=length(x);
nz=length(z);
Vp=zeros(nz,nx);
Vs=zeros(nz,nx);
for j=1:nx
    for i=1:nz
        if x(j)<=290 && z(i)<=74
            Vp(i,j)=6.28;
            Vs(i,j)=3.55;
        elseif x(j)>393 && z(i)<=63
            Vp(i,j)=6.30;
            Vs(i,j)=3.50;
        elseif x(j)>290 && x(j)<=393
            k=(63-74)/(393-290);
            ztmp=k*(x(j)-290)+74;
            if z(i)<=ztmp
                kvp=(6.30-6.28)/(393-290);
                kvs=(3.50-3.55)/(393-290);
                Vp(i,j)=kvp*(x(j)-290)+6.28;
                Vs(i,j)=kvs*(x(j)-290)+3.55;
            else
                Vp(i,j)=8.1;
                Vs(i,j)=4.6;
            end
        else
            Vp(i,j)=8.1;
            Vs(i,j)=4.6;
        end
    end
end

% smooth the velocity model
N=5;
[Vp1,~]=moving_avg(Vp,N,'constant',2);
[Vp2,~]=moving_avg(Vp1,N,'constant');
[Vs1,~]=moving_avg(Vs,N,'constant',2);
[Vs2,~]=moving_avg(Vs1,N,'constant',2);

vel=Vp2;
vel_s=Vs2;

figure;
set(gcf,'Position',[100 100 1000 600],'Color','w')
imagesc(x,z,Vp2);
xlabel('Distance (km)')
ylabel('Depth km')
colorbar
xlim([0 xmax])
set(gca,'FontSize',18)
export_fig('model.png');