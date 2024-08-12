function [mig,dp] = rflsm_migration(d,t,vel,vel_s,src,pos,tshift,param)
% This function applies migration to receiver functions
% 
% input:      d -- time-shifted receiver functions
%             t -- time axis
%             vel -- P velocity model
%             vel_s -- S velocity model
%             src -- source time function
%             pos -- positions of point sources
%             tshift -- time lag between point sources
%             param -- struct contains the parameters for radon transform
% output:     mig -- migration image
%             dp -- predicted receiver functions
%             
% August, 2024, Yunfeng Chen, write the function

global ishot

disp('Migration begins')
% parameters for forward propagation
fpeak = 10.; %not used
dt = t(2)-t(1);
nt = length(t);
dx = param.dx;
dz = param.dz;
x = param.x;
z = param.z;
nx = param.nx;
nz = param.nz;
xmax=param.xmax;
ph=abs(t(1));
bc = 1;
flow = param.flow;
fhigh = param.fhigh;
save_wavefield=0;

tic;
[mig] = kvssmadj_rf(d,vel,vel_s,pos,[], ...
    tshift-ph,nt,dt,dx,dz,flow,fhigh,bc,'P',src,save_wavefield);
toc;

[dp,~,~] = kvssmfor_rf(mig,vel,vel_s,pos,[], ...
    tshift-ph,nt,dt,dx,dz,flow,fhigh,bc,'P',src,save_wavefield);
disp('Migration done')

% plot results
figdir = param.figdir;
fig=figure;
set(gcf,'Position',[100 100 1800 600],'color','w')
subplot(131)
imagesc(x,t, reshape(d/max(d(:)),nt,nx));
xlim([0 xmax])
xlabel('Distance (km)');
ylabel('Time (sec)');
title('Time-shifted RF')
set(gca,'fontsize',14)
colorbar
caxis([-1 1]);
text(-0.2,0.98,'(a)','Units','normalized','FontSize',18)

subplot(132)
imagesc(x,z, mean(mig,3)); hold on;
xlim([0 xmax])
colormap(seismic(3));
xlabel('Distance (km)');
ylabel('Depth (km)');
title('Migration image')
set(gca,'fontsize',14)
colorbar
cmax=rms(abs(mig(:)));
caxis([-3*cmax 3*cmax]);
text(-0.2,0.98,'(b)','Units','normalized','FontSize',18)

subplot(133)
imagesc(x,t, reshape(dp/max(dp(:)),nt,nx));
xlim([0 xmax])
colorbar
caxis([-1 1]);
xlabel('Distance (km)');
ylabel('Depth (km)');
title('Predicted data')
set(gca,'fontsize',14)
text(-0.2,0.98,'(c)','Units','normalized','FontSize',18)
figname=['mig.',num2str(ishot),'.png'];
export_fig(fig,fullfile(figdir,figname));