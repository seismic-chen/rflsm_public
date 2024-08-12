function [dshift] = rflsm_shift_rfs(d,t,vel,vel_s,src,pos,tshift,param)
% This function applies time shift to receiver functions, which is required
% to restore the travel times of the P waves
% seismograms.
% input:      d -- receiver functions
%             t -- time axis
%             vel -- P velocity model
%             vel_s -- S velocity model
%             src -- source time function
%             pos -- positions of point sources
%             tshift -- time lag between point sources
%             param -- struct contains the parameters for radon transform
% output:     dshift -- time-shifted receiver functions
%             
% August, 2024, Yunfeng Chen, write the function
global ishot

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
save_wavefield=1;

[nz,nx] = size(vel);
img=zeros(nz,nx);
tic;
[~,mod_source,~] = kvssmfor_rf(img,vel,vel_s,pos,fpeak,tshift-ph,nt,dt,dx,dz,flow,fhigh,bc,'P',src,save_wavefield);
toc;
%% apply time shift to RF to obtain receiver side wavefield
% cross correlate P wave with RF
Z=squeeze(mod_source(1,:,:));
tdelay=zeros(1,nx);
dshift=zeros(size(d));

for i=1:nx
    if any(d(:,i))
        temp1=Z(:,i);
        temp2=src;
        xc=xcorr(temp1,temp2);
        tax=[-(nt-1):(nt-1)]*dt;
        [~,ind]=max(xc);
        tdelay(i)=tax(ind);
        % remove source time function from RF and filter
        % d(:,i) = bandpassSeis(d(:,i)-src,dt, f1, f2, 4);
        dshift(:,i) = fftShift(d(:,i),t,tdelay(i));

        % plot(d(:,i)); hold on;
        % plot(dshift(:,i))
        % pause;
    end
end
% plot the results
figdir = param.figdir;
fig2=figure;
set(gcf,'Position',[100 100 1200 600],'Color','w');
subplot(121)
imagesc(x,t,Z);
xlabel('Distance (km)');
ylabel('Time (sec)');
title('P wave')
set(gca,'fontsize',14)
colorbar
caxis([-0.2 0.2])
xlim([0,xmax])

subplot(122)
imagesc(x,t,dshift);
xlabel('Distance (km)');
ylabel('Time (sec)');
title('Time shifted RF')
set(gca,'fontsize',14)
colorbar
colormap(seismic(3));
caxis([-0.1 0.1])
xlim([0,xmax])
figname=['rf_tshift.',num2str(ishot),'.png'];
export_fig(fig2,fullfile(figdir,figname));