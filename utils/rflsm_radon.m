function [dp_z,dp_r,itr_denoise] = rflsm_radon(d_z,d_r,t,h,param)

% This function applies Radon transform filters to vertical and radial 
% component seismograms.
% input:      d_z -- vertical component seismogram
%             d_r -- radiao component seismogram
%             t -- time axis
%             h -- offset axis
%             param -- struct contains the parameters for radon transform
% output:     dp_z -- filtered vertical component seismogram 
%             dp_r -- filtered radial component seismogram 
%
% August, 2024, Yunfeng Chen, write the function

global ishot

disp('Radon transform begins')

% load parameters
pmax = param.pmax;
gauss = param.gauss;
itmax = param.itmax;
minderr = param.minderr;
ph = param.ph; % phase delay
N1 = param.N1;
N2 = param.N2;
rx = param.rx;
dp = param.dp;
VB = 1; % Allow verbose output

p = [-pmax:dp:pmax];
dt = t(2)-t(1);
nt=length(t);
nx=length(h);
np =length(p);

Param.h=h;
Param.v=1./p;
Param.nt=nt;
Param.dt=dt;
Param.type=1;
ma=zeros(nt,np);
% Apply LS radon to Z component
[mi,misfit] = yc_pcg(@radon_op,Param,d_z,zeros(size(ma)),N1,N2,VB);

% Prediction (Use inverted velocity gather to estimate a preduction of the
% data)
dp_z = radon_op(mi,Param,1);

% Apply LS radon to R component
[mi,misfit] = yc_pcg(@radon_op,Param,d_r,zeros(size(ma)),N1,N2,VB);
dp_r = radon_op(mi,Param,1);

% deconvolution
itr = zeros(nt,nx);
for n=1:size(dp_r,2)
    R=d_r(:,n);
    Z=d_z(:,n);
    [itr(:,n),itrms] = makeRFitdecon_la_norm(R,Z,dt,nt,ph,gauss,itmax,minderr);
end

itr_denoise = zeros(nt,nx);
for n=1:size(dp_r,2)
    R=dp_r(:,n);
    Z=dp_z(:,n);
    [itr_denoise(:,n),itrms] = makeRFitdecon_la_norm(R,Z,dt,nt,ph,gauss,itmax,minderr);
end
disp('Radon transform done')
%% plot the raw and denoised data
figdir = param.figdir;
plot_type='wiggle';
ittax = t-ph;
fig=figure;
set(gcf,'Position',[0 0 1200 700],'Color','w')
subplot(231)
switch plot_type
    case 'wiggle'
        wigb(d_z,2,rx,t)
    case 'imagesc'
        imagesc(rx,t,d_z)
        cmax=3*rms(d_z(:));
        caxis([-cmax cmax])
        colorbar
end
title('Vertical');xlabel('Distance (km)');ylabel('Time (sec)')
set(gca,'fontsize',14)
text(-0.2,0.98,'(a)','Units','normalized','FontSize',18)
subplot(232)
switch plot_type
    case 'wiggle'
        wigb(d_r,2,rx,t)
    case 'imagesc'
        imagesc(rx,t,d_r)
        caxis([-cmax cmax])
        colorbar
end
title('Horizontal');xlabel('Distance (km)');ylabel('Time (sec)')
set(gca,'fontsize',14)
text(-0.2,0.98,'(b)','Units','normalized','FontSize',18)
subplot(233)
switch plot_type
    case 'wiggle'
        wigb(itr,4,rx,ittax)
    case 'imagesc'
        imagesc(rx,ittax,itr)
        caxis([-0.3 0.3])
        colorbar
        colormap(seismic(3));
end
title('Receiver function');xlabel('Distance (km)');ylabel('Time (sec)')
set(gca,'fontsize',14)
text(-0.2,0.98,'(c)','Units','normalized','FontSize',18)
ylim([-5 30])

subplot(234)
switch plot_type
    case 'wiggle'
        wigb(dp_z,2,rx,t)
    case 'imagesc'
        imagesc(rx,t,dp_z)
        caxis([-cmax cmax])
        colorbar
end
title('Vertical');xlabel('Distance (km)');ylabel('Time (sec)')
set(gca,'fontsize',14)
text(-0.2,0.98,'(d)','Units','normalized','FontSize',18)
subplot(235)
switch plot_type
    case 'wiggle'
        wigb(dp_r,2,rx,t)
    case 'imagesc'
        imagesc(rx,t,dp_r)
        caxis([-cmax cmax])
        colorbar
end
title('Horizontal');xlabel('Distance (km)');ylabel('Time (sec)')
set(gca,'fontsize',14)
text(-0.2,0.98,'(e)','Units','normalized','FontSize',18)

subplot(236)
switch plot_type
    case 'wiggle'
        wigb(itr_denoise,4,rx,ittax)
    case 'imagesc'
        imagesc(rx,ittax,itr_denoise)
        caxis([-0.3 0.3])
        colorbar
        colormap(seismic(3));
end
title('Receiver function');xlabel('Distance (km)');ylabel('Time (sec)')
set(gca,'fontsize',14)
text(-0.2,0.98,'(f)','Units','normalized','FontSize',18)
ylim([-5 30])
figname=['z_r_rf.',num2str(ishot),'_',plot_type,'.png'];
export_fig(fig,fullfile(figdir,figname));