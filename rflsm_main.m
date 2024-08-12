% Yunfeng Chen, Global Seismology Group, Zhejiang University
%
% Receiver function least-squares migration package 
% 
% If you use this code, please consider citing:
% Chen et al. (2024), IEEE TGRS.
% "Least-squares migration imaging of receiver functions"
% 
% Please contact yunfeng_chen@zju.edu.cn if you need more assistance with this code. 
% 
% Copyright 2024 Yunfeng Chen
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
clear; clc; close all;
addpath utils/
addpath ssa/
addpath radon/
addpath OneWay_RF/
addpath export_fig/
% global variable
global ishot
%% load paramerters
param = rflsm_init_param();
datadir = param.datadir;
figdir = param.figdir;

is_radon=param.is_radon;
is_ssa=param.is_ssa;

% location of profile
lon1 = param.lon1;
lat1 = param.lat1;
lon2 = param.lon2;
lat2 = param.lat2;

% grid info
dx = param.dx;
nx = param.nx;
xpad=param.xpad;
xmax=param.xmax;

minsnr = param.minsnr;
mintrace = param.mintrace;
%% load data
load(fullfile(datadir,'decon_xf_0.02_3_gauss2.5_src.mat'));
%% prepare velocity model for migration
[vel,vel_s,x,z] = rflsm_create_initial_model(param);
%% loop over all shots (events)
eventid=unique({log.id});
nshot=length(eventid);

for ishot=1:length(eventid)
    disp(['Processing shot ', num2str(ishot)]);
    keep = strcmp({log.id},eventid{ishot});
    sub=log(keep);

    % remove bad traces
    remove= [sub.snr] < minsnr;
    sub(remove)=[];
    if length(sub) < mintrace
        continue
    end

    % project the stations onto a profile 
    slat=[sub.slat];
    slon=[sub.slon];
    slatp=zeros(size(slat));
    slonp=zeros(size(slon));
    dist=zeros(size(slat));
    for i=1:length(slat)
        [slatp(i),slonp(i),dist(i)] = proj_point_to_gcp(lat1,lon1,lat2,lon2,slat(i),slon(i));
    end

    % calculate the receiver location along the profile
    [deg0,~]= distance(lat1,lon1,slatp,slonp);
    rx = deg0*2*pi*6371/360;
    param.rx = rx;
    %% Preprocessing with Radon transform (Zhang et al., 2022, GJI)
    if is_radon
        ittax=sub(1).ittax;
        dt=ittax(2)-ittax(1);
        nt=length(ittax);
        t=(0:nt-1)*dt;
        h=[sub.dist];
        d_z=[sub.Z];
        d_r=[sub.R];
        [dp_z,dp_r,itr] = rflsm_radon(d_z,d_r,t,h,param);
    else
        itr = [sub.itr];
    end
    %% Binning
    d = zeros(nt,nx);
    for i=1:nx
        idx=rx>=x(i)-dx/2 & rx<=x(i)+dx/2;
        if sum(idx)>0
            d(:,i) = mean(itr(:,idx),2);
        end
    end
    %% Parameters for SSA/CAZDOW reconstruction
    if is_ssa
        i1=xpad/dx+1;
        i2=nx-xpad/dx-1;
        din = d(:,i1:i2);
        [dout] = rflsm_ssa(din,ittax,param);
        d(:,i1:i2) = dout;
    end
    %% forward propagation of P wavefield
    vp = mean(vel(end,:));
    rayp = mean([sub.rayp]);
    baz = mean([sub.baz]);
    [src,pos,tshift] = rflsm_create_src(dt,nt,rayp,baz,vp,param);

    % taper RF to remove later conversions
    wl=floor(((40-ittax(1))/dt+1));
    w = [tukeywin(wl,0.75); zeros(size(d,1)-wl,1)];
    win=w*ones(1,size(d,2));
    d=d.*win;

    rayps = [sub.rayp];
    bazs = [sub.baz];
    dshift = rflsm_shift_rfs(d,ittax,vel,vel_s,src,pos,tshift,param);
    %% Migration
    [mig,~] = rflsm_migration(dshift,ittax,vel,vel_s,src,pos,tshift,param);
    dmig(:,:,ishot) = mig;
    %% LSM
    [migls,~] = rflsm_lsm(dshift,ittax,vel,vel_s,src,pos,tshift,param);
    dmigls(:,:,ishot) = migls;
    %% close all figures
     close all
end
%% save migration results
save(fullfile(figdir,'mig.mat'), 'dmig', 'x', 'z')
save(fullfile(figdir,'migcg.mat'), 'dmigls', 'x', 'z')
d2d=mean(dmig,3);
d2dls=mean(dmigls,3);
%% plot migration results
figure();
set(gcf,'Position',[100 100 800 800],'color','w')
subplot(211)
imagesc(x,z, d2d/max(d2d(:))); hold on;
axis([0 xmax 0 200])
xlabel('Distance (km)');
ylabel('Depth (km)');
title('Migration image')
set(gca,'fontsize',14)
cmax=0.3;
caxis([-cmax cmax]);
colorbar
text(-0.12,0.98,'a)','Units','normalized','FontSize',18)

subplot(212)
imagesc(x,z, d2dls/max(d2dls(:))); hold on;
axis([0 xmax 0 200])
xlabel('Distance (km)');
ylabel('Depth (km)');
title('LSM')
set(gca,'fontsize',14)
caxis([-cmax cmax]);
colorbar
colormap(seismic(3));
text(-0.12,0.98,'b)','Units','normalized','FontSize',18)
figname='mig_compare_real.png';
export_fig(fullfile(figdir,figname));