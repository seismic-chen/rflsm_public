function [mod,mod_source,mod_receiver] = kvssmfor_rf(img,vel_p,vel_s,pwave,fpeak,tshift,nt,dt,dx,dz,f1,f2,bc,src_type,src,save_wavefield,if_cg)

% this function applies shot profile split step migration
% input:      rmig -- the image (same size as vel)
%             vel_p -- P velocity model
%             vel_s -- S velocity model
%             dx,dz,dt -- x, z and time interval
%             nt -- number of time samples
%             pwave -- cell array one for each event plane wave simulated by a series of point sources
%             fpeak -- peak/central frequency of ricker wavelet
%             tshift -- cell array one for each event wavelet shift in time
%             bc =1 apply boundary condition
%             src_type -- 'rikcer' or a user defined source time function
%             src -- cell array one for each event source time function, nt samples
%             save_wavefield=1 to save wavefield
%             if_cg=1 to take a vector input 
% output:     mod -- modeled data
%             mod_source -- modeled source wavefield
%             mod_receiver -- modeled receiver wavefield
%
% written by Jinkun Cheng
% modified by Yunfeng Chen
% Oct. 1, 2020, reverse the direction of source side
% wavefield to model the plane wave incidence from below
% add option to use a user defined source-time function
% add option to save the wavefield
% add option to take input data as a vector
% Apr. 10, 2022, Yunfeng Chen, modify the names of variables to be more
% readable

if ~exist('if_cg','var')
  if_cg=0;
end
if iscell(pwave)
    nshot = length(pwave);
else
    nshot = 1;
end
[nz,nx] = size(vel_p);

if if_cg==1
   img=reshape(img,nz,nx,nshot);
end

mod1 = zeros(nt,nx,nshot);

%zero padding in f and kx directions (important)
nf = 2^nextpow2(nt);
ns = 2^nextpow2(nx);

%define f-k axis with zero-padding in time axis
w = 2*pi/(nf*dt)*[0:(nf/2-1),-(nf/2):-1];
ind= (w==0);
w(ind)=1E-30;


%define loop over frequencies
iw1 = floor(f1*dt*nf)+1;
iw2 = floor(f2*dt*nf)+1;

% if frequency exceeding the nyquest freq
if iw2 > floor(nf/2)+1
    iw2=floor(nf/2)+1;
end

%get vavg and slowness purtabation du
du_p = zeros(nz,nx);
vavg_p = zeros(nz,1);
du_s = zeros(nz,nx);
vavg_s = zeros(nz,1);

for iz = 1:nz
    vavg_p(iz) = mean(vel_p(iz,:));
    vavg_s(iz) = mean(vel_s(iz,:));
    for ix = 1:nx
        du_p(iz,ix) = (1/vel_p(iz,ix)-1/vavg_p(iz));
        du_s(iz,ix) = (1/vel_s(iz,ix)-1/vavg_s(iz));
    end
end

%loop over shots
for is = 1:nshot
    if nshot==1
        sx=pwave{1};
        ts=tshift;
        stf=src;
    else
        sx=pwave{is};
        ts=tshift{is};
        stf=src{is};
    end
    %sum all point sources
    dsrc = zeros(nt,nx);
    for ips=1:length(sx)
        %generate source wavefield
        indx = floor(sx(ips)/dx)+1;
        t=0:dt:(nt-1)*dt;
        if strcmp(src_type,'ricker')
            par=pi*fpeak*(t-ts(ips)-25*dt);
            dsrc(:,indx) = 10*exp(-par.*par).*(1-2*par.*par);  % ricker wavelet
        else
            dsrc(:,indx) = fftShift(stf,t,ts(ips)); % user defined wavelet
        end
    end
    dsc = fft(dsrc,nf,1);
    
    outf = zeros(nf,ns);
    outf_source=zeros(nz,nf,ns);
    outf_receiver=zeros(nz,nf,ns);
    
    % do frequency loop in parallel
    simg = img(:,:,is);
    parfor iw = iw1:iw2        
        % for source side wavefield
        [swave,swave_full] = sspropog1_rf(dsc(iw,:),vavg_p,du_p,nx,dx,nz,dz,w(iw),1,bc,'source',save_wavefield);
        
        % for receiver side wavefield
        % upward P-wave propagation, which is only used for RF simulation
%         pimg = img(:,:,is)*1.;  % fake P reflectivity
%         [rwave_p,rwave_full_p] = sspropog1_rf((pimg.*swave),vavg,du,nx,dx,nz,dz,w(iw),1,bc,'receiver',save_wavefield);
        % upward S-wave propagation using shear velocities
        [rwave_s,rwave_full_s] = sspropog1_rf((simg.*swave),vavg_s,du_s,nx,dx,nz,dz,w(iw),1,bc,'receiver',save_wavefield);
        
        outf(iw,:) = rwave_s;
        
        if save_wavefield
            outf_source(:,iw,:)=swave_full;
            outf_receiver(:,iw,:)=rwave_full_s;
            % outf_receiver(:,iw,:)=rwave_full_s+rwave_full_p; % only used
            % when simulating P phase in RF
        end
    end

    %transform back to time domain
    r = real(ifft2(outf));
    mod1(:,:,is) = r(1:nt,1:nx);
    if save_wavefield
        mod_source = zeros(nz,nt,nx);
        mod_receiver = zeros(nz,nt,nx);
        for iz=1:nz
            s_full=real(ifft2(squeeze(outf_source(iz,:,:))));
            r_full=real(ifft2(squeeze(outf_receiver(iz,:,:))));
            mod_source(iz,:,:) = s_full(1:nt,1:nx);
            mod_receiver(iz,:,:) = r_full(1:nt,1:nx);
        end
    else
        mod_source=[];
        mod_receiver=[];
    end
   % mod = sum(mod1,3);
   mod = mod1;
end

if if_cg==1
    mod=mod1(:);
end
return
