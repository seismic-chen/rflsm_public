function [dmig] = kvssmadj_rf(din,vel_p,vel_s,pwave,fpeak,tshift,nt,dt,dx,dz,f1,f2,bc,src_type,src,save_wavefield,if_cg)

% this function applies shot profile split step migration
% input:      din -- a common shot gather
%             vel_p -- P velocity model
%             vel_s -- S velocity model
%             dx,dz,dt -- x, z and time interval
%             nt -- number of time samples
%             pwave -- plane wave simulated by a series of point sources
%             fpeak -- peak/central frequency of ricker wavelet
%             tshift -- wavelet shift in time
%             bc =1 apply boundary condition
%             src_type -- 'rikcer' or a user defined source time function
%             src -- source time function, nt samples
%             save_wavefield=1 to save wavefield
%             if_cg=1 to take a vector input 
% output:     dmig -- migrated image for all shots
%
% written by Jinkun Cheng
% Oct. 1, 2020, Yunfeng Chen, reverse the direction of source side
% wavefield to modeling the transition wave from below
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
   din=reshape(din,nt,nx,nshot);
end

dmig = zeros(nz,nx,nshot);

%define f-k axis
nf = 2^nextpow2(nt);
w = 2*pi/(nf*dt)*[0:(nf/2-1),-(nf/2):-1];
ind= (w==0);
w(ind)=1E-30;

%define loop over frequencies
iw1 = floor(f1*dt*nf)+1;
iw2 = floor(f2*dt*nf)+1;

if iw2 > floor(nf/2)+1;
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
        sx=pwave{is};
        ts=tshift;
        stf=src;
    else
        sx=pwave{is};
        ts=tshift{is};
        stf=src{is};
    end
    %sum all point source
    dsrc = zeros(nt,nx);
    for ips=1:length(sx)
        %generate source wavefield
        indx = floor(sx(ips)/dx)+1;
        t=0:dt:(nt-1)*dt;
        if strcmp(src_type,'ricker')
            par=pi*fpeak*(t-ts(ips)-25*dt);
            dsrc(:,indx) = 10*exp(-par.*par).*(1-2*par.*par);  % ricker wavelet
        else
            dsrc(:,indx) = fftShift(stf,t,ts(ips));
        end
    end
    dsc = fft(dsrc,nf,1);
    
    dfx = fft(din(:,:,is),nf,1);
    
    img = zeros(nz,nx);
    
    %loop over frequency
    parfor iw = iw1:iw2
        % source side forward progation
        [swave,~] = sspropog1_rf(dsc(iw,:),vavg_p,du_p,nx,dx,nz,dz,w(iw), 1,bc,'source',  save_wavefield);
        [rwave,~] = sspropog1_rf(dfx(iw,:),vavg_s,du_s,nx,dx,nz,dz,w(iw),-1,bc,'receiver',save_wavefield);
        
        %apply imaging condition
        img = img + real(rwave.*conj(swave));
        
    end
    
    dmig(:,:,is) = img./nf;
    
end

if if_cg==1
   dmig=dmig(:); 
end
return
