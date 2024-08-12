function [wavf,wavf_full] = sspropog1_rf(in,vavg,du,nx,dx,nz,dz,w,iflag,bc,wavf_type,save_wavefield)

% phase shift wave propogator with split-step correction for 1 frequency
% note that it is only used to extrapolate from 1 to nz (not nz to 1)

% if iflag == -1 do migration propogation
% if iflag == 1  do modeling propogation

% if dz > 0 apply downward propogation
% if dz < 0 apply upward propogation

% input:  in   -- data/model to be propogated
%          w   -- angular frequency
%       vavg   -- average velocity of each depth
%         du   -- purtubation of slowness (size of vel model)
%         dx   -- x spatial interval
%         dz   -- z spatial interval
%         dt   -- time interval
%  wavf_type   -- wavefiled type (source or receiver)
%  save_wavefield -- set to 1 to save wavefields at each depth
% output:
%         wavf -- propogated wavefield at a given frequency
%                with split step corrections
%      wav_full -- 
% written by Jinkun Cheng
% modified by Yunfeng Chen for receiver function migration
% add option for upward progation of source or receiver wavefields
% source wavefield is used to simulate plane P wave from bottom to up
% receiver wavefield is used to simulate converted S wave from bottom to up

we = ones(1,nx);
if bc==1;
    Lx=20;
    ix = 1:1:Lx;
    tap1 = exp(-(0.015*(Lx-ix)).^2);
    we = [tap1,ones(1,nx-2*Lx),fliplr(tap1)];
end

%define kx axis
ns = 2^nextpow2(nx);
Kx = 2*pi/(dx*ns)*[0:(ns/2-1),-(ns/2):-1];

if save_wavefield
    wavf_full=zeros(nz,ns);
else
    wavf_full=[];
end

if iflag == -1                 %do migration
    
    M0_old = fft(in,ns,2);
    wavf = zeros(nz,nx);
    
    for iz = 1:nz
        
        S0 = zeros(1,ns);
        for ikx = 1:ns
            arg = w^2/vavg(iz)^2 - Kx(ikx)^2;
            
            if arg >= 0
                S0(ikx) = exp(1i*sqrt(arg)*dz);
            else
                S0(ikx) = exp(-1*(sqrt(-arg))*abs(dz));
            end
            
        end
        M0 = M0_old.*S0;
        
        %apply correction in f-x domain
        P0 = ifft(M0,ns,2);
        P1 = P0(1:nx).*exp(1i*w*dz*du(iz,:));     
        P1 = P1.*we;
        wavf(iz,:) = P1;
        M0_old = fft(P1,ns,2);
        if save_wavefield
            wavf_full(iz,:)=M0_old;
        end
    end
end

if iflag == 1   % do modeling
    
    mfk = fft(in,ns,2);
    switch wavf_type
        case 'source'
            wavf = zeros(nz,nx);
            mfk_old=mfk;
        case 'receiver'
            wavf = zeros(1,ns);
    end
    
    for iz = nz:-1:1      %from bottom to up
        
        S0 = zeros(1,ns);
        for ikx = 1:ns
            arg = w^2/vavg(iz)^2 - Kx(ikx)^2;
            
            if arg >= 0
                S0(ikx) = exp(-1i*sqrt(arg)*dz);
            else
                S0(ikx) = exp(-1*(sqrt(-arg))*abs(dz));
            end
            
        end
        switch wavf_type
            case 'source' % forward propagating source wavefield output in fx domain
                % deal with the special case at the bottom of the model
                % i.e., iz=nz, where the wavefield is just the source time
                % function (wavelet), do not propagate the wavefield
                if iz == nz
                    P0 = ifft(mfk,ns,2);
                    P1 = P0(1:nx);
                    P1 = P1.*we;  
                    mfk_old = fft(P1,ns,2);
                end
                wavf(iz,:)=P1;
                mfk=mfk_old.*S0;
                P0 = ifft(mfk,ns,2);
                P1 = P0(1:nx).*exp(-1i*w*dz*du(iz,:));              % adj ss correction
                P1 = P1.*we;
                mfk_old = fft(P1,ns,2);
                if save_wavefield
                    wavf_full(iz,:)=mfk_old;
                end
                
%                 mfk=mfk_old.*S0;
%                 P0 = ifft(mfk,ns,2);
%                 P1 = P0(1:nx).*exp(-1i*w*dz*du(iz,:));              % adj ss correction
%                 P1 = P1.*we;
%                 wavf(iz,:)=P1;
%                 mfk_old = fft(P1,ns,2);
%                 if save_wavefield
%                     wavf_full(iz,:)=mfk_old;
%                 end
            case 'receiver' % forward propagating receiver wavefield output in fk domain
                P0 = ifft((wavf + mfk(iz,:)),ns,2);
                P1 = P0(1:nx).*exp(-1i*w*dz*du(iz,:));              % adj ss correction
                P1 = P1.*we;
                wavf = fft(P1,ns,2).*S0;                       %upward continuation
                if save_wavefield
                    wavf_full(iz,:)=wavf;
                end
        end
        
    end
    
end

return
