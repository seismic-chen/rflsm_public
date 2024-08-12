%% Modified from Dr. Jeff Gu's frequency domain timeshift code %%
% Yunfeng Chen, Jul 3, 2013.
function datashift = fftShift(data,tax,tdiff)
x = tax;
y = data;
dt = x(2) - x(1);
% calculate frequency
nzeros = floor(abs(tdiff)/dt);
% computes all frequecies
if tdiff > 0
    y = [y;zeros(nzeros,1)];
else
    y = [zeros(nzeros,1);y];
end

n = length(y);
df=1/(dt*n);
spec = fft(y, n);
z = spec;
wtau = [1:n/2+1];


%yf   z_new = spec;
for i=1:n/2+1
    zreal=real(z(i));
    zimag=imag(z(i));
    wtau(i)=2.0*pi*(i-1)*df*tdiff;
    zr=zreal*cos(wtau(i))+zimag*sin(wtau(i));
    zi=zimag*cos(wtau(i))-zreal*sin(wtau(i));
    %yf     z_new(i) = spec(i)*exp(-1i*wtau(i));
    z(i)=complex(zr,zi);
    if(i ~= 1)
        z(n-i+2)=complex(zr,-zi);
        %yf         z_new(end-i+2) = conj(z_new(i));
    end
end

ispec=ifft(z, n);
if tdiff > 0
    ispecnew = ispec(1:length(x));
else
    ispecnew = ispec((end-length(x)+1):end);
end
length(ispecnew);
datashift = real(ispecnew);
%    plot(x, data, x, real(datashift),'r')
