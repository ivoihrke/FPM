function [ psi] = Proj_Fourier_v2( psi0, I, I0, c, F )
%PROJ_FOURIER projection based on intensity measurement in the fourier
%domain, replacing the amplitude of the Fourier transform by measured
%amplitude, sqrt(I)
% last modified by Lei Tian, lei_tian@alum.mit.edu, 3/1/2014
%
% psi0 - FT of low-res image with estimated phase, possible a 3D stack thereof for the multiplex case
% I    - measured image (real)
% I0   - estimated image (real) 
% c    - relative intensity/brightness of illumination LED (also an array of m values for the multiplex case)
% F    - function handle to FT


[n1,n2,r] = size(psi0);

if r == 1
%     psi = F(sqrt(I/c).*psi0./(sqrt(I0)+eps));  #this seems to be the old attempt
    %understand this one -- keep angle, but replace amplitude
    psi = F(sqrt(I/c).*exp(1i .* angle(psi0)));
else
    %the multiplex version uses the old definition, and applies it per LED (channel)
    psi = zeros(n1,n2,r);
    for m = 1:r
        psi(:,:,m) = F(sqrt(I/c(m)).*psi0(:,:,m)./sqrt(I0+eps));
    end
end

end

