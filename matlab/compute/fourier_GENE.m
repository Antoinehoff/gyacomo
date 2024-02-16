function [ field_c ] = fourier_GENE( field_r )

%fft, symmetric as we are using real data
spectrumKxYZ=nx*fft(fieldxyz,[],1);
field_c=ny*fft(spectrumKxYZ,[],2);
clear spectrumKxKyZ 
 
end

