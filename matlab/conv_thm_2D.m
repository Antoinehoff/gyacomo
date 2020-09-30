function [ conv ] = conv_thm_2D( F, G, Pad )
%conv_thm_2D computes the convolution between F and G with a zero pading Pad
% kspace -> real -> product -> kspace

[Nr, Nz] = size(F);

f  = ifft2(F,Nr*Pad, Nz*Pad);
g  = ifft2(G,Nr*Pad, Nz*Pad);

conv_pad = fft2(f.*g); % convolution becomes product

conv = conv_pad(1:Nr,1:Nz); % remove padding
end

