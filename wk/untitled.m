N = 100;
L = 20;

A0 = 2;
N0 =  4;
K0 = 2*pi/L*N0;

f = @(x_) A0 * sin(K0*x_);


x = linspace(0,L, N);

dk= 2*pi/L;
k = dk*(-N/2:N/2-1);

F = zeros(1,N);

[~,ik0p] = min(abs(k-K0));
[~,ik0m] = min(abs(k+K0));

F(ik0p) = -A0/2 * 1i/dk;
F(ik0m) =  A0/2 * 1i/dk;


f_ = ifft(fftshift(F));
figure
plot(x,f(x)); hold on;
plot(x,N*dk*f_)

%%
figure
plot(k,imag(F))