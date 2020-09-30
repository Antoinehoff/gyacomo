function kernel_=kernel(n,x)
% Evaluate the kernel function of order n. The normalized argument of the kernel function 
% is b_a = k_\perp \sigma_a \sqrt{2 \tau_a}.
%
% Note: - if n<0, then return 0.

% Retrieve physical parameters

if n >= 0
    kernel_=(x./2).^(2*n).*exp(-x.^2/4)./factorial(n);
else 
    kernel_ =0;
end

end
