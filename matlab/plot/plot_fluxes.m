figure
subplot(211)
plot(data.Ts0D,data.PGAMMA_RI'*data.scale);
subplot(212)
hold on
if(data.KIN_E)
plot(data.Ts0D,data.HFLUX_XE*data.scale);
end
plot(data.Ts0D,data.HFLUX_X*data.scale);