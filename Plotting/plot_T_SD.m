function plot_T_SD(t,T,linestyle1)

SD = sqrt(T(2,:));
plot(t, SD, linestyle1{:});
ylabel('SD, Temperature (K)')
xlabel('Time (h)')

end