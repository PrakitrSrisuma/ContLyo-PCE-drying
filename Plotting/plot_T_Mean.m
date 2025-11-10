function plot_T_Mean(t,T,linestyle1)

mean = T(1,:);
plot(t, mean, linestyle1{:});
ylabel('Mean, Temperature (K)')
xlabel('Time (h)')

end