function plot_S_SD(t,c,linestyle1)

SD = sqrt(c(2,:))*100;
plot(t, SD, linestyle1{:});
ylabel({'SD, Sublimation front (cm)'})
xlabel('Time (h)')

end