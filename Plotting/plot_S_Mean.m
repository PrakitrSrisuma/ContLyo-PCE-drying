function plot_S_Mean(t,S,linestyle1)

mean = S(1,:)*100;
plot(t, mean, linestyle1{:});
ylabel({'Mean, Sublimation front (cm)'})
xlabel('Time (h)')

end