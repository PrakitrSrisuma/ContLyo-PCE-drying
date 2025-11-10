function plot_cw_Mean(t,c,linestyle1)

mean = c(1,:);
plot(t, mean, linestyle1{:});
ylabel({'Mean, Concentration';'(kg water/kg solid)'})
xlabel('Time (h)')

end