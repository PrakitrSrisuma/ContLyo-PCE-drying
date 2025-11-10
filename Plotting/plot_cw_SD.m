function plot_cw_SD(t,c,linestyle1)

SD = sqrt(c(2,:));
plot(t, SD, linestyle1{:});
ylabel({'SD, Concentration';'(kg water/kg solid)'})
xlabel('Time (h)')

end