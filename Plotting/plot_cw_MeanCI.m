function plot_cw_MeanCI(t,c,c_low,c_high,linestyle1,linestyle2)

mean = c;

curve1 = c_high;
curve2 = c_low;
a = [t, fliplr(t)];
area = [curve1, fliplr(curve2)];
plot(t, mean, linestyle1{:});
hold on
fill(a, area, linestyle2{:});
% ylabel({'Concentration';'(kg water/kg solid)'})
ylabel({'Concentration (wt/wt)'})
xlabel('Time (h)')

end