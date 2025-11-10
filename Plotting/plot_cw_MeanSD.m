function plot_cw_MeanSD(t,c,linestyle1,linestyle2)

mean = c(1,:);
SD = sqrt(c(2,:));

curve1 = mean + SD;
curve2 = mean - SD;
a = [t, fliplr(t)];
area = [curve1, fliplr(curve2)];
plot(t, mean, linestyle1{:});
hold on
fill(a, area, linestyle2{:});
% ylabel({'Concentration';'(kg water/kg solid)'})
ylabel({'Concentration (wt/wt)'})
xlabel('Time (h)')

end