function plot_S_MeanCI(t,S,S_low,S_high,linestyle1,linestyle2)

mean = S;

curve1 = S_high;
curve2 = S_low;
a = [t, fliplr(t)];
area = [curve1, fliplr(curve2)];
plot(t, mean, linestyle1{:});
hold on; 
fill(a, area, linestyle2{:});
ylabel({'Sublimation front (cm)'})
xlabel('Time (h)')

end