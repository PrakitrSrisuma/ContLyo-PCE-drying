function plot_cw_distrb(c,Nhis,color)

switch color
case '1a'
cc = [0 148 0]/256;
fa = .5;

case '1b'
cc = [173,255,17]/256;
fa = .5;

case '2a'
cc = [0 .2 .5];
fa = .5;


end

h = histogram(c,Nhis,'Normalization','probability','FaceColor',cc,'FaceAlpha',fa); hold on

ylabel('Probability')
xlabel('Concentration (wt/wt)')
ax = ancestor(h, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.3f')

end