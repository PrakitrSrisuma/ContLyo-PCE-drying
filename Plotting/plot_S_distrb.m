function plot_S_distrb(S,Nhis,color)

switch color
case '1a'
cc = [0 148 0]/256;
fa = .5;

case '1b'
cc = [173,255,17]/256;
fa = .5;


end

h = histogram(S,Nhis,'Normalization','probability','FaceColor',cc,'FaceAlpha',fa); hold on

ylabel('Probability')
xlabel('Sublimation front (cm)')
ax = ancestor(h, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.2f')

end