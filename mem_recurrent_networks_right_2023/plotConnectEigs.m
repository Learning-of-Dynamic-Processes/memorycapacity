function [fig, eigs] = plotConnectEigs(A)
%

eigs = eig(A);

fig = figure();
tiledlayout(1, 1, 'TileSpacing', 'none')
nexttile
p1 = plot(eigs, '.', 'MarkerSize', 7);
hold on
p2 = fimplicit(@(x,y) x.^2 + y.^2 - 1, 'Color', [.5,.5,.5]);
hold off
axis equal
axis padded
grid
legend(p1, "eig(A)")

end

