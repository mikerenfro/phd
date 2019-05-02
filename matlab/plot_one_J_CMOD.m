function plot_one_J_CMOD(result, index)
plot(result(index).fea.CMOD, result(index).fea.Jtotal_Avg(end,:), 'x', ...
    result(index).fea.CMOD, result(index).fea.Jtotal_Avg(17,:), 'o')
title(result(index).fea.FileName,'interpreter', 'none')
xlabel('CMOD');
ylabel('J');
legend('\phi=90', '\phi=30')
grid on;
