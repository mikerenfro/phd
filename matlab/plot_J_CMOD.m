clear
close all
sample_fraction = 1;

load('interp_solution_database.mat');
n_models = length(result(:));
if sample_fraction==1
    indices = 1:n_models;
else
    indices = randperm(n_models, sample_fraction*n_models);
end
slope_ratios = zeros(1, n_models);
pct_changes = zeros(1, n_models);
for i=indices
%     figure(2);
    CMOD = result(i).fea.CMOD;
    CMOD = CMOD/max(CMOD);
    % J_criteria = mean(result(i).fea.Jtotal_Avg(result(i).fea.Phi>45,:));
    J_criteria = result(i).fea.Jtotal_Avg(end,:);
    J_criteria = J_criteria/max(J_criteria);

    CMOD_i = linspace(0,max(CMOD),100);
    J_i = interp1(CMOD, J_criteria, CMOD_i, 'spline');
    
    last_20pct = find(CMOD_i>=0.8*max(CMOD_i));
    prev_20pct = find(CMOD_i>=0.6*max(CMOD_i) & CMOD_i<0.8*max(CMOD_i));
    remainder = find(CMOD_i<0.6*max(CMOD_i));
    tbl = table(CMOD_i.', J_i.', 'VariableNames', {'CMOD', 'J'});
    lm1 = fitlm(tbl(last_20pct,:), 'linear');
    lm2 = fitlm(tbl(prev_20pct,:), 'linear');
    pct_change = (1-lm2.Coefficients.Estimate(2)/lm1.Coefficients.Estimate(2))*100;
    slopes = diff(J_i)./diff(CMOD_i);
    slope_ratio = slopes(end)/slopes(1);
    slope_ratios(i) = slope_ratio;
    pct_changes(i) = pct_change;
    fprintf('Model %d: slope change= %.1f%%, slope ratio= %.1f\n', i, pct_change, slope_ratio);

    plot([0 CMOD], [0 J_criteria], 'o');
    hold on;
    plot([0 CMOD_i(remainder)], [0 J_i(remainder)], ...
        CMOD_i(prev_20pct), J_i(prev_20pct), ...
        CMOD_i(last_20pct), J_i(last_20pct))
    title(sprintf('Model %d: %s', i, result(i).fea.FileName), 'interpreter', 'none');
    ylabel('|J_{Total}(\phi=90)|');
    xlabel('CMOD');
    
    plot(CMOD_i, lm2.Coefficients.Estimate(2)*CMOD_i+lm2.Coefficients.Estimate(1));
    plot(CMOD_i, lm1.Coefficients.Estimate(2)*CMOD_i+lm1.Coefficients.Estimate(1),'--');
    legend('Data', '0-60% of range', '61-80% of range', ...
        '81-100% of range', 'Linear fit (61-80% of range)', ...
        'Linear fit (81-100% of range)', 'Location', 'Southeast');
    title(sprintf('Model %d: %s, %.1f%% slope change, %.1f slope ratio', ...
        i, result(i).fea.FileName, pct_change, slope_ratio), ...
        'interpreter', 'none');
    if abs(pct_change)<5
        modelfilename = sprintf('model%03d.png',i);
        print('-dpng', modelfilename);
    end
    hold off;
%     pause;
end