clear
close all
sample_fraction=1;

load('interp_solution_database.mat');
n_models = length(result(:));
if sample_fraction==1
    indices = 1:n_models;
else
    indices = randperm(n_models, sample_fraction*n_models);
end
slope_ratios = zeros(1, n_models);
pct_changes = zeros(1, n_models);
M_values = zeros(4, n_models);
J_30_values = zeros(1, n_models);
J_90_values = zeros(1, n_models);
for i=indices
    values = sscanf(result(i).fea.FileName, ...
        'ac%f_aB%f_n%d_E%d_wrp_bpf_res.out');
    ac = values(1);
    aB = values(2);
    n = values(3);
    E = values(4);
    Phi = result(i).fea.Phi;

    CMOD = result(i).fea.CMOD;
    
    J_last = result(i).fea.Jtotal_Avg(:,end);
    J_90 = result(i).fea.Jtotal_Avg(end,:);
    J_30_last = interp1(Phi, J_last, 30);
    J_90_last = J_last(end);
    J_30_values(i) = J_30_last;
    J_90_values(i) = J_90_last;

    M_epfm_a = result(i).fea.M_epfm_a(:,end);
    M_epfm_b = result(i).fea.M_epfm_b(:,end);
    M_lefm_a = result(i).fea.M_lefm_a(:,end);
    M_lefm_b = result(i).fea.M_lefm_b(:,end);
    M_epfm_a_30 = interp1(Phi, M_epfm_a, 30);
    M_epfm_b_30 = interp1(Phi, M_epfm_b, 30);
    M_lefm_a_30 = interp1(Phi, M_lefm_a, 30);
    M_lefm_b_30 = interp1(Phi, M_lefm_b, 30);
    M_epfm_a_90 = M_epfm_a(end);
    M_epfm_b_90 = M_epfm_b(end);
    M_lefm_a_90 = M_lefm_a(end);
    M_lefm_b_90 = M_lefm_b(end);
    M_values(:, i) = [M_epfm_a_30; M_epfm_a_90; M_epfm_b_30; M_epfm_b_90];

    CMOD_i = linspace(0,max(CMOD),100);
    J_i = interp1(CMOD, J_90, CMOD_i, 'spline');
    
    last_20pct = find(CMOD_i>=0.8*max(CMOD_i));
    prev_20pct = find(CMOD_i>=0.6*max(CMOD_i) & CMOD_i<0.8*max(CMOD_i));
    remainder = find(CMOD_i<0.6*max(CMOD_i));

    slopes = diff(J_90)./diff(CMOD);
    slope_ratio = slopes(end)/slopes(1);
    slope_ratios(i) = slope_ratio;

    tbl = table(CMOD_i.', J_i.', 'VariableNames', {'CMOD', 'J'});
    lm1 = fitlm(tbl(last_20pct,:), 'linear');
    lm2 = fitlm(tbl(prev_20pct,:), 'linear');
    pct_change = (1-lm2.Coefficients.Estimate(2)/lm1.Coefficients.Estimate(2))*100;
    pct_changes(i) = pct_change;
end