clear;
load('interp_solution_database');
close all;
for i=1:600
    Phi = result(i).fea.Phi;
    M_epfm_a = result(i).fea.M_epfm_a(:,end);
    M_epfm_b = result(i).fea.M_epfm_b(:,end);
    M_lefm_a = result(i).fea.M_lefm_a(:,end);
    M_lefm_b = result(i).fea.M_lefm_b(:,end);
    M_epfm_a_30(i) = interp1(Phi, M_epfm_a, 30); %#ok<*SAGROW>
    M_epfm_b_30(i) = interp1(Phi, M_epfm_b, 30);
    M_lefm_a_30(i) = interp1(Phi, M_lefm_a, 30);
    M_lefm_b_30(i) = interp1(Phi, M_lefm_b, 30);
    M_epfm_a_90(i) = M_epfm_a(end);
    M_epfm_b_90(i) = M_epfm_b(end);
    M_lefm_a_90(i) = M_lefm_a(end);
    M_lefm_b_90(i) = M_lefm_b(end);
end
semilogy(1:600, [M_epfm_a_30; M_epfm_b_30; M_epfm_a_90; M_epfm_b_90], '+');
legend('a, \phi=30', 'b, \phi=30', 'a, \phi=90', 'b, \phi=90', 'location', 'northwest');
xlabel('Model index');
ylabel('M');
grid

figure;
M_min = min([M_epfm_a_30; M_epfm_b_30; M_epfm_a_90; M_epfm_b_90]);
% semilogy(1:600, [ min([M_epfm_a_30; M_epfm_b_30; M_epfm_a_90; M_epfm_b_90]);
%     max([M_epfm_a_30; M_epfm_b_30; M_epfm_a_90; M_epfm_b_90]) ], '+');
% legend('Min M value (a or b, \phi=30 or 90', 'Max M value (a or b, \phi=30 or 90', 'location', 'southeast');
semilogy(1:600, M_min, '+');
grid
hold on;
semilogy([1 600],[20 20], 'k');
xlabel('Model index');
ylabel('min(M)');
%axis([0 600 1 400]);

j=0;
for i=1:600
    if M_min(i)>50
        filename = result(i).fea.FileName;
        fprintf('%s: %f\n',filename, M_min(i));
        j=j+1;
    end
end

figure; hist(M_min,0:10:300)
xlabel('M');
ylabel('Number of models in range');