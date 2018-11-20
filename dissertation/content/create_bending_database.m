clear;
load('interp_solution_database');
result_bending = struct([]);
at_list = 0.2:0.2:0.8;
ac_list = 0.2:0.2:1.0;
n_list = [3 4 6 10 20];
E_list = [100 200 300 500 700 1000];
num_mat = 1;
thickness = 1;
nu = 0.3;
matfilepattern = 'bending_mat_files/bend_ac%.1f_at%.1f_*_E%04d_n%02d_wrp.mat';
for i = 1:length(at_list)
    for j = 1:length(ac_list)
        for k = 1:length(n_list)
            for l = 1:length(E_list)
                
                matfile = sprintf(matfilepattern, ac_list(j), ...
                    at_list(i), E_list(l), n_list(k));
                d = dir(matfile);
                load(fullfile(d.folder, d.name)); % loads BCvalue, FileName, J, base_se, cmod, moment_arm, phi, reaction
                if isempty(reaction)
                    warning('Missing reaction data for a/t=%.1f, a/c=%.1f, n=%g, E=%g\n',...
                        at_list(i), ac_list(j), n_list(k), E_list(l));
                else
                    a = at_list(i)*thickness;
                    c = a/ac_list(j);
                    width = max(10*c, 10*thickness); % total width
                    reaction = -4*reaction;
                    % WARP3D reaction is for quarter plate in bending.
                    % Reaction is for half the length of one of a pair of
                    % rollers. Each pair of rollers takes up a total force
                    % measured as P, equal to 4*(WARP3D reaction).
                    M = reaction.*moment_arm; % moment_arm is approximately (S_outer-S_inner)/2
                    I = (width*thickness^3)/12;
                    num_steps = size(J, 1);
                    num_crack_nodes = size(J, 2);
                    S_inner = width; % total distance
                    S_outer = 2*S_inner; % total distance
                    plate_length = max(2*width, 1.1*S_outer); % total length
                    phi = phi*180/pi;
                    M_ratio = repmat(M'./M(1),1,num_crack_nodes);
                    J_step1 = repmat(J(1,:),num_steps,1);
                    
                    % s.reac_force = reaction;
                    s.reac_force = M/((S_outer-S_inner)/2); % as would be measured by instrument
                    s.moment = (s.reac_force/2)*(S_outer-S_inner)/2;
                    s.S_bend = (s.moment*thickness/2)/I;
                    s.Jtotal_Avg = J';
                    s.Kavg = [];
                    s.TstressAvg = [];
                    s.Jel_Avg = [];
                    s.Jel_EPFM = (J_step1.*(M_ratio.^2))';
                    s.inp_exists = 'yes';
                    s.moment_flag = 'yes';
                    s.CharStress = 0;
                    s.half_CMOD = 0.5*cmod;
                    s.CMOD = cmod;
                    s.CMOD_node = 1;
                    s.NameString = sprintf('a/t =%.1f, a/c =%.1f',at_list(i),...
                        ac_list(j));
                    s.FileName = FileName;
                    s.analy_type = 'EPFM';
                    s.num_steps = num_steps;
                    s.a = a;
                    s.c = c;
                    s.width = width;
                    s.B = thickness;
                    s.S_inner = S_inner;
                    s.S_outer = S_outer;
                    s.length = plate_length;
                    s.num_mat = num_mat;
                    s.crack_E = E_list(l);
                    s.crack_nu = nu;
                    s.Phi = phi;
                    s.BCvalue = BCvalue;
                    s.BCstring = 'Traction';
                    s.base_E_fea = s.crack_E;
                    s.base_nu_fea = nu;
                    s.haz_E_fea = [];
                    s.haz_nu_fea = [];
                    s.weld_E_fea = [];
                    s.weld_nu_fea = [];
                    s.base_se_fea = base_se';
                    s.haz_se_fea = [];
                    s.weld_se_fea = [];
                    s.r_phi_a = [];
                    s.r_phi_b = [];
                    s.M_lefm_a = [];
                    s.M_lefm_b = [];
                    s.M_epfm_a = [];
                    s.M_epfm_b = [];
                    
                    result_bending(i, j, k, l).fea = s;
                    fprintf('Wrote entry for a/t=%.1f, a/c=%.1f, n=%g, E=%g\n',...
                        at_list(i), ac_list(j), n_list(k), E_list(l));
                end
            end
        end
    end
end
save('interp_solution_database', 'input', 'result', 'result_bending');
