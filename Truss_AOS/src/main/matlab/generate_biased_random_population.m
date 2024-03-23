function [bool_des_run,obj,obj_pen,constr,heur] = generate_biased_random_population(pop_size,model,truss_problem,gen_mode,sidenum,sel,r,E,c_ratio,collapsibilityBiasFac,biasFactor)
% Generate a biased population based on the problem type
% Truss problem -> bias towards low satisfaction of partColl and conn
% Artery problem -> bias towards low satisfaction of conn

n_total_members = nchoosek(sidenum^2,2);
n_repeated_members = 2*nchoosek(sidenum,2);
n_variables = n_total_members - n_repeated_members;

NC = generateNC(sel,sidenum);
CA_all = get_CA_all(sidenum);

member_prob = 0.2; % for conn
member_prob_partcoll = 0.2; % for partcoll
partcoll_frac = 0.25; % fraction of designs biased against partial collapsibility heuristic
conn_frac = 0.25; % fraction of designs biased against connectivity constraint
n_des = 1;

bool_des_run = zeros(n_total_members,n_des);
bool_des_map = containers.Map;

constr = zeros(pop_size,3); %[feas, conn, stiffrat]
if ~truss_problem
    constr = zeros(pop_size,2); %[feas, conn]
end
obj = zeros(pop_size,2);
obj_pen = zeros(pop_size,2);
heur = zeros(pop_size,4); %[partcoll, nodalprop, orient, inters]

switch gen_mode
    case "Random" % No population biasing
        while n_des <= pop_size
            bool_des = randi([0,1],n_variables,1);
            complete_bool_des = get_complete_boolean_array(bool_des, sidenum);
            
            CA_des = CA_all(complete_bool_des~=0,:);
            rvar_des = r.*ones(1,size(CA_des,1));
            switch model
                case "Fibre"
                    if truss_problem
                        [C11,C22,volfrac_des] = fiberStiffnessModel_rVar_V3(sel,rvar_des,E,CA_des,nucFac,sidenum);
                        C_des = zeros(6);
                        C_des(1,1) = C11;
                        C_des(2,2) = C22;
                    else
                        disp("Fiber stiffness model not suitable for artery problem")
                        exit
                    end
                case "Truss"
                    [C_des, ~] = trussMetaCalc_NxN_1UC_rVar_AVar(sidenum,sel,rvar_des,E,CA_des);
                    volfrac_des = calcVF_NxN_feasOnly(CA_des,r,sel,sidenum);
                case "Beam"
                    C_des = Beam_2D_NxN_PBC(sel,sidenum,r,E,CA_des);
                    volfrac_des = calcVF_NxN_feasOnly(CA_des,r,sel,sidenum);
            end
            valid_design = true;
            if (any(isnan(C_des),'all'))
                valid_design = false;
            end
            if (any(C_des >= E,'all') || any(C_des <= 1,'all'))
                valid_design = false;
            end
            if ~valid_design
                continue
            end
            if (map_contains_design(bool_des_map,bool_des))
                continue
            else
                bool_des_run(:,n_des) = complete_bool_des;
                % Constraints
                constr(n_des,1) = feasibility_checker_nonbinary_V5(NC,CA_des,sel,sidenum);
                constr(n_des,2) = connectivityConstraint_PBC_2D(sidenum,NC,CA_des,sel,biasFactor);
                
                % True Objectives
                if truss_problem
                    obj(n_des,1) = C_des(2,2);
                    constr(n_des,3) = abs((C_des(2,2)/C_des(1,1)) - c_ratio);
                    obj(n_des,2) = volfrac_des;
                else
                    obj(n_des,1) = C_des(1,1)/volfrac_des;
                    obj(n_des,2) = (abs((C_des(2,2)/C_des(1,1)) - c_ratio) + abs((C_des(2,1)/C_des(1,1)) - 0.0745) + abs((C_des(1,2)/C_des(1,1)) - 0.0745) + abs(C_des(3,1)/E) + abs(C_des(1,3)/E) + abs(C_des(3,2)/E) + abs(C_des(2,3)/E) + abs((C_des(3,3)/C_des(1,1)) - 5.038))/8;
                end
                
                % Penalized Objectives
                if truss_problem
                    obj_pen(n_des,1) = -obj(n_des,1)./E;
                    obj_pen(n_des,2) = obj(n_des,2);
                else
                    obj_pen(n_des,1) = -(obj(n_des,1) - 2e5)./1e6;
                    obj_pen(n_des,2) = ((abs((C_des(2,2)/C_des(1,1)) - c_ratio))/6 + abs((C_des(2,1)/C_des(1,1)) - 0.0745) + abs((C_des(1,2)/C_des(1,1)) - 0.0745) + abs(C_des(3,1)/1.5e5) + abs(C_des(1,3)/9e4) + abs(C_des(3,2)/1.5e5) + abs(C_des(2,3)/9.5e4) + (abs((C_des(3,3)/C_des(1,1)) - 5.038) - 4.5)/0.5)/8;
                end
                
                % Heuristics
                heur(n_des,1) = partCollapseHeuristic_2D(sidenum,CA_des,NC,sel,collapsibilityBiasFac);
                heur(n_des,2) = connectivityHeuristic_2D(sidenum,NC,CA_des,sel,biasFactor);
                %[orient_des, ~] = orientationHeuristic(NC,CA_des,c_ratio);
                [orient_norm_des, ~] = orientationHeuristic_V2(NC,CA_des,c_ratio);
                %orient_run(des_count) = orient_des;
                heur(n_des,3) = orient_norm_des;
                %orient_avg_run(des_count) = avg_angle_des;
                %orient_avg_norm_run(des_count) = avg_angle_norm_des;
                heur(n_des,4) = intersectHeuristic(NC,CA_des);
                
                bool_des_map(strcat('des',num2str(n_des))) = bool_des_run;
                n_des = n_des + 1;
            end
            
        end
        
    case "PartcollAndConn" % Bias against both partial collapsibility and connectivity
        is_diag_member = identify_diagonal_members(sidenum);
        is_repeated_edge_members = identify_repeated_edge_members(sidenum);
        n_des_partcoll = 1;
        n_des_conn = 1;
        while n_des_partcoll <= ceil(partcoll_frac*pop_size)
            full_bool_des_pc = zeros(n_total_members,1);
            for i = 1:n_total_members
                if is_diag_member(i)
                    if rand < member_prob_partcoll
                        full_bool_des_pc(i) = 1;
                    end
                else
                    full_bool_des_pc(i) = randi([0 1],1,1);
                end
            end
            bool_des_pc = full_bool_des_pc(~is_repeated_edge_members);
            complete_bool_des_pc = get_complete_boolean_array(bool_des_pc, sidenum);
            
            CA_des = CA_all(complete_bool_des_pc~=0,:);
            rvar_des = r.*ones(1,size(CA_des,1));
            switch model
                case "Fibre"
                    if truss_problem
                        [C11,C22,volfrac_des] = fiberStiffnessModel_rVar_V3(sel,rvar_des,E,CA_des,nucFac,sidenum);
                        C_des = zeros(6);
                        C_des(1,1) = C11;
                        C_des(2,2) = C22;
                    else
                        disp("Fiber stiffness model not suitable for artery problem")
                        exit
                    end
                case "Truss"
                    [C_des, volfrac_des] = trussMetaCalc_NxN_1UC_rVar_AVar(sidenum,sel,rvar_des,E,CA_des);
                case "Beam"
                    C_des = Beam_2D_NxN_PBC(sel,sidenum,r,E,CA_des);
                    volfrac_des = calcVF_NxN_feasOnly(CA_des,r,sel,sidenum);
            end
            valid_design = true;
            if (any(isnan(C_des),'all'))
                valid_design = false;
            end
            if (any(C_des >= E,'all') || any(C_des <= 1,'all'))
                valid_design = false;
            end
            if ~valid_design
                continue
            end
            if (map_contains_design(bool_des_map,bool_des_pc))
                continue
            else
                bool_des_run(:,n_des_partcoll) = complete_bool_des_pc;
                % Constraints
                constr(n_des_partcoll,1) = feasibility_checker_nonbinary_V5(NC,CA_des,sel,sidenum);
                constr(n_des_partcoll,2) = connectivityConstraint_PBC_2D(sidenum,NC,CA_des,sel,biasFactor);
                
                % True Objectives
                if truss_problem
                    obj(n_des_partcoll,1) = C_des(2,2);
                    constr(n_des_partcoll,3) = abs((C_des(2,2)/C_des(1,1)) - c_ratio);
                    obj(n_des_partcoll,2) = volfrac_des;
                else
                    obj(n_des_partcoll,1) = C_des(1,1)/volfrac_des;
                    obj(n_des_partcoll,2) = (abs((C_des(2,2)/C_des(1,1)) - c_ratio) + abs((C_des(2,1)/C_des(1,1)) - 0.0745) + abs((C_des(1,2)/C_des(1,1)) - 0.0745) + abs(C_des(3,1)/E) + abs(C_des(1,3)/E) + abs(C_des(3,2)/E) + abs(C_des(2,3)/E) + abs((C_des(3,3)/C_des(1,1)) - 5.038))/8;
                end
                
                % Penalized Objectives
                if truss_problem
                    obj_pen(n_des_partcoll,1) = -obj(n_des_partcoll,1)./E;
                    obj_pen(n_des_partcoll,2) = obj(n_des_partcoll,2);
                else
                    obj_pen(n_des_partcoll,1) = -(obj(n_des_partcoll,1) - 2e5)./1e6;
                    obj_pen(n_des_partcoll,2) = ((abs((C_des(2,2)/C_des(1,1)) - c_ratio))/6 + abs((C_des(2,1)/C_des(1,1)) - 0.0745) + abs((C_des(1,2)/C_des(1,1)) - 0.0745) + abs(C_des(3,1)/1.5e5) + abs(C_des(1,3)/9e4) + abs(C_des(3,2)/1.5e5) + abs(C_des(2,3)/9.5e4) + (abs((C_des(3,3)/C_des(1,1)) - 5.038) - 4.5)/0.5)/8;
                end
                
                % Heuristics
                heur(n_des_partcoll,1) = partCollapseHeuristic_2D(sidenum,CA_des,NC,sel,collapsibilityBiasFac);
                heur(n_des_partcoll,2) = connectivityHeuristic_2D(sidenum,NC,CA_des,sel,biasFactor);
                %[orient_des, ~] = orientationHeuristic(NC,CA_des,c_ratio);
                [orient_norm_des, ~] = orientationHeuristic_V2(NC,CA_des,c_ratio);
                %orient_run(n_des_partcoll) = orient_des;
                heur(n_des_partcoll,3) = orient_norm_des;
                %orient_avg_run(n_des_partcoll) = avg_angle_des;
                %orient_avg_norm_run(n_des_partcoll) = avg_angle_norm_des;
                heur(n_des_partcoll,4) = intersectHeuristic(NC,CA_des);
                
                bool_des_map(strcat('des',num2str(n_des_partcoll))) = bool_des_run;
                n_des_partcoll = n_des_partcoll + 1;
            end
        end
        n_des_partcoll = n_des_partcoll - 1;
        while n_des_conn <= ceil(conn_frac*pop_size)
            bool_des_conn = zeros(n_variables,1);
            for j = 1:n_variables
                if rand < member_prob
                    bool_des_conn(j) = 1;
                end
            end
            complete_bool_des_conn = get_complete_boolean_array(bool_des_conn, sidenum);
            
            CA_des = CA_all(complete_bool_des_conn~=0,:);
            rvar_des = r.*ones(1,size(CA_des,1));
            switch model
                case "Fibre"
                    if truss_problem
                        [C11,C22,volfrac_des] = fiberStiffnessModel_rVar_V3(sel,rvar_des,E,CA_des,nucFac,sidenum);
                        C_des = zeros(6);
                        C_des(1,1) = C11;
                        C_des(2,2) = C22;
                    else
                        disp("Fiber stiffness model not suitable for artery problem")
                        exit
                    end
                case "Truss"
                    [C_des, volfrac_des] = trussMetaCalc_NxN_1UC_rVar_AVar(sidenum,sel,rvar_des,E,CA_des);
                case "Beam"
                    C_des = Beam_2D_NxN_PBC(sel,sidenum,r,E,CA_des);
                    volfrac_des = calcVF_NxN_feasOnly(CA_des,r,sel,sidenum);
            end
            valid_design = true;
            if (any(isnan(C_des),'all'))
                valid_design = false;
            end
            if (any(C_des >= E,'all') || any(C_des <= 1,'all'))
                valid_design = false;
            end
            if ~valid_design
                continue
            end
            if (map_contains_design(bool_des_map,bool_des_conn))
                continue
            else
                bool_des_run(:,n_des_partcoll+n_des_conn) = complete_bool_des_conn;
                % Constraints
                constr(n_des_partcoll+n_des_conn,1) = feasibility_checker_nonbinary_V5(NC,CA_des,sel,sidenum);
                constr(n_des_partcoll+n_des_conn,2) = connectivityConstraint_PBC_2D(sidenum,NC,CA_des,sel,biasFactor);
                
                % True Objectives
                if truss_problem
                    obj(n_des_partcoll+n_des_conn,1) = C_des(2,2);
                    constr(n_des_partcoll+n_des_conn,3) = abs((C_des(2,2)/C_des(1,1)) - c_ratio);
                    obj(n_des_partcoll+n_des_conn,2) = volfrac_des;
                else
                    obj(n_des_partcoll+n_des_conn,1) = C_des(1,1)/volfrac_des;
                    obj(n_des_partcoll+n_des_conn,2) = (abs((C_des(2,2)/C_des(1,1)) - c_ratio) + abs((C_des(2,1)/C_des(1,1)) - 0.0745) + abs((C_des(1,2)/C_des(1,1)) - 0.0745) + abs(C_des(3,1)/E) + abs(C_des(1,3)/E) + abs(C_des(3,2)/E) + abs(C_des(2,3)/E) + abs((C_des(3,3)/C_des(1,1)) - 5.038))/8;
                end
                
                % Penalized Objectives
                if truss_problem
                    obj_pen(n_des_partcoll+n_des_conn,1) = -obj(n_des_partcoll+n_des_conn,1)./E;
                    obj_pen(n_des_partcoll+n_des_conn,2) = obj(n_des_partcoll+n_des_conn,2);
                else
                    obj_pen(n_des_partcoll+n_des_conn,1) = -(obj(n_des_partcoll+n_des_conn,1) - 2e5)./1e6;
                    obj_pen(n_des_partcoll+n_des_conn,2) = ((abs((C_des(2,2)/C_des(1,1)) - c_ratio))/6 + abs((C_des(2,1)/C_des(1,1)) - 0.0745) + abs((C_des(1,2)/C_des(1,1)) - 0.0745) + abs(C_des(3,1)/1.5e5) + abs(C_des(1,3)/9e4) + abs(C_des(3,2)/1.5e5) + abs(C_des(2,3)/9.5e4) + (abs((C_des(3,3)/C_des(1,1)) - 5.038) - 4.5)/0.5)/8;
                end
                
                % Heuristics
                heur(n_des_partcoll+n_des_conn,1) = partCollapseHeuristic_2D(sidenum,CA_des,NC,sel,collapsibilityBiasFac);
                heur(n_des_partcoll+n_des_conn,2) = connectivityHeuristic_2D(sidenum,NC,CA_des,sel,biasFactor);
                %[orient_des, ~] = orientationHeuristic(NC,CA_des,c_ratio);
                [orient_norm_des, ~] = orientationHeuristic_V2(NC,CA_des,c_ratio);
                %orient_run(n_des_conn) = orient_des;
                heur(n_des_partcoll+n_des_conn,3) = orient_norm_des;
                %orient_avg_run(n_des_partcoll+n_des_conn) = avg_angle_des;
                %orient_avg_norm_run(n_des_partcoll+n_des_conn) = avg_angle_norm_des;
                heur(n_des_partcoll+n_des_conn,4) = intersectHeuristic(NC,CA_des);
                
                bool_des_map(strcat('des',num2str(n_des_partcoll+n_des_conn))) = bool_des_run;
                n_des_conn = n_des_conn + 1;
            end
        end
        n_des_conn = n_des_conn - 1;
        while n_des <= (pop_size - ceil(partcoll_frac*pop_size) - ceil(conn_frac*pop_size))
            bool_des = randi([0,1],n_variables,1);
            complete_bool_des = get_complete_boolean_array(bool_des, sidenum);
            
            CA_des = CA_all(complete_bool_des~=0,:);
            rvar_des = r.*ones(1,size(CA_des,1));
            switch model
                case "Fibre"
                    if truss_problem
                        [C11,C22,volfrac_des] = fiberStiffnessModel_rVar_V3(sel,rvar_des,E,CA_des,nucFac,sidenum);
                        C_des = zeros(6);
                        C_des(1,1) = C11;
                        C_des(2,2) = C22;
                    else
                        disp("Fiber stiffness model not suitable for artery problem")
                        exit
                    end
                case "Truss"
                    [C_des, volfrac_des] = trussMetaCalc_NxN_1UC_rVar_AVar(sidenum,sel,rvar_des,E,CA_des);
                case "Beam"
                    C_des = Beam_2D_NxN_PBC(sel,sidenum,r,E,CA_des);
                    volfrac_des = calcVF_NxN_feasOnly(CA_des,r,sel,sidenum);
            end
            valid_design = true;
            if (any(isnan(C_des),'all'))
                valid_design = false;
            end
            if (any(C_des >= E,'all') || any(C_des <= 1,'all'))
                valid_design = false;
            end
            if ~valid_design
                continue
            end
            if (map_contains_design(bool_des_map,bool_des))
                continue
            else
                bool_des_run(:,n_des_partcoll+n_des_conn+n_des) = complete_bool_des;
                % Constraints
                constr(n_des_partcoll+n_des_conn+n_des,1) = feasibility_checker_nonbinary_V5(NC,CA_des,sel,sidenum);
                constr(n_des_partcoll+n_des_conn+n_des,2) = connectivityConstraint_PBC_2D(sidenum,NC,CA_des,sel,biasFactor);
                
                % True Objectives
                if truss_problem
                    obj(n_des_partcoll+n_des_conn+n_des,1) = C_des(2,2);
                    constr(n_des_partcoll+n_des_conn+n_des,3) = abs((C_des(2,2)/C_des(1,1)) - c_ratio);
                    obj(n_des_partcoll+n_des_conn+n_des,2) = volfrac_des;
                else
                    obj(n_des_partcoll+n_des_conn+n_des,1) = C_des(1,1)/volfrac_des;
                    obj(n_des_partcoll+n_des_conn+n_des,2) = (abs((C_des(2,2)/C_des(1,1)) - c_ratio) + abs((C_des(2,1)/C_des(1,1)) - 0.0745) + abs((C_des(1,2)/C_des(1,1)) - 0.0745) + abs(C_des(3,1)/E) + abs(C_des(1,3)/E) + abs(C_des(3,2)/E) + abs(C_des(2,3)/E) + abs((C_des(3,3)/C_des(1,1)) - 5.038))/8;
                end
                
                % Penalized Objectives
                if truss_problem
                    obj_pen(n_des_partcoll+n_des_conn+n_des,1) = -obj(n_des_partcoll+n_des_conn+n_des,1)./E;
                    obj_pen(n_des_partcoll+n_des_conn+n_des,2) = obj(n_des_partcoll+n_des_conn+n_des,2);
                else
                    obj_pen(n_des_partcoll+n_des_conn+n_des,1) = -(obj(n_des_partcoll+n_des_conn+n_des,1) - 2e5)./1e6;
                    obj_pen(n_des_partcoll+n_des_conn+n_des,2) = ((abs((C_des(2,2)/C_des(1,1)) - c_ratio))/6 + abs((C_des(2,1)/C_des(1,1)) - 0.0745) + abs((C_des(1,2)/C_des(1,1)) - 0.0745) + abs(C_des(3,1)/1.5e5) + abs(C_des(1,3)/9e4) + abs(C_des(3,2)/1.5e5) + abs(C_des(2,3)/9.5e4) + (abs((C_des(3,3)/C_des(1,1)) - 5.038) - 4.5)/0.5)/8;
                end
                
                % Heuristics
                heur(n_des_partcoll+n_des_conn+n_des,1) = partCollapseHeuristic_2D(sidenum,CA_des,NC,sel,collapsibilityBiasFac);
                heur(n_des_partcoll+n_des_conn+n_des,2) = connectivityHeuristic_2D(sidenum,NC,CA_des,sel,biasFactor);
                %[orient_des, ~] = orientationHeuristic(NC,CA_des,c_ratio);
                [orient_norm_des, ~] = orientationHeuristic_V2(NC,CA_des,c_ratio);
                %orient_run(n_des_partcoll+n_des_conn+n_des) = orient_des;
                heur(n_des_partcoll+n_des_conn+n_des,3) = orient_norm_des;
                %orient_avg_run(n_des_partcoll+n_des_conn+n_des) = avg_angle_des;
                %orient_avg_norm_run(n_des_partcoll+n_des_conn+n_des) = avg_angle_norm_des;
                heur(n_des_partcoll+n_des_conn+n_des,4) = intersectHeuristic(NC,CA_des);
                
                bool_des_map(strcat('des',num2str(n_des_partcoll+n_des_conn+n_des))) = bool_des_run;
                n_des = n_des + 1;
            end
        end
        
    case "Conn" % Only bias against connectivity
        n_des_conn = 1;
        while n_des_conn <= ceil(conn_frac*pop_size)
            bool_des_conn = zeros(n_variables,1);
            for j = 1:n_variables
                if rand < member_prob
                    bool_des_conn(j) = 1;
                end
            end
            complete_bool_des_conn = get_complete_boolean_array(bool_des_conn, sidenum);
            
            CA_des = CA_all(complete_bool_des_conn~=0,:);
            rvar_des = r.*ones(1,size(CA_des,1));
            switch model
                case "Fibre"
                    if truss_problem
                        [C11,C22,volfrac_des] = fiberStiffnessModel_rVar_V3(sel,rvar_des,E,CA_des,nucFac,sidenum);
                        C_des = zeros(6);
                        C_des(1,1) = C11;
                        C_des(2,2) = C22;
                    else
                        disp("Fiber stiffness model not suitable for artery problem")
                        exit
                    end
                case "Truss"
                    [C_des, volfrac_des] = trussMetaCalc_NxN_1UC_rVar_AVar(sidenum,sel,rvar_des,E,CA_des);
                case "Beam"
                    C_des = Beam_2D_NxN_PBC(sel,sidenum,r,E,CA_des);
                    volfrac_des = calcVF_NxN_feasOnly(CA_des,r,sel,sidenum);
            end
            valid_design = true;
            if (any(isnan(C_des),'all'))
                valid_design = false;
            end
            if (any(C_des >= E,'all') || any(C_des <= 1,'all'))
                valid_design = false;
            end
            if ~valid_design
                continue
            end
            if (map_contains_design(bool_des_map,bool_des_conn))
                continue
            else
                bool_des_run(:,n_des_conn) = complete_bool_des_conn;
                % Constraints
                constr(n_des_conn,1) = feasibility_checker_nonbinary_V5(NC,CA_des,sel,sidenum);
                constr(n_des_conn,2) = connectivityConstraint_PBC_2D(sidenum,NC,CA_des,sel,biasFactor);
                
                % True Objectives
                if truss_problem
                    obj(n_des_conn,1) = C_des(2,2);
                    constr(n_des_conn,3) = abs((C_des(2,2)/C_des(1,1)) - c_ratio);
                    obj(n_des_conn,2) = volfrac_des;
                else
                    obj(n_des_conn,1) = C_des(1,1)/volfrac_des;
                    obj(n_des_conn,2) = (abs((C_des(2,2)/C_des(1,1)) - c_ratio) + abs((C_des(2,1)/C_des(1,1)) - 0.0745) + abs((C_des(1,2)/C_des(1,1)) - 0.0745) + abs(C_des(3,1)/E) + abs(C_des(1,3)/E) + abs(C_des(3,2)/E) + abs(C_des(2,3)/E) + abs((C_des(3,3)/C_des(1,1)) - 5.038))/8;
                end
                
                % Penalized Objectives
                if truss_problem
                    obj_pen(n_des_conn,1) = -obj(n_des_conn,1)./E;
                    obj_pen(n_des_conn,2) = obj(n_des_conn,2);
                else
                    obj_pen(n_des_conn,1) = -(obj(n_des_conn,1) - 2e5)./1e6;
                    obj_pen(n_des_conn,2) = ((abs((C_des(2,2)/C_des(1,1)) - c_ratio))/6 + abs((C_des(2,1)/C_des(1,1)) - 0.0745) + abs((C_des(1,2)/C_des(1,1)) - 0.0745) + abs(C_des(3,1)/1.5e5) + abs(C_des(1,3)/9e4) + abs(C_des(3,2)/1.5e5) + abs(C_des(2,3)/9.5e4) + (abs((C_des(3,3)/C_des(1,1)) - 5.038) - 4.5)/0.5)/8;
                end
                
                % Heuristics
                heur(n_des_conn,1) = partCollapseHeuristic_2D(sidenum,CA_des,NC,sel,collapsibilityBiasFac);
                heur(n_des_conn,2) = connectivityHeuristic_2D(sidenum,NC,CA_des,sel,biasFactor);
                %[orient_des, ~] = orientationHeuristic(NC,CA_des,c_ratio);
                [orient_norm_des, ~] = orientationHeuristic_V2(NC,CA_des,c_ratio);
                %orient_run(n_des_conn) = orient_des;
                heur(n_des_conn,3) = orient_norm_des;
                %orient_avg_run(n_des_conn) = avg_angle_des;
                %orient_avg_norm_run(n_des_conn) = avg_angle_norm_des;
                heur(n_des_conn,4) = intersectHeuristic(NC,CA_des);
                
                bool_des_map(strcat('des',num2str(n_des_conn))) = bool_des_run;
                n_des_conn = n_des_conn + 1;
            end
            
        end
        n_des_conn = n_des_conn - 1;
        while n_des <= (pop_size - ceil(conn_frac*pop_size))
            bool_des = randi([0,1],n_variables,1);
            complete_bool_des = get_complete_boolean_array(bool_des, sidenum);
            
            CA_des = CA_all(complete_bool_des~=0,:);
            rvar_des = r.*ones(1,size(CA_des,1));
            switch model
                case "Fibre"
                    if truss_problem
                        [C11,C22,volfrac_des] = fiberStiffnessModel_rVar_V3(sel,rvar_des,E,CA_des,nucFac,sidenum);
                        C_des = zeros(6);
                        C_des(1,1) = C11;
                        C_des(2,2) = C22;
                    else
                        disp("Fiber stiffness model not suitable for artery problem")
                        exit
                    end
                case "Truss"
                    [C_des, volfrac_des] = trussMetaCalc_NxN_1UC_rVar_AVar(sidenum,sel,rvar_des,E,CA_des);
                case "Beam"
                    C_des = Beam_2D_NxN_PBC(sel,sidenum,r,E,CA_des);
                    volfrac_des = calcVF_NxN_feasOnly(CA_des,r,sel,sidenum);
            end
            valid_design = true;
            if (any(isnan(C_des),'all'))
                valid_design = false;
            end
            if (any(C_des >= E,'all') || any(C_des <= 1,'all'))
                valid_design = false;
            end
            if ~valid_design
                continue
            end
            if (map_contains_design(bool_des_map,bool_des))
                continue
            else
                bool_des_run(:,n_des_conn+n_des) = complete_bool_des;
                % Constraints
                constr(n_des_conn+n_des,1) = feasibility_checker_nonbinary_V5(NC,CA_des,sel,sidenum);
                constr(n_des_conn+n_des,2) = connectivityConstraint_PBC_2D(sidenum,NC,CA_des,sel,biasFactor);
                
                % True Objectives
                if truss_problem
                    obj(n_des_conn+n_des,1) = C_des(2,2);
                    constr(n_des_conn+n_des,3) = abs((C_des(2,2)/C_des(1,1)) - c_ratio);
                    obj(n_des_conn+n_des,2) = volfrac_des;
                else
                    obj(n_des_conn+n_des,1) = C_des(1,1)/volfrac_des;
                    obj(n_des_conn+n_des,2) = (abs((C_des(2,2)/C_des(1,1)) - c_ratio) + abs((C_des(2,1)/C_des(1,1)) - 0.0745) + abs((C_des(1,2)/C_des(1,1)) - 0.0745) + abs(C_des(3,1)/E) + abs(C_des(1,3)/E) + abs(C_des(3,2)/E) + abs(C_des(2,3)/E) + abs((C_des(3,3)/C_des(1,1)) - 5.038))/8;
                end
                
                % Penalized Objectives
                if truss_problem
                    obj_pen(n_des_conn+n_des,1) = -obj(n_des_conn+n_des,1)./E;
                    obj_pen(n_des_conn+n_des,2) = obj(n_des_conn+n_des,2);
                else
                    obj_pen(n_des_conn+n_des,1) = -(obj(n_des_conn+n_des,1) - 2e5)./1e6;
                    obj_pen(n_des_conn+n_des,2) = ((abs((C_des(2,2)/C_des(1,1)) - c_ratio))/6 + abs((C_des(2,1)/C_des(1,1)) - 0.0745) + abs((C_des(1,2)/C_des(1,1)) - 0.0745) + abs(C_des(3,1)/1.5e5) + abs(C_des(1,3)/9e4) + abs(C_des(3,2)/1.5e5) + abs(C_des(2,3)/9.5e4) + (abs((C_des(3,3)/C_des(1,1)) - 5.038) - 4.5)/0.5)/8;
                end
                
                % Heuristics
                heur(n_des_conn+n_des,1) = partCollapseHeuristic_2D(sidenum,CA_des,NC,sel,collapsibilityBiasFac);
                heur(n_des_conn+n_des,2) = connectivityHeuristic_2D(sidenum,NC,CA_des,sel,biasFactor);
                %[orient_des, ~] = orientationHeuristic(NC,CA_des,c_ratio);
                [orient_norm_des, ~] = orientationHeuristic_V2(NC,CA_des,c_ratio);
                %orient_run(n_des_conn+n_des) = orient_des;
                heur(n_des_conn+n_des,3) = orient_norm_des;
                %orient_avg_run(n_des_conn+n_des) = avg_angle_des;
                %orient_avg_norm_run(n_des_conn+n_des) = avg_angle_norm_des;
                heur(n_des_conn+n_des,4) = intersectHeuristic(NC,CA_des);
                
                bool_des_map(strcat('des',num2str(n_des_conn+n_des))) = bool_des_run;
                n_des = n_des + 1;
            end
        end
end

end

