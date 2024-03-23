%% Testing
clear
clc

%% Parameter definitions (printable designs)
E = 1.8162e6; % Young's Modulus for polymeric material (example: 1.8162 MPa for SIL material)
sel = 10e-3; % Unit square side length (NOT individual truss length) (in m)
r = 250e-6; % Radius for cross-sectional area of (assumed circular) truss members (in m)
A = pi*(r^2); % Cross-sectional area of truss member
NC = sel.*[0,0;0,0.5;0,1;0.5,0;0.5,0.5;0.5,1;1,0;1,0.5;1,1]; 
%CA_all = [1,2; 1,3; 1,4; 1,5; 1,6; 1,7; 1,8; 1,9; 2,3; 2,4; 2,5; 2,6; 2,7; 2,8; 2,9; 3,4; 3,5; 3,6; 3,7; 3,8; 3,9; 4,5; 4,6; 4,7; 4,8; 4,9; 5,6; 5,7; 5,8; 5,9; 6,7; 6,8; 6,9; 7,8; 7,9; 8,9];
c_ratio = 0.421;
sidenum = 3;
CA_all = get_CA_all(sidenum);
nucFac = 1;
biasFactor = 1;
collapsibilityBiasFac = 0.5;

%% constanr radius
des_bool = '101000000010000011000100001000100100';
x_des_bool = zeros(36,1);
for i = 1:36
    x_des_bool(i) = str2num(des_bool(i));
end
CA_des = CA_all(x_des_bool~=0,:);
rvar = r.*ones(1,size(CA_des,1));

[C_des,volfrac_des] = trussMetaCalc_NxN_1UC_rVar_AVar(sidenum,sel,rvar,E,CA_des);
feas_des = feasibility_checker_nonbinary_V5(NC,CA_des,sel,sidenum)
conn_des = connectivityConstraint_PBC_2D(sidenum,NC,CA_des,sel,biasFactor)
stiffrat_des = abs(C_des(2,2)/C_des(1,1) - c_ratio);
partcoll_des = partCollapseHeuristic_2D(sidenum,CA_des,NC,sel,collapsibilityBiasFac)
nodalprop_des = connectivityHeuristic_2D(sidenum,NC,CA_des,sel,biasFactor)
orient_des = orientationHeuristic_V2(NC,CA_des,c_ratio)
vf_des = calcVF_NxN_feasOnly(CA_des,rvar,sel,sidenum)

obj1_art = ((C_des(1,1)/vf_des) - 2e5)/1e6
obj2_art = (abs((C_des(2,2)/C_des(1,1)) - c_ratio) + abs((C_des(1,2)/C_des(1,1)) - 0.0745) + abs((C_des(2,1)/C_des(1,1)) - 0.0745) + abs((C_des(3,3)/C_des(1,1)) - 5.038) + abs(C_des(1,3)/E) + abs(C_des(3,1)/E) + abs(C_des(2,3)/E) + abs(C_des(3,2)/E))/8

%% variable radii
CA_des = [1,3;1,4;1,5;1,8;1,9;2,3;2,5;2,6;2,7;3,5;3,6;3,8;4,6;4,7;5,6;6,7;6,8;6,9;7,9;8,9];
rvar = [1.390751998940366E-4,1.433558029369192E-4,2.036450039025971E-4,1.7311172285120586E-4,1.3482808637984263E-4,2.2915938132555047E-4,1.444245375660846E-4,...
    2.139872317690465E-4,2.179114424393426E-4,1.6684650790068785E-4,1.433558029369192E-4,2.105158236030546E-4,1.9639881630862405E-4,2.2288319144105962E-4,2.0682590781032818E-4,...
    1.5350597306906777E-4,1.9461516007848772E-4,2.2288319144105962E-4,1.390751998940366E-4,2.2915938132555047E-4];

[C_des,volfrac_des] = trussMetaCalc_NxN_rVar_AVar(nucFac,sidenum,sel,rvar,E,CA_des)
feas_des = feasibility_checker_nonbinary_V2(NC,CA_des)
conn_des = connectivityConstraint_NPBC_2D(sidenum,NC,CA_des,sel,biasFactor)
stiffrat_des = abs(C_des(2,2)/C_des(1,1) - c_ratio)
partcoll_des = partCollapseHeuristic_2D(sidenum,CA_des,NC,sel,collapsibilityBiasFac)
nodalprop_des = connectivityHeuristic_2D(sidenum,NC,CA_des,sel,biasFactor)
orient_des = orientationHeuristicNorm(NC,CA_des,sel,c_ratio)

%% functions
function CA_all = get_CA_all(sidenum)
    n_member = 0;
    n_total_members = nchoosek(sidenum^2,2);
    CA_all = zeros(n_total_members,2);
    for i = 1:(sidenum^2)
        for j = i+1:(sidenum^2)
            CA_all(n_member+1,:) = [i,j];
            n_member = n_member + 1;
        end
    end
end