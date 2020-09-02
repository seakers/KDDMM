function [full_filepath] = extract_filename(fibre_stiffness,eps_moea,feas_and_stab,run_num)
% This function outputs the full location of the csv data file, given the
% boolean values for fibre_stiffness and eps_moea and the case number for
% feas_and_stab

filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\';
if eps_moea
    fileloc = 'Epsilon MOEA Runs\\';
    if fibre_stiffness
        filename = strcat('Fibre Stiffness code run results\\EpsilonMOEA_emoea',num2str(run_num),'_fibrestiffness.csv');
    else
        filename = strcat('Truss code run results\\EpsilonMOEA_emoea',num2str(run_num),'_trussstiffness.csv');
    end
else 
    switch feas_and_stab
        case 0
            fileloc = 'AOS MOEA Runs\\Feas and Stab False\\';
            if fibre_stiffness
                filename = strcat('Fibre Stiffness code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_fibrestiffness.csv');
                credit_filename = strcat('Fibre Stiffness code run results\\constraint_adaptive',num2str(run_num),'.credit');
                quality_filename = strcat('Fibre Stiffness code run results\\constraint_adaptive',num2str(run_num),'.qual');
            else
                filename = strcat('Truss code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_trussstiffness.csv');
                credit_filename = strcat('Truss code run results\\constraint_adaptive',num2str(run_num),'.credit');
                quality_filename = strcat('Truss code run results\\constraint_adaptive',num2str(run_num),'.qual');
            end
        case 1
            fileloc = 'AOS MOEA Runs\\Feas True\\';
            if fibre_stiffness
                filename = strcat('Fibre Stiffness code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_fibrestiffness.csv');
                credit_filename = strcat('Fibre Stiffness code run results\\constraint_adaptive',num2str(run_num),'.credit');
                quality_filename = strcat('Fibre Stiffness code run results\\constraint_adaptive',num2str(run_num),'.qual');
                
            else
                filename = strcat('Truss code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_trussstiffness.csv');
                credit_filename = strcat('Truss code run results\\constraint_adaptive',num2str(run_num),'.credit');
                quality_filename = strcat('Truss code run results\\constraint_adaptive',num2str(run_num),'.qual');
            end
        case 2
            fileloc = 'AOS MOEA Runs\\Feas and Stab True\\';
            if fibre_stiffness
                filename = strcat('Fibre Stiffness code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_fibrestiffness.csv');
                credit_filename = strcat('Fibre Stiffness code run results\\constraint_adaptive',num2str(run_num),'.credit');
                quality_filename = strcat('Fibre Stiffness code run results\\constraint_adaptive',num2str(run_num),'.qual');
            else
                filename = strcat('Truss code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_trussstiffness.csv');
                credit_filename = strcat('Truss code run results\\constraint_adaptive',num2str(run_num),'.credit');
                quality_filename = strcat('Truss code run results\\constraint_adaptive',num2str(run_num),'.qual');
            end
    end
end
full_filepath = strcat(filepath,fileloc,filename);

end

