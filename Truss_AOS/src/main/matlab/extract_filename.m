function [full_filepath] = extract_filename(fibre_stiffness,eps_moea,feas_and_stab,bias_init,run_num)
% This function outputs the full location of the csv data file, given the
% boolean values for fibre_stiffness and eps_moea and the case number for
% feas_and_stab

filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\';
if eps_moea
    fileloc = 'Epsilon MOEA Runs\\';
    if fibre_stiffness
        filename_model = 'Fibre Stiffness code run results\\';
        if bias_init
            filename = strcat('Biased Initialization\\EpsilonMOEA_emoea',num2str(run_num),'_biasedinit_fibrestiffness.csv');
        else
            filename = strcat('Random Initialization\\EpsilonMOEA_emoea',num2str(run_num),'_fibrestiffness.csv');
        end
    else
        filename_model = 'Truss code run results\\';
        if bias_init
            filename = strcat('Biased Initialization\\EpsilonMOEA_emoea',num2str(run_num),'_biasedinit_trussstiffness.csv');
        else
            filename = strcat('Random Initialization\\EpsilonMOEA_emoea',num2str(run_num),'_trussstiffness.csv');
        end
    end
else 
    switch feas_and_stab
        case 0
            fileloc = 'AOS MOEA Runs\\Feas and Stab False\\';
            if fibre_stiffness
                filename_model = strcat('Fibre Stiffness code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_fibrestiffness.csv');
            else
                filename_model = strcat('Truss code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_trussstiffness.csv');
            end
        case 1
            fileloc = 'AOS MOEA Runs\\Feas True\\';
            if fibre_stiffness
                filename_model = strcat('Fibre Stiffness code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_fibrestiffness.csv');
                
            else
                filename_model = strcat('Truss code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_trussstiffness.csv');
            end
        case 2
            fileloc = 'AOS MOEA Runs\\Feas and Stab True\\';
            if fibre_stiffness
                filename_model = strcat('Fibre Stiffness code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_fibrestiffness.csv');
            else
                filename_model = strcat('Truss code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_trussstiffness.csv');
            end
    end
    filename = '';
end
full_filepath = strcat(filepath,fileloc,filename_model,filename);

end

