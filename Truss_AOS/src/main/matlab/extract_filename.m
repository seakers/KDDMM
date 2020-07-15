function [full_filepath] = extract_filename(fibre_stiffness,eps_moea,feas_and_stab)
% This function outputs the full location of the csv data file, given the
% boolean values for fibre_stiffness and eps_moea and the case number for
% feas_and_stab

filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\';
if eps_moea
    fileloc = 'Epsilon MOEA Runs\\';
    if fibre_stiffness
        filename = 'Fibre Stiffness code run results\\EpsilonMOEA_emoea0_fibrestiffness.csv';
    else
        filename = 'Truss code run results\\EpsilonMOEA_emoea0_trussstiffness.csv';
    end
else 
    switch feas_and_stab
        case 0
            fileloc = 'AOS MOEA Runs\\';
            if fibre_stiffness
                filename = 'Fibre Stiffness code run results\\Feas and Stab False\\AOSMOEA_constraint_adaptive0_fibrestiffness.csv';
                credit_filename = 'Fibre Stiffness code run results\\Feas and Stab False\\constraint_adaptive0.credit';
                quality_filename = 'Fibre Stiffness code run results\\Feas and Stab False\\constraint_adaptive0.qual';
            else
                filename = 'Truss code run results\\Feas and Stab False\\AOSMOEA_constraint_adaptive0_trussstiffness.csv';
                credit_filename = 'Truss code run results\\Feas and Stab False\\constraint_adaptive0.credit';
                quality_filename = 'Truss code run results\\Feas and Stab False\\constraint_adaptive0.qual';
            end
        case 1
            fileloc = 'AOS MOEA Runs\\';
            if fibre_stiffness
                filename = 'Fibre Stiffness code run results\\Feas True\\AOSMOEA_constraint_adaptive0_fibrestiffness.csv';
                credit_filename = 'Fibre Stiffness code run results\\Feas True\\constraint_adaptive0.credit';
                quality_filename = 'Fibre Stiffness code run results\\Feas True\\constraint_adaptive0.qual';
                
            else
                filename = 'Truss code run results\\Feas True\\AOSMOEA_constraint_adaptive0_trussstiffness.csv';
                credit_filename = 'Truss code run results\\Feas True\\constraint_adaptive0.credit';
                quality_filename = 'Truss code run results\\Feas True\\constraint_adaptive0.qual';
            end
        case 2
            fileloc = 'AOS MOEA Runs\\';
            if fibre_stiffness
                filename = 'Fibre Stiffness code run results\\Feas and Stab True\\AOSMOEA_constraint_adaptive0_fibrestiffness.csv';
                credit_filename = 'Fibre Stiffness code run results\\Feas and Stab True\\constraint_adaptive0.credit';
                quality_filename = 'Fibre Stiffness code run results\\Feas and Stab True\\constraint_adaptive0.qual';
            else
                filename = 'Truss code run results\\Feas and Stab True\\AOSMOEA_constraint_adaptive0_trussstiffness.csv';
                credit_filename = 'Truss code run results\\Feas and Stab True\\constraint_adaptive0.credit';
                quality_filename = 'Truss code run results\\Feas and Stab True\\constraint_adaptive0.qual';
            end
    end
end
full_filepath = strcat(filepath,fileloc,filename);

end

