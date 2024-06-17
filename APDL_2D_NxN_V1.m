% -----------------------------------------------------------------
% ANSYS APDL 2D 3x3 TRUSS CODE
% This code considers only 3x3 2D trusses
% Function for calculating metamaterial properties for a 2D NxN truss
% All lengths are in [m], all stresses and moduli are in [Pa] 
% Each node has two degrees of freedom (DOF), x and y.  The (node number*2) 
%    gives the number for that node's yDOF, and (yDOF - 1) is the xDOF
% -----------------------------------------------------------------
% This code uses MATLAB to call ANSYS APDL for FEM anlysis.
% First, MATLAB writes the design parameters into input file 
% (e.g. para_in.txt)
% Then, MATLAB calls ANSYS APDL to excute the APDL file 
% The APDL file reads the parameters from the input file and writes the 
% analysis result to the output file (e.g. para_out.txt)
% Finally, MATLAB reads the results from the output file
% -----------------------------------------------------------------
% CREDITS FOR ORIGINAL CODE PACKAGE:
% https://github.com/zhandawei/Call_ANSYS_in_MATLAB
% Dawei Zhan
% zhandawei{at}swjtu.edu.cn
% -----------------------------------------------------------------
% Modified/adapted for KDDM project by Srikar Srivatsa
% -----------------------------------------------------------------
% Sample values for Case 1 (for testing script):
%{
sel = 0.05; r = 50*(10^-6); E = 10000; sidenum = 3;
CA = [1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;1,5;2,5;3,5;4,5;5,6;5,7;5,8;5,9];
%}

function C = APDL_2D_NxN_V1(sel,sidenum,r,E,CA)
    C = [0,0,0;0,0,0;0,0,0];
    % Iterate through once for each strain component
    for y = 1:1:3
        % Define strain vector: [e11, e22, e12]'
        strainvec = [0;0;0];

        % Set that component equal to a dummy value (0.01 strain), set all 
        % other values to zero
        strainvec(y) = 0.01; 
        strainvec(3) = strainvec(3)*2;
        
        % Populate Input Vector
        qsel = sel./(10^-3); qr = r./(10^-6);
        input_vec=[qsel;qr;E;size(CA,1);y;sidenum;CA(:,1);CA(:,2)];

        % Define paths and input/output files (PATH INCORRECT)
        ANSYS_path=...
 'C:\Program Files\ANSYS Inc\v202\ansys\bin\winx64\ANSYS202';
        APDL_name='C:\Users\Public\APDL_Runs\Truss2D_NxN.txt';
        input_file_name='C:\Users\Public\APDL_Runs\para_in.txt';
        output_file_name='C:\Users\Public\APDL_Runs\para_out.txt';
        ANSYS_path=strcat('"',ANSYS_path,'"');
        disp('Current stage: ');disp(y);

        % Write input vector to para_in.txt
        %writematrix(input_vec,input_file_name);
        dlmwrite(input_file_name,input_vec,'delimiter',' ',...
           'precision','%8.3f','newline','pc');
        % Call ANSYS APDL, perform FEA
        status = system(sprintf('%s -b -p aa_r -i %s -o out.txt',...
            ANSYS_path,APDL_name));
        disp('Simulation status: ');disp(status);
        % Read the results from para_out.txt
        output_vec=load(output_file_name)';
        F_x = sum(output_vec(1:sidenum));
        F_y = sum(output_vec((sidenum+1):(2*sidenum)));
        F_xy = sum(output_vec(((2*sidenum)+1):(3*sidenum)));
        FBank(y,:) = [F_x,F_y,F_xy];

        % Calculate stress vector
        stressvec = (1/(sel*2*r)).*[F_x;F_y;F_xy];

        % Use strain and stress vectors to solve for the corresponding row 
        % of the C matrix
        Cdummy = stressvec/strainvec;
        C(:,y) = Cdummy(:,y);
    end
    C = C*10;
end






