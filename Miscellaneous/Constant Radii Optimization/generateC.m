    function [C,uBasket,FBasket] = generateC(sel,r,NC,CA,A,E,C)
    % Initialize outputs
    uBasket = []; FBasket = [];
    
    % Iterate through once for each strain component
    for y = 1:1:3
    %   Define strain vector: [e11, e22, e12]'
        strainvec = [0;0;0];

    %   set that component equal to a dummy value (0.01 strain), 
    %       set all other values to zero
        strainvec(y) = 0.01; 
        strainvec(3) = strainvec(3)*2;

    %   use strain relations, BCs, and partitioned K-matrix to 
    %       solve for all unknowns
        e11 = strainvec(1); e22 = strainvec(2); e12 = strainvec(3);
        if (e11 ~= 0) || (e22 ~= 0) % when testing non-zero e11 or e22
            % Vector used to rearrange degrees of freedom
            qrvec = [4,7,9,10,11,16,1,2,3,5,6,8,12,13,14,15,17,18];
            % Rearranging K for rearranged degrees of freedom
            K = formK(NC,CA,A,E);
            newK = [K([4,7,9,10,11,16],:);K([1:3,5,6,8,12:15,17,18],:)];
            newK = [newK(:,[4,7,9,10,11,16]),...
                    newK(:,[1:3,5,6,8,12:15,17,18])];
            K = newK;
            % Defining known displacement, force boundary conditions
            u_r = [0;0;0;0;e22*sel;0;e22*sel;e11*sel;...
                   0;e11*sel;e11*sel;e22*sel];
            F_q = [0;0;0;0;0;0];
            % Partitioning K
            K_qq = K(1:6,1:6);
            K_rq = K(7:18,1:6);
            K_qr = K(1:6,7:18);
            K_rr = K(7:18,7:18);
            % Solving for unknown forces, displacements
            u_q = K_qq\(F_q-(K_qr*u_r)); 
            F_r = (K_rq*u_q)+(K_rr*u_r);
            % Assembling F,u vector in original order
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end
        else % when testing non-zero e12
            % Vector used to rearrange degrees of freedom
            qrvec = [4,9,10,16,1,2,3,5,6,7,8,11,12,13,14,15,17,18];
            % Rearranging K for rearranged degrees of freedom
            K = formK(NC,CA,A,E);            
            newK = [K([4,9,10,16],:);K([1:3,5:8,11:15,17,18],:)];
            newK = [newK(:,[4,9,10,16]),newK(:,[1:3,5:8,11:15,17,18])];
            K = newK;
            % Defining known displacement, force boundary conditions
            u_r = [0;0;0.5*e12*sel;e12*sel;0;0;0;e12*sel;0;0;...
                   0;0.5*e12*sel;e12*sel;0];
            F_q = [0;0;0;0];
            % Partitioning K
            K_qq = K(1:4,1:4);
            K_rq = K(5:18,1:4);
            K_qr = K(1:4,5:18);
            K_rr = K(5:18,5:18);
            % Solving for unknown forces, displacements
            u_q = K_qq\(F_q-(K_qr*u_r));
            F_r = (K_rq*u_q)+(K_rr*u_r);
            % Assembling F,u vector in original order
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end 
        end
        
    %   use F vector and geometry to solve for stress vector [s11,s22,s12]'
        F_x = F(13)+F(15)+F(17);
        F_y = F(6)+F(12)+F(18);
        F_xy = F(5)+F(11)+F(17);
        stressvec = (1/(sel*2*r)).*[F_x;F_y;F_xy];

    %   use strain and stress vectors to solve for the corresponding
    %       row of the C matrix
        Cdummy = stressvec/strainvec;
        C(:,y) = Cdummy(:,y);
        FBasket(:,y) = F;
        uBasket(:,y) = u;
    end
    
    %%% Write infeasible design CA into file
    %K_f = formK(NC,CA,A,E);
    %inf_count_new = inf_count;
    %CA_infeasible_up = CA_infeasible;
    %if (det(K_f) == 0)
        %inf_count = inf_count + 1;
        %CA_infeasible{inf_count} = CA;
    %end
    
    %disp(CA)
    
end
