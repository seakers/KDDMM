% FUNCTION TO CALCULATE VOLUME FRACTION 
function volFrac = calcVF(NC,CA,rvar,sel,sidenum)
    totalTrussVol = 0;
    for i = 1:size(CA,1)
        % Finding element length from nodal coordinates
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        L = sqrt(((x2-x1)^2)+((y2-y1)^2));
        % Adding current element to total volume of trusses
        totalTrussVol = totalTrussVol + (L*pi*(rvar(i)^2));
    end
    
    % Finding average side "thickness" due to differing element radii
    horizrads = [];
    for i = 1:1:size(CA,1)
        if ((CA(i,1) + sidenum) == CA(i,2)) && (NC(CA(i,1),2) == sel)
            horizrads = [horizrads,rvar(i)];
        end
    end
    vertrads = [];
    for i = 1:1:size(CA,1)
        if ((CA(i,1) + 1) == CA(i,2)) && (NC(CA(i,1),1) == sel)
            vertrads = [vertrads,rvar(i)];
        end
    end
    thick = mean([mean(horizrads),mean(vertrads)]);
    
    % Calculating volume fraction (using a solid square with 2*(avg 
    %   thickness) as a baseline)
    volFrac = totalTrussVol/(2*thick*(sel^2)); 
end