function [] = visualize_truss_varrad(NC,CA,rad_array,max_rad,acc_rad_fac)
% This function provides a graphic visualization of the truss 
% represented by CA (variable radii)
figure
% Plot node positions
labels = {'1','2','3','4','5','6','7','8','9'};
for i = 1:size(NC,1)
    plot(NC(i,1),NC(i,2),'*r')
    hold on
    text(NC(i,1),NC(i,2),labels{i},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15)
    hold on
end
% Plot truss elements one-by-one (line width based on member radius)
% smallest radius members have 0.5 linewidth, largest radius members have
% 15 linewidth and everything in between is linearly interpolated
min_rad = acc_rad_fac*max_rad;
for i = 1:size(CA,1)
    % Finding Positions of truss element end points
    x1 = NC(CA(i,1),1); 
    y1 = NC(CA(i,1),2);
    x2 = NC(CA(i,2),1);
    y2 = NC(CA(i,2),2);
    % Plotting line between the two end points
    %x_val = linspace(x1,x2);
    %y_val = linspace(y1,y2);
    lw = (((15-0.5)/(max_rad - min_rad))*(rad_array(i) - min_rad)) + 0.5; % linewidth calculation
    plot([x1,x2],[y1,y2],'-b','LineWidth',lw)
    hold on
end
hold off
%xlabel('X-position')
%ylabel('Y-position')
end

