function [] = visualize_truss_3x3_gui(NC,CA,axes_handle)
% This function provides a graphic visualization of the 3x3 truss 
% represented by CA

axes(axes_handle)

% Plot node positions
labels = {'1','2','3','4','5','6','7','8','9'};

for i = 1:size(NC,1)
    plot(NC(i,1),NC(i,2),'*r');
    hold on
    text(NC(i,1),NC(i,2),labels{i},'VerticalAlignment','bottom','HorizontalAlignment','right')
    hold on
end
% Plot truss elements one-by-one
for i = 1:size(CA,1)
    % Finding Positions of truss element end points
    x1 = NC(CA(i,1),1); 
    y1 = NC(CA(i,1),2);
    x2 = NC(CA(i,2),1);
    y2 = NC(CA(i,2),2);
    % Plotting line between the two end points
    %x_val = linspace(x1,x2);
    %y_val = linspace(y1,y2);
    plot([x1,x2],[y1,y2],'-b','LineWidth',2);
    hold on
end
hold off
%set(axes_handle,'CurrentAxes',gca);
set(get(axes_handle, 'xlabel'), 'string', 'X-position')
set(get(axes_handle, 'ylabel'), 'string', 'Y-position')

end

