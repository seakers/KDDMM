function [] = visualize_truss_fromx_3x3(NC,CA_all,x)
% This function provides a graphic visualization of the 3x3 truss 
% represented by x

% Obtain the CA matrix from x and CA_all
% x is a double vector, so converting to logical vector
%x_vec = x>0.5;
%x_vec = [0, 1, 0, 0, x(1:2), 1, x(3:4), 0, x(5:18), x(19:end), 0, 1, 0]; % for forcing truss elements in 1,3; 1,7 and
%7;9
%x_des = x; % general case
%x_des = [x(1), 0, x(2:3), 0, 0, 0, 0, x(4:7), 0, 0, 0, 0, x(8:9), 0, 0, 0, ...
        %x(10), 0, x(11:12), 0, x(13:16), 0, x(17:19), 0, x(20)]; % only adjacent node connections allowed
x_des = [x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), x(10), ...
        x(11), x(12), x(13), x(14), x(15), x(16), x(17), x(3), x(18), x(19), ...
        x(20), x(21), x(22), x(23), x(24), x(25), x(26), x(27), x(28), x(29), ...
        x(30), x(31), x(23), x(1), x(32), x(9)]; 
% repeatable designs
        
CA = CA_all(x_des~=0,:);

figure
% Plot node positions
labels = {'1','2','3','4','5','6','7','8','9'};
for i = 1:size(NC,1)
    plot(NC(i,1),NC(i,2),'*r')
    hold on
    text(NC(i,1),NC(i,2),labels{i},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15,'FontWeight','bold')
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
    plot([x1,x2],[y1,y2],'-b','LineWidth',2)
    hold on
end
hold off
xlabel('X-position','FontSize',15,'FontWeight','bold')
ylabel('Y-position','FontSize',15,'FontWeight','bold')
end

