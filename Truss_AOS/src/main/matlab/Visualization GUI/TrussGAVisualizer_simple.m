function TrussGAVisualizer
% TrussGAvisualizer displays the pareto front, along with insights into
% specific designs in the pareto front for GA runs using either Epsilon
% MOEA or AOS, using the fibre stiffness or truss model. This app is
% specifically for the 2D 3x3 node matrix case with target stiffness ratio
% as 1.
    
    %  Create and then hide the UI as it is being constructed.
    f = figure('Visible','off','Position',[250,250,800,500]);
    
    % Set initial values for shared data
    setappdata(f,'pen_fac',0);
    setappdata(f,'run_number',0);
    setappdata(f,'feas_score',0.0);
    setappdata(f,'stab_score',0.0);
    setappdata(f,'true_objs',[0.0,0.0]);
    setappdata(f,'C_mat',zeros(3));
    setappdata(f,'filename','none');
    setappdata(f,'filename2','none');
    setappdata(f,'filepath2','none');
    setappdata(f,'filepath3','none');
    
    % Construct the components
    H.hEditRunNumText = uicontrol('Style','text','String','Enter Run Number (0-29)',...
                    'Position',[490,450,100,30]);
       
    H.hEditRunNum = uicontrol('Parent',f,'BackgroundColor',[1 1 1], 'Position',[490,425,100,20],...
                'Style','edit','Tag','run_number','Callback',@runNumberCallback);
    
    H.hAlgPopupText = uicontrol('Style','text','String','Algorithm Selection',...
           'Position',[30,450,200,20]);
    
    H.hAlgPopUp = uicontrol('Parent',f,'Style','popupmenu', ...
                'String',{'Select Algorithm','Epsilon MOEA','AOS - Constraint Independent','AOS - Feasibility Dependent','AOS - Constraint Dependent'},...
                'Position',[30,425,200,20],'Callback',@algorithmPopupCallback);
            
    H.hModelPopupText = uicontrol('Parent',f,'Style','text','String','Model Selection',...
                    'Position',[260,450,200,20]);
       
    H.hModelPopUp = uicontrol('Style','popupmenu',...
                'String',{'Select Model','Fibre Stiffness','Truss Stiffness'},...
                'Position',[260,425,200,20],'Callback',@modelPopupCallback);
            
    H.hButPlotPareto = uicontrol('Parent',f,'Style','pushbutton','FontSize',10,...
                'FontWeight','bold', 'Position',[620,425,175,20],...
                'String','Plot Pareto Front','Callback',@plotParetoFrontCallback);
            
    H.hParetoFigText = uicontrol('Style','text','String','Pareto Front',...
                   'FontWeight','bold', 'Position',[200,375,200,20]);
               
    H.hParetoFig = axes('Parent',f,'Units','Pixels','Position',[100,125,400,250],...
                    'Tag','pareto_axes');
    
    H.hTrussFigText = uicontrol('Style','text','String','Selected Design',...
                   'FontWeight','bold', 'Position',[565,375,200,20]);
    
    H.hTrussFig = axes('Parent',f,'Units','Pixels','Position',[575,200,175,175],...
                    'Tag','truss_axes');
    
    H.hFeasScoreText = uicontrol('Style','text','String','Feasibility Score',...
                   'FontWeight','bold', 'Position',[40,50,100,20]);
    
    H.hEditFeasScore = uicontrol('Parent',f,'BackgroundColor',[1 1 1],...
                    'Position',[40,30,100,20],'Style','edit','Tag','feas_score');
                
    H.hStabScoreText = uicontrol('Style','text','String','Stability Score',...
                   'FontWeight','bold', 'Position',[180,50,100,20]);
                
    H.hEditStabScore = uicontrol('Parent',f,'BackgroundColor',[1 1 1],...
                    'Position',[180,30,100,20], 'Style','edit','Tag','stab_score');
                
    H.hCMatText = uicontrol('Style','text','String','Truss Stiffness Matrix',...
                   'FontWeight','bold', 'Position',[560,130,200,20]);
                
    %H.hEditCMat = uicontrol('Parent',f,'BackgroundColor',[1 1 1],...
                  %'Position',[560,30,200,60],'Style','edit','Tag','C_mat');
                
    H.hEditCMat = uitable('Parent',f,'Position',[560,30,205,100],...
                    'Data',cell(3,3),'Tag','C_mat');
    
    H.hTrueObjsText = uicontrol('Style','text','String','True Objectives',...
                   'FontWeight','bold', 'Position',[320,50,200,20]);
                
    H.hEditTrueObjs = uicontrol('Parent',f,'BackgroundColor',[1 1 1],...
                    'Position',[320,30,200,20], 'Style','edit','Tag','true_objs');     
                
    % Change units to normalized so components resize automatically
    f.Units = 'normalized';
    H.hEditRunNumText.Units = 'normalized';
    H.hEditRunNum.Units = 'normalized';
    H.hAlgPopupText.Units = 'normalized';
    H.hAlgPopUp.Units = 'normalized';
    H.hModelPopupText.Units = 'normalized';
    H.hModelPopUp.Units = 'normalized';
    H.hButPlotPareto.Units = 'normalized';
    H.hParetoFigText.Units = 'normalized';
    H.hParetoFig.Units = 'normalized';
    H.hTrussFigText.Units = 'normalized';
    H.hTrussFig.Units = 'normalized';
    H.hFeasScoreText.Units = 'normalized';
    H.hEditFeasScore.Units = 'normalized';
    H.hStabScoreText.Units = 'normalized';
    H.hEditStabScore.Units = 'normalized';
    H.hCMatText.Units = 'normalized';
    H.hEditCMat.Units = 'normalized';
    H.hTrueObjsText.Units = 'normalized';
    H.hEditTrueObjs.Units = 'normalized';
    
    % Assign a name to appear on the figure title
    f.Name = 'Truss GA Visualizer';
    
    % Move the window to the centre of the screen
    movegui(f,'center')
                
    % Make the UI Visible
    f.Visible = 'on';
    
    % Run Number textbox callback function
    function runNumberCallback(source, eventdata)
        run_number = floor(str2double(get(source, 'String')));
        setappdata(source.Parent,'run_number',run_number);
    end

    % Algorithm popup menu callback function
    function algorithmPopupCallback(source, eventdata)
        % Determine the selected option
        str = get(source,'String');
        val = get(source,'val');
        algorithm_choice = str{val};
        % Set filepath2 based on the selected value
        switch algorithm_choice
            case 'Epsilon MOEA'
                filepath2 = 'Epsilon MOEA Runs\\';
                filename = 'EpsilonMOEA_emoea';
            case 'AOS - Constraint Independent'
                filepath2 = 'AOS MOEA Runs\\Feas and Stab False\\';
                filename = 'AOSMOEA_constraint_adaptive';
            case 'AOS - Feasibility Dependent'
                filepath2 = 'AOS MOEA Runs\\Feas True\\';
            case 'AOS - Constraint Dependent'
                filepath2 = 'AOS MOEA Runs\\Feas and Stab True\\';
        end
        setappdata(source.Parent,'filepath2',filepath2);
        setappdata(source.Parent,'filename',filename);
    end
    
    % Model popup menu callback function
    function modelPopupCallback(source, eventdata)
        % Determine the selected option
        str = get(source,'String');
        val = get(source,'val');
        modelChoice = str{val};
        % Set filepath3 based on the selected value
        switch modelChoice
            case 'Fibre Stiffness'
                filepath3 = 'Fibre Stiffness code run results\\';
                filename2 = '_fibrestiffness.csv';
                pen_fac = 1.5;
            case 'Truss Stiffness'
                filepath3 = 'Truss code run results\\';
                filename2 = '_trussstiffness.csv';
                pen_fac = 1;
        end
        setappdata(source.Parent,'filepath3',filepath3);
        setappdata(source.Parent,'filename2',filename2);
        setappdata(source.Parent,'pen_fac',pen_fac);
    end
    
    % Plot pareto front button callback
    function plotParetoFrontCallback(source, eventdata)
        filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\';
        
        run_number = getappdata(source.Parent,'run_number');
        filename = getappdata(source.Parent,'filename');
        filename2 = getappdata(source.Parent,'filename2');
        filepath2 = getappdata(source.Parent,'filepath2');
        filepath3 = getappdata(source.Parent,'filepath3');
        pen_fac = getappdata(source.Parent,'pen_fac');
        
        full_filename = strcat(filename,num2str(run_number),filename2);
        full_filepath = strcat(filepath,filepath2,filepath3,full_filename);
        data_table = readtable(full_filepath,'Format','%s%f%f%f%f','HeaderLines',1);
        pop_size =  size(data_table,1);
        csv_data = zeros(pop_size,4);
        designs = strings(pop_size);
        csv_data = data_table(:,2:end);
        designs = data_table(:,1);
        
        % Parameters
        E = 10000; % Young's Modulus for polymeric material (example: 10000 Pa)
        sel = 0.05; % Unit square side length (NOT individual truss length) (example: 5 cm)
        r = 50*(10^-6); % Radius for cross-sectional area of (assumed circular) truss members (example: 50 micrometers)
        A = pi*(r^2); % Cross-sectional area of truss member
        NC = sel.*[0,0;0,0.5;0,1;0.5,0;0.5,0.5;0.5,1;1,0;1,0.5;1,1]; 
        CA_all = [1,2; 1,3; 1,4; 1,5; 1,6; 1,7; 1,8; 1,9; 2,3; 2,4; 2,5; 2,6; 2,7; 2,8; 2,9; 3,4; 3,5; 3,6; 3,7; 3,8; 3,9; 4,5; 4,6; 4,7; 4,8; 4,9; 5,6; 5,7; 5,8; 5,9; 6,7; 6,8; 6,9; 7,8; 7,9; 8,9];
        %c_ratio = 1;
        %sidenum = 3;
        %nucFac = 1;
        
        designs_array = designs{:,:};
        f_penalized = csv_data{:,1:2};
        feas_scores = csv_data{:,3}; 
        stab_scores = csv_data{:,4};
        
        %axes(H.hParetoFig)
        pareto_bool = paretofront(f_penalized);
        f_pareto = f_penalized(pareto_bool==1,:);
        designs_pareto = designs_array(pareto_bool==1);
        feas_pareto = feas_scores(pareto_bool==1);
        stab_pareto = stab_scores(pareto_bool==1); 
        f_pareto_true = [f_pareto(:,1), -f_pareto(:,2)]; % objective 1 is to be
        % minimized while objective 2 is to be maximized  
        gui_object = source.Parent;
        plot(f_pareto_true(:,1),f_pareto_true(:,2),'*','ButtonDownFcn',{@plot_pareto_callback_gui,NC,CA_all,designs_pareto,feas_pareto,stab_pareto,f_pareto,sel,r,E,pen_fac,H,gui_object},'Parent',H.hParetoFig)
        set(get(H.hParetoFig, 'xlabel'), 'string', '$\frac{\left|\frac{C_{22}}{C_{11}} - c_{target}\right|}{15}$ + penalty','Interpreter','latex') 
        set(get(H.hParetoFig, 'ylabel'), 'string', '$\frac{\frac{C_{22}}{v_f}}{8500}$ + penalty','Interpreter','latex') 
        
        set(H.hEditFeasScore,'string',[]);
        set(H.hEditStabScore,'string',[]);
        set(H.hEditCMat,'Data',cell(3,3));
        set(H.hEditTrueObjs,'string',[]);
        cla(H.hTrussFig);
    end

    
end
     
    
        
            
      