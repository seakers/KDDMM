function [ParetoFront,ParetoSet,eMOGA]=evMOGA(eMOGA)
%
% Multiobjective Optimization Algorithm evMOGA
%
% [ParetoFront,ParetoSet,eMOGA]=evMOGA(eMOGA)
%    ParetoFront : Pareto front
%    ParetoSet   : Pareto set
%    eMOGA       : Data structure for algorithm parameters
%
%
% %% Parameters that have to (or can) be tuned by user
% %% Minimal parameters set: MOP description (have to be tuned by user)
% eMOGA.objfun           % m-function name for objectives computation
% eMOGA.objfun_dim       % Objective space dimension
% eMOGA.searchspaceUB    % Search space upper bound
% eMOGA.searchspaceLB    % Search space lower bound
%
% %% Additional MOP parameters (can be tuned by user)
% eMOGA.param            % Objectives additional parameter
%                        %   empty if not used
%                        %   default value: []
% %% Algorithm parameters (can be tuned by user)
% eMOGA.Nind_P           % Individuals for the P population
%                        %   defaultvalue:: 100
% eMOGA.Generations      % Number of generations have to be  >= 2
% eMOGA.n_div            % Number of division for each dimension
%                        % to build the objective space grid
%                        %   default value: 50 divisions per dim
% eMOGA.Nind_GA          % Individuals for aux population GA
%                        %   have to be an even number
%                        %   defaultvalue:: 8
% eMOGA.subpobIni        % Initial subpopulation
%                        %   empty if not used
%                        %   default value: []
% eMOGA.subpobIni_obj    % subpobIni objective values
%                        %   empty if not used
%                        %   default value: []
% eMOGA.dd_ini           % Parameter used for crossover
%                        %   defaultvalue:: 0.25
% eMOGA.dd_fin           % Parameter used for crossover
%                        %   defaultvalue:: 0.1
% eMOGA.Pm               % Mutation Probability
%                        %   defaultvalue:: 0.2
% eMOGA.Sigma_Pm_ini     % Gaussian Mutation Sigma_Pm
%                        %   defaultvalue:: 20
% eMOGA.Sigma_Pm_fin     % Gaussian Mutation Sigma_Pm
%                        %   defaultvalue:: 0.1
% %% Internal parameters
% eMOGA.searchSpace_dim  % Search space dimension
% eMOGA.epsilon          % box size
% eMOGA.max_f            % maximum objective functions values
% eMOGA.min_f            % minimum objective functions values
% eMOGA.ele_P:           % individuals of population P
% eMOGA.coste_P          % objective function values of P
% eMOGA.mod              % dominance parameter
% eMOGA.ele_A            % individuals of population A
% eMOGA.coste_A          % objective function values of A
% eMOGA.box_A            % box for individuals of A
% eMOGA.Nind_A           % size of pareto front
% eMOGA.ele_GA           % individuals of population GA
% eMOGA.coste_GA         % objective function values of GA
% eMOGA.box_GA           % box for individuals of GA
% eMOGA.gen_counter      % generation counter
%
%
% (c) 2005-2014 CPOH - Universitat Politecnica de Valencia
%    cpoh.upv.es
%



%global eMOGA
disp('------ evMOGA Algorithm -------')
% Inicializaciones ---------------------------------------------
disp('######### evMOGA initialization')
% if nargin==0
%     disp('######### using InievMOEA.m')
%     eMOGA=InievMOEA;
% end

% If the parameter doesn't exist in the data structure it is created with the default value
if ~isfield(eMOGA,'searchSpace_dim')
    eMOGA.searchSpace_dim=size(eMOGA.searchspaceUB,2);         % Search space size
    disp(['### Value assigned: eMOGA.searchSpace_dim=', num2str(eMOGA.searchSpace_dim)])
end
if ~isfield(eMOGA,'dd_ini')
    eMOGA.dd_ini= 0.25;         % Parameter used for crossover
    disp(['### Default value assigned: eMOGA.dd_ini=', num2str(eMOGA.dd_ini)])
end
if ~isfield(eMOGA,'dd_fin')
    eMOGA.dd_fin= 0.1;          % Parameter used for crossover
    disp(['### Default value assigned: eMOGA.dd_fin=', num2str(eMOGA.dd_fin)])
end
if ~isfield(eMOGA,'Pm')
    eMOGA.Pm= 0.2;              % Mutation Probability
    disp(['### Default value assigned: eMOGA.Pm=', num2str(eMOGA.Pm)])
end
if ~isfield(eMOGA,'Sigma_Pm_ini')
    eMOGA.Sigma_Pm_ini= 20.0;   % Gaussian Mutation Sigma_Pm
    disp(['### Default value assigned: Sigma_Pm_ini=', num2str(eMOGA.Sigma_Pm_ini)])
end
if ~isfield(eMOGA,'Sigma_Pm_fin')
    eMOGA.Sigma_Pm_fin= 0.1;    % Gaussian Mutation Sigma_Pm
    disp(['### Default value assigned: Sigma_Pm_fin=', num2str(eMOGA.Sigma_Pm_fin)])
end
if ~isfield(eMOGA,'Nind_GA')
    eMOGA.Nind_GA= 8;    % Individuals for the GA population
    disp(['### Default value assigned: Nind_GA=', num2str(eMOGA.Nind_GA)])
else
    eMOGA.Nind_GA=round(eMOGA.Nind_GA/2)*2; % Have to be an even number
end
if ~isfield(eMOGA,'Nind_P')
    eMOGA.Nind_P= 100;    % Individuals for the P population
    disp(['### Default value assigned: Nind_P=', num2str(eMOGA.Nind_P)])
end
if ~isfield(eMOGA,'n_div')
    eMOGA.n_div= 50*ones(1,eMOGA.objfun_dim);    % Number of division for each dimension
    % to build the objective space grid
    disp(['### Default value assigned: n_div=[', num2str(eMOGA.n_div),']'])
elseif size(eMOGA.n_div,2)==1
    eMOGA.n_div= eMOGA.n_div*ones(1,eMOGA.objfun_dim);
    disp(['### Value assigned for every objective: n_div=[', num2str(eMOGA.n_div),']'])
end
if ~isfield(eMOGA,'Generations')
    eMOGA.Generations= 100;    % Number of generations
    disp(['### Default value assigned: Generations=', num2str(eMOGA.Generations)])
end
if ~isfield(eMOGA,'param')
    eMOGA.param = [];     % Additional parameter of objective function
    %   empty if not used
    disp('### Additional objective function parameters empty')
end
if ~isfield(eMOGA,'subpobIni')
    eMOGA.subpobIni = [];    % (optional) Initial subpopulation
    %   empty if not required
end
if isfield(eMOGA,'subpobIni_obj')
    if (not(size(eMOGA.subpobIni,1)==size(eMOGA.subpobIni_obj,1)))
        disp('### subpobIni_obj dimensions is not consistent with subpobIni then it is set empty')
        eMOGA.subpobIni_obj = [];    % (optional) Initial subpopulation
        %   empty if not required
    end
    if (not(size(eMOGA.subpobIni_obj,2)==eMOGA.objfun_dim))
        disp('### subpobIni_obj dimensions is not consistent with objfun_dim then it is set empty')
        eMOGA.subpobIni_obj = [];    % (optional) Initial subpopulation
        %   empty if not required
    end
else
    eMOGA.subpobIni_obj = [];    % (optional) Initial subpopulation
    %   empty if not required
end


eMOGA.epsilon=ones(1,eMOGA.objfun_dim);
eMOGA.max_f=ones(1,eMOGA.objfun_dim);
eMOGA.min_f=ones(1,eMOGA.objfun_dim);

% Polupation P
eMOGA.ele_P=NaN(eMOGA.Nind_P,eMOGA.searchSpace_dim);
eMOGA.coste_P=NaN(eMOGA.Nind_P,eMOGA.objfun_dim);
eMOGA.mod=zeros(eMOGA.Nind_P,1); 

% Polupation A
eMOGA.ele_A=NaN(1,eMOGA.searchSpace_dim);
eMOGA.coste_A=NaN(1,eMOGA.objfun_dim);
eMOGA.box_A=NaN(1,eMOGA.objfun_dim);
eMOGA.Nind_A=0; % 0 means empty set

% Poblacion GA
eMOGA.ele_GA=NaN(eMOGA.Nind_GA,eMOGA.searchSpace_dim);
eMOGA.coste_GA=NaN(eMOGA.Nind_GA,eMOGA.objfun_dim);
eMOGA.box_GA=NaN(eMOGA.Nind_GA,eMOGA.objfun_dim);
eMOGA.gen_counter=0;

% auxiliar
if eMOGA.Generations<1 % Number of generations should be >= 1
    eMOGA.Generations=1;
end 
eMOGA.dd_fin_aux=((eMOGA.dd_ini/eMOGA.dd_fin)^2-1)/(eMOGA.Generations-1); 
eMOGA.Sigma_Pm_fin_aux=((eMOGA.Sigma_Pm_ini/eMOGA.Sigma_Pm_fin)^2-1)/(eMOGA.Generations-1);
% End initialization ----------------------------------------

% Create initial population 
disp('######### Generating initial population')
eMOGA=CrtpReal(eMOGA);

if isempty(eMOGA.subpobIni_obj)
    i=1;
else
    i=size(eMOGA.subpobIni_obj,1)+1;
    eMOGA.coste_P(1:i-1,:)=eMOGA.subpobIni_obj;
    disp(['######### Adding values of the objective functions for the subpopulation of ' num2str(i-1) ' individuals' ])
end
% Objective function evaluation
if isempty(eMOGA.param)
    for i=i:eMOGA.Nind_P
        eMOGA.coste_P(i,:)=feval(eMOGA.objfun,eMOGA.ele_P(i,:));
    end
else
    for i=i:eMOGA.Nind_P
        eMOGA.coste_P(i,:)=feval(eMOGA.objfun,eMOGA.ele_P(i,:),eMOGA.param);
    end
end
% GENERATION 0
eMOGA=Iteracion_MOEA(eMOGA,0);
% GENERATIONS 1 to end
for i=1:eMOGA.Generations
    % Create GA with individuals of  P and A
    % Crossover and mutation (Pm mutation probability)
    eMOGA=rellena_GA2(eMOGA);
    % Objective function evaluation of GA
    if isempty(eMOGA.param)
        for jj=1:eMOGA.Nind_GA
            eMOGA.coste_GA(jj,:)=feval(eMOGA.objfun,eMOGA.ele_GA(jj,:));
        end 
    else
        for jj=1:eMOGA.Nind_GA
            eMOGA.coste_GA(jj,:)=feval(eMOGA.objfun,eMOGA.ele_GA(jj,:),eMOGA.param);
        end
    end
    eMOGA=Iteracion_MOEA(eMOGA,i);
    % user action in each iteration
    eMOGA=evMOGAiteration(eMOGA); 
end
% Pareto Optimal set and front
ParetoSet=eMOGA.ele_A(1:eMOGA.Nind_A,:);
ParetoFront=eMOGA.coste_A(1:eMOGA.Nind_A,:);
eMOGA=evMOGAresults(eMOGA);

function eMOGA=CrtpReal(eMOGA)
% eMOGA=CrtpReal(eMOGA)
% Create random population for P  (variable eMOGA.ele_P)
[nind,nsearchSpace_dim]=size(eMOGA.subpobIni);
k=1;
if nind>0
    if nsearchSpace_dim~=eMOGA.searchSpace_dim
        disp(' #### Ignoring initial subpopulation, incorrect searchSpace_dim')
    else
        disp(['######### Adding subpopulation of ' num2str(nind) ' individuals' ])
        eMOGA.ele_P(1:nind,:)=eMOGA.subpobIni;
        k=nind+1;
    end
else
    disp('######### Empty initial subpopulation')
end

for j=1:eMOGA.searchSpace_dim
    for i=k:eMOGA.Nind_P
        eMOGA.ele_P(i,j)=(eMOGA.searchspaceUB(j)-eMOGA.searchspaceLB(j))*rand+eMOGA.searchspaceLB(j);
    end
end

function eMOGA=Iteracion_MOEA(eMOGA,gen_counter)
% eMOGA=Iteracion_MOEA(eMOGA,gen_counter)
%  execute one iteration/generation of the algorithm
%  gen_counter==0
%       Update of A from P
%  gen_counter>0
%       Update A from  GA

eMOGA.gen_counter=gen_counter;
if (gen_counter==0)
    % Objective function values of P have been calculated
    % Update of A from P (box belonging is calculated)
    % Individuals of P are marked according to dominance
    % mod used as support to make the non-dominated
    % mod(i)=0 not yet inspected
    % mod(i)=1 the individual is dominated
    % mod(i)=2 the individual is non-dominanted
    for i=1:(eMOGA.Nind_P-1)
        if eMOGA.mod(i)==1
            continue
        end % i is a dominated individual
        for j=(i+1):eMOGA.Nind_P
            if eMOGA.mod(j)==1
                continue
            end % j is a dominated individual
            a=domina(eMOGA.coste_P(i,:),eMOGA.coste_P(j,:));
            if (a==0)
                % Nothing to do (very frequent)
            elseif (a==3 || a==1)
                %i dominates j or they are the same, mark j as dominated and 
                % I follow with next j
                eMOGA.mod(j)=1;
                continue
            elseif (a==2)
                %j dominates i, mark I as dominated and  I follow with next i
                eMOGA.mod(i)=1;
                break;
            end
        end
        if (eMOGA.mod(i)==0) % if nobody dominates him 
            eMOGA.mod(i)=2;  % it is marked as non-dominated
        end
    end
    if (eMOGA.mod(eMOGA.Nind_P)==0)% If last individual is not marked
        eMOGA.mod(eMOGA.Nind_P)=2; % it is marked as non-dominated
    end                            
    % All non-dominated individuals have been marked
    % Next steps are: obtaining bounds, grid and inserting non-dominated if required
    % Initial values for max and min values: eMOGA.max_f=-inf  and  eMOGA.min_f=inf
    for j=1:eMOGA.objfun_dim
        eMOGA.max_f(j)=-inf;
        eMOGA.min_f(j)=inf;
    end
    % updating bounds
    for i=1:eMOGA.Nind_P 
        if (eMOGA.mod(i)==2) % if non-dominated, bounds are updated
            for j=1:eMOGA.objfun_dim
                if (eMOGA.coste_P(i,j)>eMOGA.max_f(j))
                    eMOGA.max_f(j)=eMOGA.coste_P(i,j);
                end
                if (eMOGA.coste_P(i,j)<eMOGA.min_f(j))
                    eMOGA.min_f(j)=eMOGA.coste_P(i,j);
                end
            end
        end
    end
    
    % Obtaining grid.
    for j=1:eMOGA.objfun_dim
        eMOGA.epsilon(j)=(eMOGA.max_f(j)-eMOGA.min_f(j))/eMOGA.n_div(j);
    end
    % Obtaining Box and inserting elements in A
    for i=1:eMOGA.Nind_P
        if(eMOGA.mod(i)==2) % if non-dominante, checking inclusion in A
            % box
            box=calcula_box(eMOGA,eMOGA.coste_P(i,:));
            % inserting in A
            eMOGA=Archivar(eMOGA,eMOGA.ele_P(i,:),eMOGA.coste_P(i,:),box);
        end
    end
    
else
    % GENERATIONS > 0
    % Obtaining box of individuals of GA and checking bounds (changing them
    % if necessary)
    % updating A with individuals non e-dominated of GA
    for i=1:eMOGA.Nind_GA  % FOR generations > 0
        % Checking if it is inside e-dominance area or not
        a=0;
        b=0;
        for j=1:eMOGA.objfun_dim
            if (eMOGA.coste_GA(i,j)>=eMOGA.max_f(j))
                a=a+1;
            end
            if (eMOGA.coste_GA(i,j)<=eMOGA.min_f(j))
                b=b+1;
            end
        end
        if (a==0) && (b==0)
            %Case Z1/CASE A,
            %Bounds of A are unchanged and checking if it is e-non-dominated
            eMOGA.box_GA(i,:)=calcula_box(eMOGA,eMOGA.coste_GA(i,:));
            eMOGA=Archivar(eMOGA,eMOGA.ele_GA(i,:),eMOGA.coste_GA(i,:),eMOGA.box_GA(i,:));
        elseif (a==eMOGA.objfun_dim) && (b==0)
            %Case Z3/CASe B,
            % e-dominated individual (rejected)
        elseif (a==0) && (b==eMOGA.objfun_dim)
            %Case Z2/CASE C,
            %All individuals of A disappear and only e-non-dominated individual is maintained
            eMOGA.Nind_A=1;
            % Inserted in first position
            eMOGA.ele_A(1,:)=eMOGA.ele_GA(i,:);
            eMOGA.coste_A(1,:)=eMOGA.coste_GA(i,:);
            % Updating bounds
            eMOGA.max_f=eMOGA.coste_GA(i,:);
            eMOGA.min_f=eMOGA.coste_GA(i,:);
        else
            %Case Z4/CASE D,
            %Checking dominance and remove from A individuals dominated
            %Finally new bounds and epsilon are recalculated
            %and after individuals are reinserted to check e-dominance
            %First checking if somebody dominates me 
            aux2=0;
            for j=1:eMOGA.Nind_A 
                c=domina(eMOGA.coste_A(j,:),eMOGA.coste_GA(i,:));
                if(c==1)
                    aux2=1;
                    break;
                end
            end
            if (aux2>0)
                
                %ending if somebody dominates me (equivalent to case Z3
                continue;
            end
            %As nobody dominates me I check the ones I dominate to remove it
            j=1;
            while j<=(eMOGA.Nind_A-1)
                c=domina(eMOGA.coste_GA(i,:),eMOGA.coste_A(j,:));
                if(c==1)
                    % removing from A
                    eMOGA.ele_A(j:eMOGA.Nind_A-1,:)=eMOGA.ele_A(j+1:eMOGA.Nind_A,:);
                    eMOGA.box_A(j:eMOGA.Nind_A-1,:)=eMOGA.box_A(j+1:eMOGA.Nind_A,:);
                    eMOGA.coste_A(j:eMOGA.Nind_A-1,:)=eMOGA.coste_A(j+1:eMOGA.Nind_A,:);
                    eMOGA.Nind_A=eMOGA.Nind_A-1;
                    j=j-1;
                end
                j=j+1;
            end
            % Checking the last one
            if (domina(eMOGA.coste_GA(i,:),eMOGA.coste_A(eMOGA.Nind_A,:)))
                eMOGA.Nind_A=eMOGA.Nind_A-1;
            end
            % Updating bounds
            eMOGA.max_f=eMOGA.coste_GA(i,:);
            eMOGA.min_f=eMOGA.coste_GA(i,:);
            for j=1:eMOGA.Nind_A
                for k=1:eMOGA.objfun_dim
                    if (eMOGA.coste_A(j,k)>eMOGA.max_f(k))
                        eMOGA.max_f(k)=eMOGA.coste_A(j,k);
                    end
                    if (eMOGA.coste_A(j,k)<eMOGA.min_f(k))
                        eMOGA.min_f(k)=eMOGA.coste_A(j,k);
                    end
                end
            end
            %Obtaining epsilon
            for j=1:eMOGA.objfun_dim
                eMOGA.epsilon(j)=(eMOGA.max_f(j)-eMOGA.min_f(j))/eMOGA.n_div(j);
            end
            % Recalculating its box and marking empty_A=2 and pass them to
            % function 'archivar'
            % I store the number of individuals to analyse  to avoid using  empty_A=2
            
            for j=1:eMOGA.Nind_A
                eMOGA.box_A(j,:)=calcula_box(eMOGA,eMOGA.coste_A(j,:));
            end
            
            Nind_A_temp=eMOGA.Nind_A;
            % Not necessary to analyse the first individual
            % Analysing the rest 
            copiaele=eMOGA.ele_A(1:eMOGA.Nind_A,:);
            copiacos=eMOGA.coste_A(1:eMOGA.Nind_A,:);
            copiabox=eMOGA.box_A(1:eMOGA.Nind_A,:);
            eMOGA.Nind_A=1;
            for j=2:Nind_A_temp
                eMOGA=Archivar(eMOGA,copiaele(j,:),copiacos(j,:),copiabox(j,:));
            end
            % Storing in A the individual checked
            eMOGA.box_GA(i,:)=calcula_box(eMOGA,eMOGA.coste_GA(i,:));
            eMOGA=Archivar(eMOGA,eMOGA.ele_GA(i,:),eMOGA.coste_GA(i,:),eMOGA.box_GA(i,:));
        end %END  else cases Z1 Z2 Z3 and Z4
        eMOGA=Actualizar_P(eMOGA,eMOGA.ele_GA(i,:),eMOGA.coste_GA(i,:));
    end % FOR generation>0
end % ELSE generation>0

function eMOGA=Archivar(eMOGA,x,coste_f,box_f)
% eMOGA=Archivar(eMOGA,x,coste_f,box_f)
% Try to insert individual in A

j=1;
while (j<=eMOGA.Nind_A) 
    % Checking box-dominance
    a=box_domina(box_f,eMOGA.box_A(j,:));
    if (a==0) %none dominate the other
        %Nothing to do
    elseif (a==2)
        % Can't enter in A because it is box-dominated
        % I have to found it
        return;
    elseif (a==3)
        % They are in the same box (but one of them have to be removed)
        b=domina(coste_f,eMOGA.coste_A(j,:));
        switch b
            case 0
                % The one nearest of the coordinates of the box is
                % maintained
                dist1=norm((coste_f-((box_f-0.5).*eMOGA.epsilon+eMOGA.min_f))./eMOGA.epsilon,2);
                dist2=norm((eMOGA.coste_A(j,:)-((box_f-0.5).*eMOGA.epsilon+eMOGA.min_f))./eMOGA.epsilon,2);
               if (dist1<dist2)
                    eMOGA.ele_A(j,:)=x;
                    eMOGA.coste_A(j,:)=coste_f;
                end
                return;
            case 1
                % I substitute g by f
                eMOGA.ele_A(j,:)=x;
                eMOGA.coste_A(j,:)=coste_f;
                return;
            case {2,3}
                %If they are identical or the one in A dominates the new one, 
                %the one in A is maintained
                return;
            otherwise
                error('error switch b')
        end
    else
        % box-dominate a individual of A then this individual is removed
        if j<eMOGA.Nind_A %If it is in the last position there is nothing to move, only modifying Nind_A
            eMOGA.ele_A(j:eMOGA.Nind_A-1,:)=eMOGA.ele_A(j+1:eMOGA.Nind_A,:);
            eMOGA.box_A(j:eMOGA.Nind_A-1,:)=eMOGA.box_A(j+1:eMOGA.Nind_A,:);
            eMOGA.coste_A(j:eMOGA.Nind_A-1,:)=eMOGA.coste_A(j+1:eMOGA.Nind_A,:);
        end
        j=j-1;
        eMOGA.Nind_A=eMOGA.Nind_A-1;%Removing
    end
    j=j+1;
end
%Adding at the end
eMOGA.ele_A(eMOGA.Nind_A+1,:)=x;
eMOGA.coste_A(eMOGA.Nind_A+1,:)=coste_f;
eMOGA.box_A(eMOGA.Nind_A+1,:)=box_f;
eMOGA.Nind_A=eMOGA.Nind_A+1;

function eMOGA=rellena_GA2(eMOGA)
% eMOGA=rellena_GA2(eMOGA)
% Creates GA from P and A (randomly)
% crossover and mutation (Pm probability of mutation)

for i=1:2:eMOGA.Nind_GA
    %Selecting two individuals
    %from P
    a=uniform_int(eMOGA.Nind_P) ; 
    eMOGA.ele_GA(i,:)=eMOGA.ele_P(a,:);
    %from A
    b=uniform_int(eMOGA.Nind_A);
    eMOGA.ele_GA(i+1,:)=eMOGA.ele_A(b,:);
    
    pm=rand;
    if (pm>eMOGA.Pm)%crossover
        eMOGA=Lxov(eMOGA,i,i+1);
    else  %mutation with Pm probability
        eMOGA=MutRealDGausPond2(eMOGA,i);
        eMOGA=MutRealDGausPond2(eMOGA,i+1);
    end
end

function x=uniform_int(xmax)
% x=uniform_int(xmax)
% random integer between 1 and xmax (uniform distribution)
x=ceil(rand*xmax);


function eMOGA=MutRealDGausPond2(eMOGA,p)
% eMOGA=MutRealDGausPond2(eMOGA,p)
% Perform mutation of an individual
% shr in percentage of the range of each variable
% Mutants are maintained inside the search space bound

padre=eMOGA.ele_GA(p,:);
sigma=eMOGA.Sigma_Pm_ini/sqrt(1+eMOGA.gen_counter*eMOGA.Sigma_Pm_fin_aux);

for j=1:eMOGA.searchSpace_dim
    padre(j)=padre(j)+gauss(1,1,0,(eMOGA.searchspaceUB(j)-eMOGA.searchspaceLB(j))*sigma/100.0);
    if(padre(j)>eMOGA.searchspaceUB(j))
        padre(j)=eMOGA.searchspaceUB(j);
    elseif (padre(j)<eMOGA.searchspaceLB(j))
        padre(j)=eMOGA.searchspaceLB(j);
    end
end
eMOGA.ele_GA(p,:)=padre;

function aux=gauss(i,j,m,s2)
% aux=gauss(i,j,m,s2)
% Generates a ixj matrix of random numbers
% normal distribution (mean 'm' and standard deviation 's2')
aux = m + s2* randn(i,j);

function box=calcula_box(eMOGA,coste)
% box=calcula_box(eMOGA,coste)
% Obtain box from cost vector and epsilon vector

box=zeros(1,eMOGA.objfun_dim);
for i=1:eMOGA.objfun_dim
    if(eMOGA.epsilon(i)==0.0)
        box(i)=0;
    else
        box(i)=ceil(round(1e6*(coste(i)-eMOGA.min_f(i))/eMOGA.epsilon(i))/1e6);
    end
end

function eMOGA=Actualizar_P(eMOGA,x,coste_f)
% eMOGA=Actualizar_P(eMOGA,x,coste_f)
% Analyse if an individual have to be inserted in P

d=0;
% Random selection a position in P to begin the search of a dominated individual
a=uniform_int(eMOGA.Nind_P);
for i=a:eMOGA.Nind_P
    if (domina(coste_f,eMOGA.coste_P(i,:))==1)%f dominates i
        d=i;
        break; %found
    end
end

if d==0 %Searching in the rest of the P
    for i=1:a-1
        if (domina(coste_f,eMOGA.coste_P(i,:))==1) %f dominates i
            d=i;
            break; %found
        end
    end
end

if d~=0
    %Substituting  dominated individual
    eMOGA.ele_P(i,:)=x;
    eMOGA.coste_P(i,:)=coste_f;
end

function aux=domina(f,g)
% aux=domina(f,g)
% Check dominance:
% 0 if neither dominates the other
% 1 if f dominates g,
% 2 if g dominates f
% 3 they are the same

a=0;
b=0;
for i=1:size(f,2)
    if (f(i)<g(i))
        a=1; %g non-dominates f
    elseif (f(i)>g(i))
        b=1; %f non-dominates g
    end 
end
aux=3-a*2-b;

function aux=box_domina(box_f,box_g)
% aux=box_domina(box_f,box_g)
% Check box-dominance
% 0 if neither e-dominates the other
% 1 if f_box dominates g_box,
% 2 if g_box dominates f_box
% 3 same box

a=0;
b=0;
for i=1:size(box_f,2)
    if (box_f(i)<box_g(i))
        a=1;
        %box_g non-dominate box_f
    elseif (box_f(i)>box_g(i))
        b=1;
    end %box_f non-dominate box_g
end
aux=3-a*2-b;

function eMOGA=Lxov(eMOGA,p1,p2)
% eMOGA=Lxov(eMOGA,p1,p2)
% Crossover of two individuals
% offsprings are maintained inside the search space bound

padre1=eMOGA.ele_GA(p1,:);
padre2=eMOGA.ele_GA(p2,:);

dd=eMOGA.dd_ini/sqrt(1+eMOGA.gen_counter*eMOGA.dd_fin_aux);

beta=(1.0+2.0*dd)*rand-dd;
hijo1=beta*padre1+(1.0-beta)*padre2;
hijo2=(1.0-beta)*padre1+beta*padre2;
for k=1:eMOGA.searchSpace_dim
    if(hijo1(k)>eMOGA.searchspaceUB(k))
        hijo1(k)=eMOGA.searchspaceUB(k);
    elseif (hijo1(k)<eMOGA.searchspaceLB(k))
        hijo1(k)=eMOGA.searchspaceLB(k);
    end
    if(hijo2(k)>eMOGA.searchspaceUB(k))
        hijo2(k)=eMOGA.searchspaceUB(k);
    elseif (hijo2(k)<eMOGA.searchspaceLB(k))
        hijo2(k)=eMOGA.searchspaceLB(k);
    end
end
eMOGA.ele_GA(p1,:)=hijo1;
eMOGA.ele_GA(p2,:)=hijo2;
