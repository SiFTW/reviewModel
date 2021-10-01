% NF-kB Model Version 6.2 RelA + RelB dimers, Run File
% Updated with following features:
% 1) p100 txn
% 2) add in multipliers for NIK and IKK -- IKKmultipli, NIK_deg_hmod
% 3) add in NEMOpulse for simulation of priming effect
% Last modified: 2010/09/09 by Allison.  Add in one NIK_TYPE:'basal_multi'.
% Add in NIK_multiplier.

function v=nfkb_run_62_passParams1_mod(loop_cell,paramsVersion)

global time; 
if nargin == 0
    tic;
end
v = struct;

%==== Simulation Parameters
v.START_TIME      = -1e6;   % seconds of basal state
v.SIM_TIME        = 100;    % seconds of stimulation phase
v.PLOT_FLAG       = false;   % Plot results?
v.PLOT_POS        = false;  % Used in plot function. default = false
v.graphs          = 0;      % Used to plot stuff, if chosen.
v.DISP            = true;   % Display simulation text?
v.SEC_RUN         = false;  % Run in seconds(true) or minutes (false)
v.CONSTRAINT_FLAG = false;  % Check results against constraints?
v.CONSTRAINT_LOOP = false;  % True to return only constraint results
v.RATE_MODS       = {};     % set by input argument (below)

%==== IKK Input Curves
v.TOTAL_IKK2      = 140; % nM  ~200k IKK dimers (G Ghosh measured)
v.BASAL_IKK2      = .01; % 1% of total
ikk_curves        = {'TNFp15' 'IGM' 'LPS45p' 'basal' 'custom' 'custom_TNFp15' 'BAFF' 'custom_LPSc2'};
ikk_curve_num     = 7;  % index of ikk curve to use (1 ... 5)
v.IKK2_TYPE       = ikk_curves{ikk_curve_num};
v.IKK2_FRACTION   = 0.5;% What fraction of total is IKK2 curve going to reach? For custom_LPSc2;
v.NEMOpulse       = 0;  % pretreatment of NEMO?
v.Pulsetime       = 0;  % pulse length of NEMO pretreatment

v.NIK_INPUT       = 1; % logical : NIK activation or not
v.NIK_TYPE        ='basal';   %'basal', 'basal_10_5hr', 'basal_1_5hr', 'basal_1_20hr'
v.NIK_deg_mod1    = 5;
v.NIK_deg_mod2    = 5; 
v.NIKdelay        = 0;  % delay time of NIK
v.NIK_multiplier  = 1;

v.LOOP_NEMO       = 0; % NEMO input is generated from outside function
v.LOOP_NIK        = 0; % NIK input is generated from outside functon

v.NIK_mod         = []; % leave blank, don't change
v.IKK2            = []; % leave blank, don't change
v.IKK2multipli    = 1 ;

%==== Number of species and # components for each
v.IKBS            = 4; % Number of IkBs a, b, e, d
v.IKBtext = {'IkBa','IkBb','IkBe','IkBd'};
v.IKB_COMP        = 3; % Number of components per IkB (mRNA,cyt,nuc)
v.DIMERS          = 7; % Number of NFkB Dimers
v.DIMERStext = {'AA','A50','A52','B50','B52','5050','5252'};
v.DIMER_COMP      = 10;% Number of components per Dimer
v.MONOMERS        = 6; % Number of Monomers  A (1), 50 (2), B(3), 100(4) 52 (5), C (6)
v.MONOMERStext    = {'RelA','p50','RelB','p100','p52', 'cRel'};
v.MONOMER_COMP    = 3; % Number of components per Monomer
v.OTHER_COMP      = 1; % currently just NIK protein (free in cyt)
v.OTHERtext    = {'NIK'};

%==== Species Knockouts  1 = present / 0 = knocked-out
v.IKB_KO          = logical([1;1;1;1]);        % IkBa ; b ; e ; d
v.MONOMER_KO      = logical([1;1;1;1;1;1]);      % RelA ; p50 ; RelB ; p100 ; p52; cRel
v.DIMER_KO        = logical([1;1;1;1;1;1;1]);  % AA; A50; B50; B52; C50; C52; B100



set_values_names(); %Set Initial Values and Names
    function set_values_names()
        v.SPECIES     = cell(2,1); % initial concentrations and names
        v.numSpecies = v.IKB_COMP*v.IKBS + v.MONOMER_COMP*v.MONOMERS + v.DIMER_COMP*v.DIMERS + v.OTHER_COMP;
        v.SPECIES{1}  = zeros(v.numSpecies,1); % all species start at zero [] by default
        
        % Oh shoot, I'm gonna have to fix this ....
        
        v.SPECIES{2}  = { ...  % Component # 1-12 , 13-27 , 28-37, 38-47, 48-57, 58-67,68
            'tIkBa','IkBa', 'IkBan','tIkBb','IkBb','IkBbn','tIkBe','IkBe','IkBen','tIkBd','IkBd','IkBdn',...
            'tRelA', 'RelA','RelAn','tP50','p50','p50n', 'tRelB', 'RelB','RelBn','tp100','p100','p100n','tp52','p52','p52n'...
            'A50','A50n','IkBaA50','IkBaA50n','IkBbA50','IkBbA50n','IkBeA50','IkBeA50n','IkBdA50','IkBdA50n',...
            'B50','B50n','IkBaB50','IkBaB50n','IkBbB50','IkBbB50n','IkBeB50','IkBeB50n','IkBdB50','IkBdB50n',...
            'B52','B52n','IkBaB52','IkBaB52n','IkBbB52','IkBbB52n','IkBeB52','IkBeB52n','IkBdB52','IkBdB52n',...
            'B100','B100n','IkBaB100','IkBaB100n','IkBbB100','IkBbB100n','IkBeB100','IkBeB100n','IkBdB100','IkBdB100n',...
            'NIK' ...
            }';
%         v.masterList = masterList;
    end


%====  Get the parameters.
eval(sprintf('v=%s(v);',paramsVersion));

%==== Deal with variables passed into the run file
%{
        Input variable is a cell array
        Each row contains a v struct field name and a replacement value
            ie: {'SIM_TIME', 3600} or {'v.IKK2_TYPE' , 'TNF15p'}
        You cannot use the input to create a new struct field
        Likewise, fields set after this section of code cannot be changed
            (CURRENT_IKK, parameter lists, species names, etc)
        A special input field 'RATE_MODS' is used in the parameter function
            {'RATE_MODS',{[param #s],mod value}}
%}


if nargin > 0
    for j = 1:size(loop_cell,1)        
        if isfield(v,loop_cell{j,1})
            v.(loop_cell{j,1}) = loop_cell{j,2};
        else
            disp(j);
            error('nfkb_run_62: You cannot create a new struct field.');
        end
    end
end

%==== Run Setup Functions
if ~v.LOOP_NEMO  % unless NEMO was set as an input, generate curve
    generate_IKK2curve(v.IKK2multipli); % Generate IKK Curve
end

    function generate_IKK2curve(IKK2multipli)
        % IKK library is in minutes so it has to be converted to sec
        if v.SEC_RUN
            convert = 60; % 60 seconds per minute
            ikk_time = v.SIM_TIME;
        else
            convert = 1;
            ikk_time = v.SIM_TIME /60;
        end
        
        if (strcmp(v.IKK2_TYPE,'TNF15p'))
            values = [1 60 100 65 50 36 21 16 10 1 1 1] / 100;
            IKK2multiarray=ones(size(values));
            IKK2multiarray(2:3)=IKK2multipli;
            if IKK2multipli>1
                IKK2multiarray(4:12)=interp1([15 2400],[IKK2multipli 1], [15 20 25 30 45 59 60 90 2400]);
            elseif IKK2multipli<1
                IKK2multiarray(4:12)=IKK2multipli;
            end
            values=values.*IKK2multiarray;
            time = [0 5 10 15 20 25 30 45 59 60 90 2400]*convert;
%         elseif (strcmp(v.IKK2_TYPE,'IGM'))
%             values  = [1 2 3 5 14 17 25 23 20  13 4  2   2   1 1] / 100;
%             time    = [0 5  10 15 20 25 30 60 90 120 240 360 480 600 1441]*convert;    
        elseif (strcmp(v.IKK2_TYPE,'IGM'))
            values  = [1 2 3 5 14 17 25 36 35  34   30  17 4*3  2*5   2*5   1*10 1*10 1*10 1*10] / 100;
            values = [1 2 3 5 14 25 30 35 36  35   30  17 4*3  2*5   2*5   1*10 1*10 1*10 1*10]/100;
            %values  = [1 5 15 25 34 45 80 100 60  34   30  17 4*3  2*5   2*5   1*10 1*10 1*10 1*10] / 100;
            values = [1 5 15 25 45 65 75 100 90  75  50  45 20  15   2*5   1*10 1*10 1*10 1*10]/100*IKK2multipli;
            %values  = [1 2 3 5 14 17 25 36 35  34   30  30 30 30 28 27 25 24 20] / 100;
            %time = time*convert;
            time    = [0 5 10 15 20 25 30 60 90 120 180 210 240 360 480 600 1440 2160 2880]*convert;   
%         elseif (strcmp(v.IKK2_TYPE,'IGM'))
%            values  = [1 5 8 15 25 60 90 60 30  25 20 20  15 15 15] / 100;
%             time    = [0 5 10 15 20 25 30 60 90 120 240 360 480 600 1441]*convert; 
        elseif (strcmp(v.IKK2_TYPE,'IGM_mod'))
            values = [0 85 100 85 45 30 25 23 18 15 12 10 10]/100*IKK2multipli;
            %values = [0 85 100 85 75 60 55 43 28 15 12 10 10]/100*IKK2multipli;
            %values = [0 100 85 10 5 4 3 2 1 1 0.5 0.5 0.5]/100*IKK2multipli;
            
            %correct curve
            time    = [0 5 10 15 18 20 300 400 450 480 500 520 2880]*convert;  
            
            %slightly extended
            %time    = [0 15 30 55 78 100 300 400 450 480 500 520 2880]*convert;  
            
            %more extended
            %time    = [0 35 60 95 148 250 300 400 450 480 500 520 2880]*convert;  
            
            %super late and long
            %time    = [0 55 90 125 188 350 380 450 480 490 500 520 2880]*convert;  
            
        elseif (strcmp(v.IKK2_TYPE,'BAFF'))
            values = [1 2 3 5 7 8 9 10  9  8 7 6  5 4 1 1 1 1 1] /100;
            values = [1 2*2 3*2 5*2 7*2 8*2 9*2 10*2  9*2  8*2 7*2 6*1.5  5 4 1 1 1 1 1] /100;
            %time = time*convert; 
           time =   [0 5  10 15 20 25 30 60 90 120 180 210 240 360 480 600 1440 2160 2880]*convert;
        elseif (strcmp(v.IKK2_TYPE,'costim'))
            values  = ([1 2 3 5 14 17 25 36 35  34   30  17 4*3  2*5   2*5   1*10 1*10 1*10 1*10]) / 100;
            values = [1 2 3 5 14 25 30 35 36  35   30  17 4*3  2*5   2*5   1*10 1*10 1*10 1*10]/100;
            %values = [1 5 15 25 34 45 80 100 60  34   30  17 4*3  2*5   2*5   1*10 1*10 1*10 1*10] / 100;
            values = [1 5 15 25 45 65 75 100 90  75  50  45 20  15   2*5   1*10 1*10 1*10 1*10]/100;

            %time = time*convert;
            time    = [0 5  10 15 20 25 30 60 90 120 180 210 240 360 480 600 1440 2160 2880]*convert; %changed from just using IGM curve
%         elseif (strcmp(v.IKK2_TYPE,'costim'))
%             values  = ([1 5 8 15 25 60 90 60 30  25 20 20  15 15 15for x = ]) / 100;
%             time    = [0 5  10 15 20 25 30 60 90 120 240 360 480 600 1441]*convert; %changed from just using IGM curve

        elseif (strcmp(v.IKK2_TYPE,'LPS45p'))
            values = [1 2 7 13 19 25 7 6 3 1 1] /100;
            IKK2multiarray=ones(size(values));
            IKK2multiarray(2:5)=IKK2multipli;
            if IKK2multipli>1
                IKK2multiarray(6:11)=interp1([90 2400],[IKK2multipli 1], [90 120 240 360 600 2400]);
            elseif IKK2multipli<1
                IKK2multiarray(1:11)=IKK2multipli;
            end
            values=values.*IKK2multiarray;
            time = [0 15 30 45 60 90 120 240 360 600 2400]*convert;
        elseif (strcmp(v.IKK2_TYPE,'basal'))
            values = [1 1] /100 * IKK2multipli;
            time = [0 2880]*convert;
        elseif (strcmp(v.IKK2_TYPE,'custom'))
            values = [1 1 2 2.5 2.5] / 100 * IKK2multipli;
            time = [0 60 480 960 2400]*convert;
        elseif (strcmp(v.IKK2_TYPE,'Basal_pulse'))
            values = [v.BASAL_IKK2*100 v.NEMOpulse v.NEMOpulse v.BASAL_IKK2*100 v.BASAL_IKK2*100 ]/100;
            time = [0 2 v.Pulsetime v.Pulsetime+2 v.SIM_TIME]*convert;
        elseif (strcmp(v.IKK2_TYPE,'custom_TNFp15'))
            values = [1 60 100 65 50 36 21 16 10 2.5 2.5 2.5] / 100;
            IKK2multiarray=ones(size(values));
            IKK2multiarray(2:3)=IKK2multipli;
            if IKK2multipli>1
                IKK2multiarray(4:12)=interp1([15 2400],[IKK2multipli 1], [15 20 25 30 45 59 60 90 2400]);
            elseif IKK2multipli<1
                IKK2multiarray(4:12)=IKK2multipli;
            end
            values=values.*IKK2multiarray;
            time = [0 5 10 15 20 25 30 45 59 60 90 2400]*convert;
            
        elseif (strcmp(v.IKK2_TYPE,'custom'))
            values = [1 2 3 3 3 5 3 3 2.5 2.5 2.5] /100;
            IKK2multiarray=ones(size(values));
            IKK2multiarray(2:5)=IKK2multipli;
            if IKK2multipli>1
                IKK2multiarray(6:11)=interp1([90 2400],[IKK2multipli 1], [90 120 240 360 600 2400]);
            elseif IKK2multipli<1
                IKK2multiarray(1:11)=IKK2multipli;
            end
            values=values.*IKK2multiarray;
            time = [0 15 30 45 60 90 120 240 360 600 2400]*convert;
        elseif (strcmp(v.IKK2_TYPE,'custom_LPSc2'))
            values = [1 [2 7 13 19 25:-1:20]*v.IKK2_FRACTION] /100;
            IKK2multiarray=ones(size(values));
            IKK2multiarray(2:5)=IKK2multipli;
            if IKK2multipli>1
                IKK2multiarray(6:11)=interp1([90 2400],[IKK2multipli 1], [90 120 240 360 600 2400]);
            elseif IKK2multipli<1
                IKK2multiarray(1:11)=IKK2multipli;
            end
            values=values.*IKK2multiarray;
            time = [0 15 30 45 60 90 120 240 360 600 2400]*convert;
            
            
        else
            error('nfkb_run_62:generate_ikk_curve() Incorrect IKK input name');
        end
        screen(['IKK ' v.IKK2_TYPE],v.DISP);
        
 
        v.IKK2 = interp1(time, values, 0:ikk_time,'pchip');
    end

%====
if ~v.LOOP_NIK  % unless NIK was set as an input, generate curve
    generate_NIKcurve(); % degradation modifier curve
end

    function generate_NIKcurve()
        % NIK library is in minutes so it has to be converted to sec
        if v.SEC_RUN
            convert = 60; % 60 seconds per minute
            nik_time = v.SIM_TIME;
        else
            convert = 1;
            nik_time = v.SIM_TIME /60;
        end
        
        if v.NIK_INPUT
            mod = v.NIK_deg_mod1; % modified at 5hr point
            mod2 = v.NIK_deg_mod2;
        else
            mod = 1;  % no modifier for NIK degradation (No NIK input)
            mod2 = 1;
        end
        
        if (strcmp(v.NIK_TYPE,'basal_10_5hr'))
            deg_mod  = [.1 .1 .1/mod2 .1/mod2 .1/mod2]; % generate basal NIK=10
            deg_time = [0 60 300 720 2880]*convert;
        elseif (strcmp(v.NIK_TYPE,'basal_0hr'))
            deg_mod= [1/mod 1/mod]./v.NIK_multiplier;
            deg_time=[0 2880]*convert;
        elseif (strcmp(v.NIK_TYPE,'basal_multi'))  %% Added on 09/09/2010
            deg_mod= [1 1]./v.NIK_multiplier;
            deg_time=[0 1440]*convert;
         elseif (strcmp(v.NIK_TYPE,'basal_multi_0hr'))
            deg_mod= [1/mod 1/mod]./v.NIK_multiplier;
            deg_time=[0 2880];%deg_time=[0 1440]*convert;
        elseif (strcmp(v.NIK_TYPE,'basal'))
            deg_mod= [1 1];
            deg_time=[0 1440]*convert;
        elseif (strcmp(v.NIK_TYPE,'basal_5hr'))
            deg_mod  = [1 1 1/mod 1/mod 1/mod 1/mod 1/mod]; % decrease deg 12 fold by 5hr point
            deg_time = [0 60 300 720 1440 2160 2880]*convert;
        elseif (strcmp(v.NIK_TYPE,'basal_5hr_2'))
            deg_mod  = [1 1 1/mod 1/mod2 1/mod2]; % decrease deg 12 fold by 5hr point
            deg_time = [0 60 420 720 1440]*convert;
        elseif (strcmp(v.NIK_TYPE,'basal_delay'))
            deg_mod  = [1 1 1/mod 1/mod2 1/mod2]; % start at NIKstarttime
            deg_time = ([0 60 300 720 2880]+v.NIKdelay)*convert;
        end
        
        v.NIK_mod = interp1(deg_time,deg_mod,0:nik_time,'pchip');
    end


%==== Adjust for Rate modifiers set in v.RATE_MODS
adjust_parameters();
    function adjust_parameters()
        %  v.RATE_MODS(:,1) = param list identifier (pi,pm,pd,pdi,pdm)
        %  v.RATE_MODS(:,2) = vector of param numbers
        %  v.RATE_MODS(:,3) = modifier value
        for i = 1:size(v.RATE_MODS,1)
            switch lower(v.RATE_MODS{i,1})
                case {'pi'}
                    v.pi(v.RATE_MODS{i,2}) = ...
                        v.pi(v.RATE_MODS{i,2}) * v.RATE_MODS{i,3};
                case {'pm'}
                    v.pm(v.RATE_MODS{i,2}) = ...
                        v.pm(v.RATE_MODS{i,2}) * v.RATE_MODS{i,3};
                case {'pd'}
                    v.pd(v.RATE_MODS{i,2}) = ...
                        v.pd(v.RATE_MODS{i,2}) * v.RATE_MODS{i,3};
                case {'pdi'}
                    v.pdi(v.RATE_MODS{i,2}) = ...
                        v.pdi(v.RATE_MODS{i,2}) * v.RATE_MODS{i,3};
                case {'pdm'}
                    v.pdm(v.RATE_MODS{i,2}) = ...
                        v.pdm(v.RATE_MODS{i,2}) * v.RATE_MODS{i,3};
                case {'ps'}
                    v.ps(v.RATE_MODS{i,2}) = ...
                        v.ps(v.RATE_MODS{i,2}) * v.RATE_MODS{i,3};
                otherwise
                    error('nfkb_run_62: improper RATE_MOD identifier');
            end
        end
        
        %====  Adjust knockout arrays and synthesis parameters
        % Set the IkB or monomer synthesis rate to zero if knocked out
        v.pi(1,:) = v.pi(1,:) .* v.IKB_KO';     %ko ikb txn
        v.pm(1,:) = v.pm(1,:) .* v.MONOMER_KO'; %ko monomer txn
        
        % If a monomer is knocked out, also knock out its dimers
        for i = 1:v.DIMERS
            if ~v.MONOMER_KO(v.pd(1,i),1) || ~v.MONOMER_KO(v.pd(2,i),1)
                v.DIMER_KO(i) = false;
            end
        end
    end

%====


%==== Run the Simulation
run_sim();
    function run_sim()
        v.OUTPUT1 = []; % phase 1 (equilibrium) numerical results
        v.OUTPUT2 = []; % phase 2 (stimulated) numerical results
        % While the simulation is in seconds, the results are per minute
        if ~v.SEC_RUN  %you want to run in minutes
            v.START_TIME = floor(v.START_TIME/60);
            v.SIM_TIME   = floor(v.SIM_TIME/60);
        end
        
        global FLUX_CELL;
        global FLUX_VEC;
        global FLUX_T;
        global timeStepCount;
        FLUX_CELL = {};
        FLUX_VEC = [];
        FLUX_T = [];
        timeStepCount = 0;
        
        %--- Call the ODE function to reset persistent variables
        nfkb_ode_62_mod([],[],[],v);
        
        %--- Equilibrium Phase
        v.PHASE = 1;
        %screen(['Phase ' num2str(v.PHASE)],v.DISP);
        
        static      = false;    % True if at equilibrium
        count       = 1;        % Iteration Counter
        threshold   = .01;      % Max Fraction Change for equil. state
        absthreshold=0.000001;
        initvalues = v.SPECIES{1};
        r1          = [];
        
        while ~static
            [t1, r1] = ode15s('nfkb_ode_62_mod', [v.START_TIME 0], initvalues,[],v);
            fprintf('running first time\n');
            % Evaluate results and return true if at equilibrium
            %   else, increment the counter and update initvalues
            %diff = max(abs( 1 - (r1(1,:) ./ r1(end,:))));
            diff = max(abs((r1(end,:)-r1(end-1,:))./r1(end,:)));
            absdiff=max(abs((r1(end,:)-r1(end-1,:))));
            
            %if ~any((diff < threshold) && (diff < 1e-20),1)
            %diff
            if (diff < threshold||absdiff<absthreshold)
                v.START_TIME = (count * v.START_TIME) - (count-1);
                static = true;
                fprintf('static found\n');

            elseif count > 100
                error('=====> maximum phase 1 steps reached (100 steps)\n');
            else
                count = count + 1;
                initvalues = r1(end,:)';
                %fprintf('looping again in phase 1\n');
            end
        end   

        
%         setInitialValues62(r1(end,:)');
        v.OUTPUT1 = r1(end,:)';
        
        %--- Stimulation Phase
        if v.SIM_TIME > 0
            v.PHASE = 2;
            FLUX_CELL = {};
            FLUX_VEC = [];
            FLUX_T = [];
            timeStepCount = 0;
            
            [t2, r2] = ode15s('nfkb_ode_62_mod', [0 v.SIM_TIME], r1(end,:)',[],v);
            FLUX_T = FLUX_T(FLUX_T ~= 0);
            FLUX_VEC = FLUX_VEC(1:length(FLUX_T),:); % Get rid of padded zeros
            
            if v.SEC_RUN  % the model ran in seconds, but results are converted to min
                v.OUTPUT2   = interp1(t2, r2(:,:),(0:60:v.SIM_TIME));
                [b, m, n] = unique(FLUX_T);
                v.FLUXES = interp1(FLUX_T(m), FLUX_VEC(m,:),(0:60:v.SIM_TIME));
            else % the model ran in minutes
                v.OUTPUT2   = interp1(t2, r2(:,:),(0:1:v.SIM_TIME));
                [b, m, n] = unique(FLUX_T);
                v.FLUXES = interp1(FLUX_T(m), FLUX_VEC(m,:),(0:1:v.SIM_TIME));
            end
             % Because of the interpolation, the first time point may be NaN, so set it to the next point
            if isnan(v.FLUXES(1))
                v.FLUXES(1,:) = v.FLUXES(2,:);
            end
            
        else % only running phase 1 to get basal levels
            v.OUTPUT2 = real(r1(end,:));
            v.FLUXES = FLUX_VEC(end,:)';
        end                
        
    end

%==== Check criteria
if v.CONSTRAINT_FLAG
    v.constraints = run_constraints();
    
    if v.CONSTRAINT_LOOP
        v = v.constraints;
    end
end
    function constraints = run_constraints()
        
        % Logical Array for constraints
        constraints = ones(8,1);  % number of constraints
        
        % 1: basal NFkBn < 5% of total NFkB (p65 het+hom)
        Active_NFkB = [20 30];  % p65 het and hom active
        NFkB        = [(19:28) ...    % p65 hom
            (29:38)];      % p65 het
        
        basal_nfkbn = sum(v.OUTPUT2(1,Active_NFkB),2);
        total_nfkb  = sum(v.OUTPUT2(1,NFkB),2);
        RelA = sum(v.OUTPUT2(1,(29:38))) + 2*sum(v.OUTPUT2(1,(19:28)));
        
        max1 = total_nfkb * .05; % threshold (5%)
        min1 = total_nfkb * .01; % min threshold (1%)
        
        if (basal_nfkbn > max1) && (basal_nfkbn < min1)
            constraints(1) = false;
            screen(['C1 failed. Basal NFkBn = ' ...
                num2str((basal_nfkbn/total_nfkb)*100) '%'],v.DISP);
        end
        
        % 8: Basal Hom is .1 to .5 of Total Het
        Hom = sum(v.OUTPUT2(1,19:28),2);
        Het = sum(v.OUTPUT2(1,29:38),2);
        
        if (Hom > Het*.5) && (Hom < Het*.1)
            constraints(8) = false;
            screen(['C8 failed. Basal Hom is = ' ...
                num2str((Hom/Het)*100) '% Het'],v.DISP);
        end
        
        
        % 2: basal total IkBa > IkBb > IkBe
        IkBa = sum(v.OUTPUT2(1,[2 3 21 22 31 32 41 42]),2);
        IkBb = sum(v.OUTPUT2(1,[5 6 23 24 33 34 43 44]),2);
        IkBe = sum(v.OUTPUT2(1,[8 9 25 26 35 36 45 46]),2);
        IkBd = sum(v.OUTPUT2(1,[11 12 27 28 37 38 47 48]),2);
        
        if ((IkBa < IkBe) || (IkBa < IkBb)) && v.IKB_KO(1)
            constraints(2) = false;
            screen(['C2 failed. IkBa/b/e = ' num2str(floor(IkBa)) ', ' ...
                num2str(floor(IkBb)) ', ' num2str(floor(IkBe)) ' nM'],v.DISP);
        end
        
        if (IkBb < IkBe) && v.IKB_KO(2)
            constraints(2) = false;
            screen(['C2 failed. IkBa/b/e = ' num2str(floor(IkBa)) ', ' ...
                num2str(floor(IkBb)) ', ' num2str(floor(IkBe)) ' nM'],v.DISP);
        end
        
        % 3: Free IkB is between 5-25% of total IkB pool
        % For each IkB
        free_IkBa = sum(v.OUTPUT2(1,[2 3]));
        free_IkBb = sum(v.OUTPUT2(1,[5 6]));
        free_IkBe = sum(v.OUTPUT2(1,[8 9]));
        free_IkBd = sum(v.OUTPUT2(1,[11 12]));
        
        min2 = 0.05; % min threshold 5%
        max2 = 0.25; % max threshold 25%
        
        
        constraints(3) = ( ...
            ((free_IkBa < IkBa*max2) && (free_IkBa > IkBa*min2) || ~v.IKB_KO(1))...
            &&((free_IkBb < IkBb*max2) && (free_IkBb > IkBb*min2) || ~v.IKB_KO(2)) ...
            &&((free_IkBe < IkBe*max2) && (free_IkBe > IkBe*min2) || ~v.IKB_KO(3)) ...
            &&((free_IkBd < IkBd*max2) && (free_IkBd > IkBd*min2) || ~v.IKB_KO(4)) ...
            );
        if ~constraints(3)
            screen(['C3 failed. Free IkB a/b/e/d = ' ...
                num2str(floor((free_IkBa/IkBa)*100)) ', ' ...
                num2str(floor((free_IkBb/IkBb)*100)) ', ' ...
                num2str(floor((free_IkBe/IkBe)*100)) ', ' ...
                num2str(floor((free_IkBd/IkBd)*100)) '% ' ...
                ],v.DISP);
        end
        
        %4 IGM or TNFp, activity must be >10fold basal in first hour
        if strcmp(v.IKK2_TYPE,'TNF15p') || strcmp(v.IKK2_TYPE,'IGM')
            nfkbn = sum(v.OUTPUT2(:,Active_NFkB),2);
            fold = nfkbn(1:61)/nfkbn(1);
            constraints(4) = any(fold > 10,1);
            
            if ~constraints(4)
                screen(['C4 failed: 1h NFkBn < 10fold (' num2str(max(fold)) ')'],v.DISP);
            end
        else
            constraints(4) = NaN;
        end
        
        %5 IGM, activity must be > 4 fold at 4hr and at least 20% of peak
        if strcmp(v.IKK2_TYPE,'IGM') && (v.SIM_TIME > 239)
            nfkbn = sum(v.OUTPUT2(:,Active_NFkB),2);
            fold = nfkbn(241)/nfkbn(1);
            peak_to_phase2 = nfkbn(241)/max(nfkbn(1:61));
            constraints(5) = (fold >= 4) && (peak_to_phase2 >= .2);
            
            if ~constraints(5)
                screen(['C5 failed: 4h NFkBn < 4fold (' num2str(max(fold)) ')'],v.DISP);
            end
        else
            constraints(5) = NaN;
        end
        
        %6 IGM, mRNA fold
        %  ikba: 5-50 , ikbb: 1-5, ikbe: 5-50
        if strcmp(v.IKK2_TYPE,'IGM') && (v.SIM_TIME > 120)
            a_fold = max(v.OUTPUT2(1:v.SIM_TIME,1)/ v.OUTPUT2(1,1));
            b_fold = max(v.OUTPUT2(1:v.SIM_TIME,4)/ v.OUTPUT2(1,4));
            e_fold = max(v.OUTPUT2(1:v.SIM_TIME,7)/ v.OUTPUT2(1,7));
            
            constraints(6) = (  (((a_fold >= 5) && (a_fold <= 50)) || ~v.IKB_KO(1))...
                &&(((b_fold >= 1) && (b_fold <= 5)) || ~v.IKB_KO(2))...
                &&(((e_fold >= 5) && (e_fold <= 50))|| ~v.IKB_KO(3))  );
            
            if ~constraints(6)
                screen(['C6 failed: mRNA fold' ' a(' num2str(a_fold) ')' ...
                    ' b(' num2str(b_fold) ')' ...
                    ' e(' num2str(e_fold) ')' ],v.DISP);
            end
            
        else
            constraints(6) = NaN;
        end
        
        %7 Basal IkB levels are +/- 25% of measured abundances
        %    using IkBa, IkBb, IkBe basal values from above
        cell_val = .7; %70nM model / 100k molecules measured = .7
        threshold = .25;% 25% error rate
        a_prot = 400*cell_val; % 400k molecules IkBa
        b_prot = 100*cell_val;
        e_prot = 25*cell_val;
        R_prot = 480*cell_val; % RelA
        
        a_min = a_prot * (1-threshold);
        a_max = a_prot * (1+threshold);
        b_min = b_prot * (1-threshold);
        b_max = b_prot * (1+threshold);
        e_min = e_prot * (1-threshold);
        e_max = e_prot * (1+threshold);
        R_min = R_prot * (1-threshold);
        R_max = R_prot * (1+threshold);
        
        constraints(7) = (  (((IkBa > a_min) && (IkBa < a_max))|| ~v.IKB_KO(1)) ...
            &&(((IkBb > b_min) && (IkBb < b_max))|| ~v.IKB_KO(2)) ...
            &&(((IkBe > e_min) && (IkBe < e_max))|| ~v.IKB_KO(3)) ...
            &&(((RelA > R_min) && (RelA < R_max))|| ~v.MONOMER_KO(1))  );
        
        if ~constraints(7)
            screen(['C7 failed: [Prot]' 'a(' num2str(IkBa/a_prot) ')' ...
                ' b(' num2str(IkBb/b_prot) ')' ...
                ' e(' num2str(IkBe/e_prot) ')' ...
                ' R(' num2str(RelA/R_prot) ')'],v.DISP);
        end
        
    end

%==== Disable the display of 'v' if this run is not part of a loop
if nargin == 0
    v = ['Elapsed time: ' num2str(toc) 'secs'];
end
end

% Display Function
%{
    Input is any variable that is compatible with the disp function
    Will only display if the flag variable is non-zero
%}
function screen(message,flag)
if(flag)
    disp(message);
end
end