function v=setDefaults(name)

% Set up initial data for simulation type
v = {'IKB_KO'      , logical([1;1;1;1]);          % IkBa ; b ; e ; d
    'MONOMER_KO'  , logical([1;1;1;1;1;1]);        %RelA ; p50 ; RelB ; p100 ; p52; cRel
    'DIMER_KO'    , logical([1;1;1;0;1;1;1;1;1]);    % A50; B50; B52; C100; C50; C52; A52; 5050; 5252
    'IKK2multipli', 1;
    'BASAL_IKK2',0.01;
    };

if(strcmp(name,'NEMO'))
    name='custom_LPSc';
end

%param1 = {'ps' 2, 10}; 

%param2 = {'ps' 2, 1}; 

%param3 = {'ps' 2, 0.1};

%param4 = {'ps' 1, 0; 'ps' 2, 0 };  

%species = {[zeros(120,1); 1]}; 

if(strcmp(name,'equil'))
    %Nemo stimulation
    v(size(v,1)+1,:) = {'NIK_INPUT',0};
    v(size(v,1)+1,:) = {'IKK2_TYPE','basal'};
    v(size(v,1)+1,:) = {'NIK_TYPE','basal'};
    v(size(v,1)+1,:) = {'BASAL_IKK2',0.01};
    v(size(v,1)+1,:) = {'IKK2multipli',1};    
    v(size(v,1)+1,:) = {'SIM_TIME',3600*24*2};
elseif(strcmp(name,'NIK'))
    %NIK stimulation
    v(size(v,1)+1,:) = {'NIK_INPUT',1};
    v(size(v,1)+1,:) = {'NIK_TYPE','basal_5hr'};
    v(size(v,1)+1,:) = {'NIK_deg_mod1',10};
    v(size(v,1)+1,:) = {'NIK_deg_mod2',100};
    v(size(v,1)+1,:) = {'IKK2_TYPE','basal'};
    v(size(v,1)+1,:) = {'IKK2multipli',1};
    v(size(v,1)+1,:) = {'BASAL_IKK2',0.03};
    v(size(v,1)+1,:) = {'SIM_TIME',3600*24*2};
elseif(strcmp(name,'BAFF'))
    %NIK stimulation
    v(size(v,1)+1,:) = {'NIK_INPUT',1};
    v(size(v,1)+1,:) = {'NIK_TYPE','basal_5hr'};
    v(size(v,1)+1,:) = {'NIK_deg_mod1', 30/2};
    v(size(v,1)+1,:) = {'NIK_deg_mod2',20};
    v(size(v,1)+1,:) = {'IKK2multipli',3};
    v(size(v,1)+1,:) = {'BASAL_IKK2',0.03};
    v(size(v,1)+1,:) = {'SIM_TIME',3600*24*2};
    v(size(v,1)+1,:) = {'IKK2_TYPE','BAFF'};      
elseif(strcmp(name,'IGM') || strcmp(name,'TNF15p')  || strcmp(name,'LPS45p') || strcmp(name,'custom') || strcmp(name,'Basal_pulse') || strcmp(name,'custom_TNFp15') || strcmp(name,'custom_LPSc') || strcmp(name,'custom_LPSc2'))
    v(size(v,1)+1,:) = {'IKK2_TYPE',name};
    v(size(v,1)+1,:) = {'NIK_TYPE','basal'};
    v(size(v,1)+1,:) = {'IKK2multipli', 3;    };  
    v(size(v,1)+1,:) = {'NIK_INPUT',0};
    v(size(v,1)+1,:) = {'BASAL_IKK2',0.03};
    v(size(v,1)+1,:) = {'SIM_TIME',3600*24*2};
elseif(strcmp(name,'IGM_nfkb1KO'))
    v(size(v,1)+1,:) = {'NIK_INPUT',1};
    v(size(v,1)+1,:) = {'MONOMER_KO'  , logical([1;0;1;1;1;1])}; 
    v(size(v,1)+1,:) = {'NIK_TYPE','basal_multi'};
    v(size(v,1)+1,:) = {'IKK2multipli', 3;    };  
    v(size(v,1)+1,:) = {'NIK_multiplier', 0.3}; 
    v(size(v,1)+1,:) = {'BASAL_IKK2',0.03};
    v(size(v,1)+1,:) = {'SIM_TIME',3600*24*2};
    v(size(v,1)+1,:) = {'IKK2_TYPE','IGM'}; 
elseif(strcmp(name,'IGM_nfkb2het'))
    v(size(v,1)+1,:) = {'NIK_INPUT',1};
    v(size(v,1)+1,:) = {'MONOMER_KO'  , logical([1;1;1;1;1;1])}; 
    v(size(v,1)+1,:) = {'IKK2multipli', 3;    };  
    v(size(v,1)+1,:) = {'NIK_TYPE','basal_multi'};
    v(size(v,1)+1,:) = {'BASAL_IKK2',0.03};
%     v(size(v,1)+1,:) = {pm(1,4), pm(1,4)/4};
    v(size(v,1)+1,:) = {'SIM_TIME',3600*24*2};
    v(size(v,1)+1,:) = {'IKK2_TYPE','IGM'}; 
elseif(strcmp(name,'primingOFF'))
    % Priming test: 'baseline' (without pretreatment)
    v(size(v,1)+1,:) = {'IKK2_TYPE'   , 'basal'; };
    v(size(v,1)+1,:) = {'IKK2multipli', 3;    };                 % Without pretreatment, basal level stays 3x
    v(size(v,1)+1,:) = {'BASAL_IKK2'  , .01*3; };
    v(size(v,1)+1,:) = {'NIK_INPUT'   , 1;};                     % logical : NIK activation or not
    v(size(v,1)+1,:) = {'NIK_TYPE'    ,'basal_delay';};          %'basal_delay', shift 'basal_5hr' to allow it start at NIKstarttime
    v(size(v,1)+1,:) = {'NIKdelay'    , 3*60+8*60;};
    v(size(v,1)+1,:) = {'SIM_TIME',3600*48};
elseif(strcmp(name,'primingON'))
    % Priming test: with pretreatment
    v(size(v,1)+1,:) = {'IKK2_TYPE'   , 'Basal_pulse';};         %'Basal_pulse' - Step input of NEMO
    v(size(v,1)+1,:) = {'NEMOpulse'   , 30;};                    % 30x basal as pulse amplitute
    v(size(v,1)+1,:) = {'Pulsetime'   , 8*60} ;                % pulse length== 8 hrs (convert to mins)
    v(size(v,1)+1,:) = {'BASAL_IKK2'  , .01*3; };                %
    v(size(v,1)+1,:) = {'NIK_INPUT'   , 1;} ;                    % logical : NIK activation or not
    v(size(v,1)+1,:) = {'NIK_TYPE'    ,'basal_delay';};          %'basal_delay', shift 'basal_5hr' to allow it start at NIKstarttime
    v(size(v,1)+1,:) = {'NIKdelay'    , 3*60+8*60;} ;            % delay time == 3 hrs (convert to mins)
    v(size(v,1)+1,:) = {'SIM_TIME',3600*48};
elseif(strcmp(name,'priming'))
    % custom: run Allison's pre-treatment and only save the output
    % as ratio of 24 hour post-stimulus. Also this will ignore the 'customParams' except for RATE_MODS, I believe
elseif(strcmp(name, 'costim'))
    %costimulation with NIK input and IKK input 
    v(size(v,1)+1,:) = {'NIK_INPUT',1};
    v(size(v,1)+1,:) = {'NIK_TYPE','basal_5hr'};
    v(size(v,1)+1,:) = {'NIK_deg_mod1',30/2};
    v(size(v,1)+1,:) = {'NIK_deg_mod2',20};
    v(size(v,1)+1,:) = {'BASAL_IKK2',0.03};
    v(size(v,1)+1,:) = {'SIM_TIME',3600*24*2};
    v(size(v,1)+1,:) = {'IKK2_TYPE','costim'};  
    
elseif(strcmp(name, 'low'))
    v(size(v,1)+1,:) = {'NIK_INPUT',1};
    v(size(v,1)+1,:) = {'NIK_TYPE','basal_5hr'};
    v(size(v,1)+1,:) = {'NIK_deg_mod1', 30/2};
    v(size(v,1)+1,:) = {'NIK_deg_mod2',20};
    v(size(v,1)+1,:) = {'IKK2_TYPE','basal'};
    v(size(v,1)+1,:) = {'BASAL_IKK2',0.01};
    v(size(v,1)+1,:) = {'IKK2multipli',0.01};    
    v(size(v,1)+1,:) = {'SIM_TIME',3600*24*2};
    
elseif(strcmp(name, 'med'))
    v(size(v,1)+1,:) = {'NIK_INPUT',1};
    v(size(v,1)+1,:) = {'NIK_TYPE','basal_5hr'};
    v(size(v,1)+1,:) = {'NIK_deg_mod1', 30/2};
    v(size(v,1)+1,:) = {'NIK_deg_mod2',20};
    v(size(v,1)+1,:) = {'IKK2_TYPE','basal'};
    v(size(v,1)+1,:) = {'BASAL_IKK2',0.01};
    v(size(v,1)+1,:) = {'IKK2multipli',1};    
    v(size(v,1)+1,:) = {'SIM_TIME',3600*24*2};
    
elseif(strcmp(name, 'high'))
    v(size(v,1)+1,:) = {'NIK_INPUT',1};
    v(size(v,1)+1,:) = {'NIK_TYPE','basal_5hr'};
    v(size(v,1)+1,:) = {'NIK_deg_mod1', 30/2};
    v(size(v,1)+1,:) = {'NIK_deg_mod2',20};
    v(size(v,1)+1,:) = {'IKK2_TYPE','basal'};
    v(size(v,1)+1,:) = {'BASAL_IKK2',0.01};
    v(size(v,1)+1,:) = {'IKK2multipli',100};    
    v(size(v,1)+1,:) = {'SIM_TIME',3600*24*2};
    
elseif(strcmp(name, 'lowNIK'))
    v(size(v,1)+1,:) = {'NIK_INPUT',1};
    v(size(v,1)+1,:) = {'NIK_TYPE','basal_0hr'};
    v(size(v,1)+1,:) = {'NIK_deg_mod1', 1};
    v(size(v,1)+1,:) = {'NIK_deg_mod2',20};
    v(size(v,1)+1,:) = {'IKK2_TYPE','IGM'};
    v(size(v,1)+1,:) = {'BASAL_IKK2',0.01};
    v(size(v,1)+1,:) = {'IKK2multipli',1};    
    v(size(v,1)+1,:) = {'SIM_TIME',3600*24*2};
    %v(size(v,1)+1,:) = {'RATE_MODS',param1};
    
elseif(strcmp(name, 'medNIK'))
    v(size(v,1)+1,:) = {'NIK_INPUT',1};
    v(size(v,1)+1,:) = {'NIK_TYPE','basal_0hr'};
    v(size(v,1)+1,:) = {'NIK_deg_mod1', 1};
    v(size(v,1)+1,:) = {'NIK_deg_mod2',20};
    v(size(v,1)+1,:) = {'IKK2_TYPE','IGM'};
    v(size(v,1)+1,:) = {'BASAL_IKK2',0.01};
    v(size(v,1)+1,:) = {'IKK2multipli',1};    
    v(size(v,1)+1,:) = {'SIM_TIME',3600*24*2};
    %v(size(v,1)+1,:) = {'RATE_MODS',param2}; 
    
elseif(strcmp(name, 'highNIK'))
    v(size(v,1)+1,:) = {'NIK_INPUT',1};
    v(size(v,1)+1,:) = {'NIK_TYPE','basal_0hr'};
    v(size(v,1)+1,:) = {'NIK_deg_mod1', 1};
    v(size(v,1)+1,:) = {'NIK_deg_mod2',20};
    v(size(v,1)+1,:) = {'IKK2_TYPE','IGM'};
    v(size(v,1)+1,:) = {'BASAL_IKK2',0.01};
    v(size(v,1)+1,:) = {'IKK2multipli',50};    
    v(size(v,1)+1,:) = {'SIM_TIME',3600*24*2};
    %v(size(v,1)+1,:) = {'RATE_MODS',param3}; 
    
elseif(strcmp(name, 'varyNIK'))
    v(size(v,1)+1,:) = {'NIK_INPUT',1};
    v(size(v,1)+1,:) = {'NIK_TYPE','basal_0hr'};
    v(size(v,1)+1,:) = {'NIK_deg_mod1', 1};
    v(size(v,1)+1,:) = {'IKK2_TYPE','IGM_mod'};
    v(size(v,1)+1,:) = {'BASAL_IKK2',0.01};
    v(size(v,1)+1,:) = {'IKK2multipli',1};    
    v(size(v,1)+1,:) = {'SIM_TIME',3600*24*2};
    %v(size(v,1)+1,:) = {'RATE_MODS',param4}; 
    %v(size(v,1)+1,:) = {'SPECIES',species};  
    
%         v(size(v,1)+1,:) = {'NIK_INPUT',1};
%     v(size(v,1)+1,:) = {'NIK_TYPE','basal_5hr'};
%     v(size(v,1)+1,:) = {'NIK_deg_mod1', 30/2};
%     v(size(v,1)+1,:) = {'NIK_deg_mod2',20};
%     v(size(v,1)+1,:) = {'IKK2_TYPE','basal'};
%     v(size(v,1)+1,:) = {'BASAL_IKK2',0.01};
%     v(size(v,1)+1,:) = {'IKK2multipli',1};    
%     v(size(v,1)+1,:) = {'SIM_TIME',3600*24*2};

elseif(strcmp(name, 'varyboth'))
    v(size(v,1)+1,:) = {'NIK_INPUT',1};
    v(size(v,1)+1,:) = {'NIK_TYPE','basal_0hr'};
    v(size(v,1)+1,:) = {'NIK_deg_mod1', 1};
    v(size(v,1)+1,:) = {'IKK2_TYPE','basal'};
    v(size(v,1)+1,:) = {'BASAL_IKK2',0.01};
    v(size(v,1)+1,:) = {'IKK2multipli',1};    
    v(size(v,1)+1,:) = {'SIM_TIME',3600*24*2};
    v(size(v,1)+1,:) = {'RATE_MODS',param4}; 
    v(size(v,1)+1,:) = {'SPECIES',species}; 

    
end

end