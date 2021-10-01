% B-Cell Model - Jon Almaden IkBd crosstalk in B-cell proliferation 
% adapted from Jeremy and Vincent's model 
% adapted from Bryce's model 031913 (Rachel) 
% adds p50 homodimer and p52 homodimer 032613 (Rachel)


%% passParams File:
% This will add in A52, C50, and C52 for Jon's B-cell simulation 
% Test version to figure out the EMSA activity situation 
% 011614 This will be modified to make RelB and p50 both NFkB inducible and
% B52 can also be bound by IkBd, assuming the same params as IkBd:B50
% 030604 Modified to make RelB more dependent on cRel than RelA dimers,
% in order to have high cRel or low cRel abundance when running simulations
% NOTE: basal has same levels as 011614 version

function v=passParams030614_tcell(v)
         v.DIMERS = 9;
         v.DIMERStext = {'A50','B50','B52','C100', 'C50', 'C52', 'A52', '5050', '5252'};
         v.numSpecies = v.IKB_COMP*v.IKBS + v.MONOMER_COMP*v.MONOMERS + v.DIMER_COMP*v.DIMERS + v.OTHER_COMP;
         v.SPECIES{1}  = zeros(v.numSpecies,1); % all species start at zero [] by default
                     
        v.SPECIES{2}  = { ...  % Component # 1-12 (IkB) , 13-30 (NFkB mon) , 31-40 (A50), 
                               %41-50 (B50), 51-60 (B52), 61-70 (C100), 71-80 (C50), 81-90 (C52), 91-100 (A52), 
                               % 101-110 (5050), 111-120 (5252), 121 (NIK)
            'tIkBa','IkBa', 'IkBan','tIkBb','IkBb','IkBbn','tIkBe','IkBe','IkBen','tIkBd','IkBd','IkBdn',...
            'tRelA', 'RelA','RelAn','tp50','p50','p50n', 'tRelB', 'RelB','RelBn','tp100','p100','p100n','tp52','p52','p52n', 'tcRel', 'cRel', 'cReln'...
            'A50','A50n','IkBaA50','IkBaA50n','IkBbA50','IkBbA50n','IkBeA50','IkBeA50n','IkBdA50','IkBdA50n',...
            'B50','B50n','IkBaB50','IkBaB50n','IkBbB50','IkBbB50n','IkBeB50','IkBeB50n','IkBdB50','IkBdB50n',...
            'B52','B52n','IkBaB52','IkBaB52n','IkBbB52','IkBbB52n','IkBeB52','IkBeB52n','IkBdB52','IkBdB52n',...
            'C100','C100n','IkBaC100','IkBaC100n','IkBbC100','IkBbC100n','IkBeC100','IkBeC100n','IkBdC100','IkBdC100n',...
            'C50','C50n','IkBaC50','IkBaC50n','IkBbC50','IkBbC50n','IkBeC50','IkBeC50n','IkBdC50','IkBdC50n',...
            'C52','C52n','IkBaC52','IkBaC52n','IkBbC52','IkBbC52n','IkBeC52','IkBeC52n','IkBdC52','IkBdC52n',...
            'A52','A52n','IkBaA52','IkBaA52n','IkBbA52','IkBbA52n','IkBeA52','IkBeA52n','IkBdA52','IkBdA52n',...
            '5050', '5050n', 'IkBa5050', 'IkBa5050n', 'IkBb5050', 'IkBb5050n', 'IkBe5050', 'IkBe5050n', 'IkBd5050', 'IkBd5050n',...
            '5252', '5252n', 'IkBa5252', 'IkBa5252n', 'IkBb5252', 'IKBb5252n', 'IkBe5252', 'IkBe5252n', 'IkBd5252', 'IkBd5252n',...
            'NIK' ...
            }';
        
        v.masterList = masterList2(v.SPECIES{2});
    function out=masterList2(SPECIES)
        
        
        out = struct('p100_total', [11 12 23 24 39 40 49 50 59 60 61:70 79 80 89 90 99 100], ... 
            'NFkB_dimer_total_RelB',[(41:70)], ...
            'IkBd_protein_RelB', [49 50 59 60 61:70], ...  % includes p100:B dimers
            'IkBa_protein_RelB', [43 44 53 54 63 64], ...            
            'IkBe_protein_RelB', [47 48 57 58 67 68], ...
            'Total_RelB_Dimers_cyt', [41:2:70], ...
            'Total_RelB_Dimers_nuc', [42:2:60], ...
            'Total_RelB_cyt',[20 (41:2:70)], ...
            'Total_RelB_nuc',[21 (42:2:70)], ...
            'Total_RelA_dimer_cyt',[31:2:40 91:2:100], ...
            'Total_RelA_dimer_nuc',[32:2:40 92:2:100], ...
            'NFkB_protein_total_RelB', [20 21 (41:70)], ...
            'NFkB_dimer_total_B50', [41:50], ...
            'NFkB_dimer_total_B52', [51:60], ...
            'NFkB_dimer_total_C100',[61:70], ...
            'IkBd_protein_total', [11 12 39 40 49 50 59 60 69 70 79 80 89 90 99 100], ... 
            'IkBa_protein_total', [2 3 33 34 43 44 53 54 63 64 73 74 83 84 93 94], ...            
            'IkBb_protein_total', [5 6 35 36 45 46 55 56 65 66 75 76 85 86 95 96], ...            
            'IkBe_protein_total', [8 9 37 38 47 48 57 58 67 68 77 78 87 88 97 98],...
            'NFkB_protein_total_p52', [26  27 51:60 81:90 91:100],...
            'NFkB_dimer_total_cRel', [71:90],...
            'IkBd_protein_cRel', [79 80 89 90],...
            'IkBe_protein_cRel', [77 78 87 88],...
            'Total_cRel_dimers_cyt', [71:2:90],...
            'Total_cRel_dimer_nuc', [72:2:90],...
            'Total_cRel_cyt', [29 71:2:90],...
            'Total_cRel_nuc', [30 72:2:90],...
            'NFkB_dimer_total_C50', [71:80],...
            'NFkB_dimer_total_C52', [81:90],...
            'NFkB_protein_total_RelA', [14 15 31:40 91:100],...
            'NFkB_protein_total_cRel', [29 30 71:90]);
          
       % For each species in SPECIES, add this to the struct:
       for s=1:length(SPECIES)
           name = SPECIES{s};  
           
           %add 'P' to the beginning of the proteins (anything not t*)
           start = name(1);
           if (~strcmp('t',start)) 
               name = strcat('P',name);
           end
           
           out.(name)=s;
       end
    end
        
         
        %{
        Parameters use seconds as their time unit
         *t can put units into minutes
            % Association : 1 / (nM * sec)
            % Dissociation : 1 / (sec)
            % mRNA and Protein Deg : 1 / (sec)
            % IKK Mediated Deg : 1 / (nM * sec)
            % Localization : 1 / (sec)
            % const txn rate: nM / sec
            % inducible rate : unitless
        To simulate in minutes, the following code scales them
        %}
        if v.SEC_RUN % Simulate in seconds
            t = 1;  % no change. Rates are by default in seconds
        else
            t = 60; % Convert to minutes (60sec/min)
        end
        
        %==== Initialize Parameter Lists
        v.pi         = zeros(16,v.IKBS);               % IkBs
        v.pm         = zeros(10,v.MONOMERS);           % Monomers
        v.pd         = zeros(12,v.DIMERS);             % Dimer-only
        v.pdi        = zeros(7,v.DIMERS,v.IKBS);       % Dimer with IKB
        v.pdm        = zeros(1,v.DIMERS,v.MONOMERS);   % Dimer with Monomer
        v.ps         = zeros(5,1);                     % Special reactions
        
        %==== IkB reactions  % 1-14 ikba, 15-28 b, 29-42 e, 43-56 d
        %           ikba      ikbb    ikbe    ikbd
        v.pi(1,:) = [2e-5*4;  1e-5*2;  3e-7*4;    0     ]*t; % const txn
        v.pi(2,:) = [4.8e-4;  4.8e-5; 6.4e-5; 0     ]*t; % mRNA deg
        v.pi(3,:) = [.2;      .2;     .2;     0     ]*t; % tsl rate
        
        v.pi(4,:) = [2e-3;    2e-3;   1.925e-4;    2e-4/4]*t; % prot deg cyt
        v.pi(5,:) = v.pi(4,:);                           % prot deg nuc
        
        v.pi(6,:) = [1e-3;    1.5e-4; 7.5e-4; 7.5e-4]*t; % free import
        v.pi(7,:) = [2e-4;    2e-4;   2e-4;   2e-4  ]*t; % free export
        
        v.pi(8,:) = [900;     2700;   2700;   0       ]/t; % i-txn delay
        v.pi(9,:) = [1.1;     1.1;      1.1;      1       ];   % Hill Coeff.
        v.pi(10,:)= [150;   150;    150;    1       ];   % NFkB Kd (nM)
        
        %v.pi(11,:)= [1.5e-5;  5e-6; 3.75e-6;   0       ]*t*1.5; % free IKK2 mediated deg
        v.pi(11,:)= [1.5e-5;  1.5e-5; 1.5e-5;   0       ]*t*1.5; % free IKK2 mediated deg
        v.pi(12,:)= v.pi(11,:);                              % bound IKK2 mediated deg
        
        v.pi(13,:)= [4e-6;  4e-6;   4e-6;   4e-6    ]*t;  % bound prot deg cyt
        v.pi(14,:)= v.pi(13,:);                           % bound prot deg nuc
        
        v.pi(15,:)= [0;     0;      0;      2e-2/2  ]*t;  % free NIK mediated deg
        v.pi(16,:)= v.pi(15,:);                         % bound NIK mediated deg
        
        %==== NFkB monomer reactions
        %              RelA         p50       RelB      p100     p52     cRel
        %figures until June May 2019
        %v.pm(1,:) 	= [1.2e-6*1.7;  5.76e-08*1.04*1.04;  1.5641e-08*10;  4.0104e-07/20;  0;  6e-9]*t.*5;%/8 to make cRel lower; % const txn relb was 3e-8
        v.pm(1,:) 	= [0.6e-6*1.7;  5.76e-08*1.04*1.04;  1.5641e-08*10;  4.0104e-07/20;  0;  6e-9]*t.*5;%/8 to make cRel lower; % const txn relb was 3e-8
        v.pm(2,:) 	= [4.8e-5;  4.8e-5; 4.8e-5; 1.6e-5; 0;      1.6e-5]*t; % mRNA degradation
        v.pm(3,:)	= [.2;      .2;     .2;     .2;     0;      .2]*t; % translation rate
        
        v.pm(4,:)	= [3.8e-4;  3.8e-4; 3.8e-4; 3.8e-4; 3.8e-4; 3.8e-4]*t; % free prot deg cyt
        v.pm(5,:)	= v.pm(4,:);             % free prot deg nuc 2hr
        
        v.pm(6,:)   = [3600;    3600;   3600;   3600*4;   0;  3600*12 ]/t; % inducible txn delay
        v.pm(7,:)   = [1;       1;      1;      1;      1;  1  ]; % Hill Coefficient
        v.pm(8,:)	= [150;     150;    150;     150;     1;  150]; % NFkB Kd for txn limit (nM)
        v.pm(9,:)	= [0;       0;      0;      5e-2;   0;  0  ]; % NIK processing rate
        v.pm(10,:)	= [0;       0;      0;      5;      0;  0  ]; % NIK processing target index
        
        %==== NFkB Dimer-only reactions
        %          A50    B50     B52     C100     C50      C52   A52  5050   5252
        %       
        v.pd(1,:)    = [1;      2;      3;      6;      6;      6;      1; 2; 5];% Mon #1 index
        v.pd(2,:)    = [2;      3;      5;      4;      2;      5;      5; 2; 5];% Mon #2 index
        %v.pd(3,:)    = [3.16e-5; 5.06e-6; 1.6e-5; 1.6e-5; 1.6e-5; 1.6e-5; 1.6e-5; 3e-5; 3e-5;]*t;  % mon asn cyt : B50 starts 10-fold less than RelA:p50 Kd. (here is 1/sqrt(10), dissoc is 1*sqrt(10))
        v.pd(3,:)    = [3.16e-5; 1.6e-5; 1.6e-5; 1.6e-5; 1.6e-5; 1.6e-5; 1.6e-5; 3e-5; 3e-5;]*t;  % mon asn cyt : B50 starts 10-fold less than RelA:p50 Kd. (here is 1/sqrt(10), dissoc is 1*sqrt(10))
    
        v.pd(4,:)    = v.pd(3,:);                   % mon asn nuc
        v.pd(5,:)    = [3.16e-4; 7.59e-5; 2.4e-4;  2.4e-5; 2.4e-5; 2.4e-5; 3.16e-4*50;9e-4; 9e-4;]*t;% dimer dsn cyt
        v.pd(6,:)    = [3.16e-4/10; 7.59e-5/10; 2.4e-4/10;  7.59e-6; 2.4e-5/10; 2.4e-5/10; 3.16e-4/10; 9e-4; 9e-4;]*t; % dimer dsn nuc
        
        v.pd(7,:)    = [9e-2;   9e-2;    9e-2;   8e-5;  9e-2;   9e-2;   9e-2; 9e-2;   9e-2;]*t;  % nuclear import % import/export of C100 is equal ... no strong NLS is the assumption!
        v.pd(8,:)    = [8e-5;   8e-5;    8e-5;   8e-5;   8e-5;   8e-5;  8e-5; 8e-5;   8e-5;]*t;  % nuclear export
        
        v.pd(9,:)    = [4e-6;   4e-6;    4e-6; 4e-6; 4e-6; 4e-6; 4e-6; 4e-6; 4e-6]*t;  % Free dimer Deg cyt 
        v.pd(10,:)   = v.pd(9,:);                   % Free dimer Deg nuc
        v.pd(11,:)   = v.pd(9,:);                   % Bound dimer Deg cyt
        v.pd(12,:)   = v.pd(9,:);                   % Bound dimer Deg nuc
        
        %==== NFkB Dimer-IKB reactions  (split up by IkB isoform)
        %           A50    B50     B52     C100     C50      C52    A52   5050   5252
         v.pdi(1,:,:) = [5.01e-5  5.01e-5 0     0     5.01e-5     5.01e-5 5.01e-5 0 0;  ...% Ikba : Cyt Asn 
                        3.55e-6  0       0     0     3.55e-6      3.55e-6  3.55e-6 0 0; ...% Ikbb  
                        2.24e-5 2.24e-5  0     0     2.24e-5     2.24e-5 2.24e-5 0 0;...% Ikbe  
                        1e-5     1e-6    1e-5     0     1e-5      3.3e-4 1e-5 0 0;]'*t;% Ikbd  
                                
        v.pdi(3,:,:) = [1e-5  1e-5    0     0    8e-5       8e-5     1e-5 0 0; ... % Ikba : Cyt Dsn 
                        2.82e-4   0      0     0    2.82e-4    2.82e-4   2.82e-4 0 0; ... % Ikbb  
                        1e-4  1e-4   0     0    4.47e-7    4.47e-7  1e-4 0 0; ... % Ikbe  
                        1.4e-5   1.4e-5  1.4e-5     0    1.4e-5     1.4e-5   1.4e-5 0 0; ]'*t;% Ikbd 
        
        v.pdi(2,:,:) = v.pdi(1,:,:); % nuclear asn
        v.pdi(4,:,:) = v.pdi(3,:,:); % nuclear dsn
        
        %           A50    B50     B52     C100     C50      C52    A52   5050   5252
        v.pdi(5,:,:) = [4.6e-3  4.6e-3   0      0       4.6e-3  4.6e-3  4.6e-3 0 0;...  % Ikba : Nuc Import
                        4.6e-4  4.6e-4   0      0       4.6e-4  4.6e-3  4.6e-3 0 0;...  % Ikbb  
                        2.3e-3  2.3e-3   0      0       2.3e-3  2.3e-3  2.3e-3 0 0;...  % Ikbe
                        4.6e-4  4.6e-4   0      0       4.6e-4  4.6e-4  4.6e-4 0 0;]'*t;% Ikbd
        
        v.pdi(6,:,:) = [1.4e-2  1.4e-2   0      0       1.4e-2   1.4e-2  1.4e-2 0 0;...  % Ikba : Nuc Export 
                        7e-3    7e-3     0      0       7e-3     7e-3    7e-3 0 0;...    % Ikbb
                        7e-3    7e-3     0      0       7e-3     7e-3    7e-3 0 0;...    % Ikbe
                        7e-3    7e-3     0      0       7e-3     7e-3    7e-3 0 0;]'*t;   % Ikbd
        
        v.pdi(7,:,:) = [200      0       0       0       0       0     200 0 0;...  % Ikba : I-Txn
                         0       0       0       0       0       0       0 0 0;...   % Ikbb
                         5      0       0       0       25    25     5 0 0;... % Ikbe
                         0       0       0       0       0       0       0 0 0]';    % Ikbd
        
        %==== NFkB Dimer-Monomer synthesis reactions
        %               A50    B50     B52     C100     C50      C52   A52
        v.pdm(1,:,:) = [0      0       0        0       0       0        0 0 0;...  % RelA
                        20    0       0        0      20      20     20 0 0;...  % p50
                        5   0       0        0       10*10       10*10    5 0 0;... % RelB
                        1000    0       0        0      1500     1500     1000 0 0;... % p100
                        0      0       0        0       0       0        0 0 0;...  % p52
                        200    0       0        0       200     200    200 0 0]'*1; % cRel
                v.pdm(1,:,end)=v.pdm(1,:,end)./10;

        %==== Special reactions
        % NIK synthesis and degradation (combined txn and tsl)
        v.ps(1) = 7e-4*t;   % NIK synthesis
        v.ps(2) = 7.7e-4*t./5; % NIK degradation (15 min half life)
        v.ps(3) = 10;       % Km of NIK activity ... rate*(NIK*p100)/(Km+p100) processing rate of p100
        v.ps(4) = 2e-4*t; % p100 homodimer/association rate to form IkBd
        v.ps(5) = 2e-4*t;   % p100 homodimer/dissocation rate from IkBd
        v.ps(6) = 10;      % Km of NIK activity ... rate*(NIK*IkBd)/(Km+IkBd) degradation of IkBd 
        
          

end