function delta = nfkb_ode_62_mod(t,x,options,v) %#ok<INUSL>
  % input arguments:
  %    options: ode15s solver options (usually an empty array [])
  %    v : struct of simulation variables
  
%--- Set the Persistent variables
%{
   The sim file calls the ode function once without arguments to reset
    the persistent variables. Once done, it returns.
   nfkb_ode_62(
%}
persistent DELAY_IKB;
persistent DELAY_MON;
persistent MON_DIMS;
persistent IKB_DIMS;

global FLUX_CELL;
global FLUX_VEC;
global FLUX_T;
global timeStepCount;

persistent FLUX_IKB;
persistent FLUX_MON;
persistent FLUX_DIMER;
persistent FLUX_TRIMER;
persistent FLUX_SPECIAL;

%v.pdm(1,:,4)=v.pdm(1,:,4)./100000;
defaultSize = 2000; % Used for increments of fluxes to keep track of so that I don't have to allocate memory every loop

if isempty(t) 
%         disp('empty');
    % Create arrays for inducible transcription delays
    %{ 
     3 columns in each cell array
        1 : array of species concentration
        2 : array of time points for each concentration
        3 : index of last element assigned (iterated in ODE)
    %}
    DELAY_IKB          = cell(v.IKBS,3); 
    DELAY_IKB(:,[1 2]) = {zeros(1000,1)}; %increase if getting errors
    DELAY_IKB(:,3)     = {1};
    
    DELAY_MON          = cell(v.MONOMERS,3);
    DELAY_MON(:,[1 2]) = {zeros(1000,1)}; %increase if getting errors
    DELAY_MON(:,3)     = {1};
   
    FLUX_IKB = zeros(v.IKBS,10);
    FLUX_MON = zeros(v.MONOMERS,7); % technically there should be a value here for induced expression, but they're all 0 in this model so not included
    FLUX_DIMER = zeros(v.DIMERS,8);
    FLUX_TRIMER = zeros(v.DIMERS * v.IKBS ,12);%  even though 5050 is all zeros
    FLUX_SPECIAL = zeros(1,4);% zeros(1,1);  %p100 dimerization,  NIK activity
       
    FLUX_VEC = zeros(defaultSize,numel(FLUX_IKB) + numel(FLUX_MON) + numel(FLUX_DIMER) + numel(FLUX_TRIMER) + numel(FLUX_SPECIAL));
    FLUX_T = zeros(1,defaultSize);
    timeStepCount = 0;
    
    % Create Arrays of which dimers induce which monomers or IkBs 
    %{
    This prevents dimers that do not induce txn from becoming repressors
    (when a dimer doesn't induce txn, its inducible rate = 0.  This
    makes the term for that dimer vanish from the numerator, but it
    still remains in the denominator [see the functional form below])
    %}
    MON_DIMS = cell(v.MONOMERS,1);
    IKB_DIMS = cell(v.IKBS,1);
    
    for i = 1:v.MONOMERS
        MON_DIMS{i} = find(v.pdm(1,:,i) > 0);
    end
    
    for i = 1:v.IKBS
        IKB_DIMS{i} = find(v.pdi(7,:,i) > 0);
    end
    
    return; 
end
  
%--- Load Previous Concentrations and set iteration changes to zero
ikb             = zeros(v.IKB_COMP,v.IKBS);
monomer         = zeros(v.MONOMER_COMP,v.MONOMERS);
dimer           = zeros(v.DIMER_COMP,v.DIMERS);
other           = zeros(v.OTHER_COMP,1);

ikb(:)          = x( 1:numel(ikb));
monomer(:)      = x((1+numel(ikb)):(numel(ikb)+numel(monomer)));
dimer(:)        = x((1+numel(ikb)+numel(monomer)):(numel(ikb)+numel(monomer)+numel(dimer)));
other(:)        = x((1+numel(ikb)+numel(monomer)+numel(dimer)):numel(x));

delta_ikb       = zeros(v.IKB_COMP,v.IKBS);
delta_monomer   = zeros(v.MONOMER_COMP,v.MONOMERS);
delta_dimer     = zeros(v.DIMER_COMP,v.DIMERS);
delta_other     = zeros(v.OTHER_COMP,1);

timeStepCount = timeStepCount + 1;

%--- Calculate NIK and IKK2 Activity 
%{
 IKK2 activity is the product of total IKK2 and percent active
    Max activity is total IKK2 (100% active * TOTAL_IKK2)
    Basal activity is defined by the v.BASAL_IKK2 variable
        Typically, we have been using 1% basal activity
%}
if v.PHASE == 1
    ikk_activity = v.BASAL_IKK2 * v.TOTAL_IKK2;
    %======Added 9/9/2010 for Vincent======================
    % Jeremy modified.
    if (length(v.NIK_TYPE) >= 11 && strcmp(v.NIK_TYPE(1:11),'basal_multi'))
        %nik_deg_mod = v.NIK_multiplier;
        nik_deg_mod=v.NIK_mod(1);
    else
    %previous:
    %nik_deg_mod = 1; % 1x degradation rate    
    nik_deg_mod=v.NIK_mod(1);
    end
else
     if t >= v.SIM_TIME
            ikk_activity = v.IKK2(end);
            nik_deg_mod  = v.NIK_mod(end);
     else
        % IKK = mx + b ... linear over 1 minute interval 
        start = floor(t)+1; % +1 to adjust for t=0 index of 1
        stop  = start + 1;
        
        ikk_activity = ((v.IKK2(stop) - v.IKK2(start)) ...
            * (t-floor(t)) + v.IKK2(start))* v.TOTAL_IKK2;
        
        nik_deg_mod  = ((v.NIK_mod(stop) - v.NIK_mod(start)) ...
            * (t-floor(t)) + v.NIK_mod(start));
     end
end
   
%---- NFkB Dimer Reactions
%{
    A for loop is used to calculate fluxes for each dimer
    A second, nested, for loop is used to calculate dimer:IkB reactions
%}
for i = 1:v.DIMERS
    if ~v.DIMER_KO(i) %Skip if dimer is Knocked-Out
        continue;
    end

    mon1 = v.pd(1,i); % index# of the 1st monomer
    mon2 = v.pd(2,i); % index# of the 2nd monomer 

    % Association  mon1 + mon2 --> NFkB
    flux1 = v.pd(3,i) * monomer(2,mon1) * monomer(2,mon2); % cytoplasm
    flux2 = v.pd(4,i) * monomer(3,mon1) * monomer(3,mon2); % nucleus

    % Dissociation NFkB --> mon1 + mon2
    flux3 = v.pd(5,i) * dimer(1,i); % cytoplasm
    flux4 = v.pd(6,i) * dimer(2,i); % nucleus

    % Localization NFkB <--> NFkBn
    flux5 = v.pd(7,i) * dimer(1,i); % import          
    flux6 = v.pd(8,i) * dimer(2,i); % export
    
    % NFkB dimer degradation
    flux7 = v.pd(9,i) * dimer(1,i); % cytoplasm    
    flux8 = v.pd(10,i)* dimer(2,i); % nucleus 
    
    % Update concentrations with fluxes
    delta_dimer(1,i)= delta_dimer(1,i)+flux1-flux3-flux5+flux6-flux7; % NFkB cyt
    delta_dimer(2,i)= delta_dimer(2,i)+flux2-flux4+flux5-flux6-flux8; % NFkB nuc
    
    delta_monomer(2,mon1)= delta_monomer(2,mon1) -flux1+flux3; % mon1 cyt
    delta_monomer(2,mon2)= delta_monomer(2,mon2) -flux1+flux3; % mon2 cyt
    
    delta_monomer(3,mon1)= delta_monomer(3,mon1) -flux2+flux4; % mon1 nuc
    delta_monomer(3,mon2)= delta_monomer(3,mon2) -flux2+flux4; % mon2 nuc
    
    FLUX_DIMER(i,:) = [flux1 flux2 flux3 flux4 flux5 flux6 flux7 flux8];
        
    % IKB:Dimer Reactions             
    for j = 1:v.IKBS  
        if ~v.IKB_KO(j) %Skip if the IkB is knocked out
            continue;
        end

        cyt = 1+(j*2);    % calculate index for each cyt IkB:NFkB
                          %  The 1st cyt IkB:Dimer is index 3 
        nuc = 1 + cyt;    % calculate index for each nuc IkB:NFkBn
                          %  The corresponding nuc IkB:Dimer is cyt+1

        % Association   IkB + NFkB --> IkB:NFkB
        flux1 = v.pdi(1,i,j) * dimer(1,i) * ikb(2,j); % cytoplasm
        flux2 = v.pdi(2,i,j) * dimer(2,i) * ikb(3,j); % nucleus

        % Dissociation  IkB:NFkB --> IkB + NFkB
        flux3 = v.pdi(3,i,j) * dimer(cyt,i); % cytoplasm
        flux4 = v.pdi(4,i,j) * dimer(nuc,i); % nucleus  

        % Localization.  IkB:NFkB <==> IkB:NFkBn  
        flux5 = v.pdi(5,i,j) * dimer(cyt,i); % import   
        flux6 = v.pdi(6,i,j) * dimer(nuc,i); % export   

        % IKK2-mediated IkB degradation.  IkB:NFkB --> NFkB
            % Using an IkB-list parameter
        flux7 = v.pi(12,j) * dimer(cyt,i) * ikk_activity;
        
        % NIK-mediated IkB degradation.  IkB:NFkB --> NFkB
            % Using an IkB-list parameter
        %flux7b = v.pi(16,j) * (dimer(cyt,i)*other(1))/(v.ps(6)+dimer(cyt,i)); % saturation
        flux7b = v.pi(16,j)*dimer(cyt,i)*other(1)/(v.ps(6)*(1+monomer(2,4)/v.ps(3))+dimer(cyt,i));
                % v.pi(16,j) *  dimer(cyt,i) * other(1); % no saturation 

        % IKK-indepedent IkB degradation.  IkB:NFkB --> NFkB
            % Using an IkB-list parameter
        flux8 = v.pi(13,j) * dimer(cyt,i); % cytoplasm
        flux9 = v.pi(14,j) * dimer(nuc,i); % nucleus
        
        % IkB-bound NFkB Dimer Deg.  IkB:Dimer --> IkB
        flux10 = v.pd(11,i) * dimer(cyt,i); % cytoplasm
        flux11 = v.pd(12,i) * dimer(nuc,i); % nucleus 
            
        % Update concentrations with fluxes
        delta_dimer(1,i)= delta_dimer(1,i) -flux1+flux3+flux8+flux7+flux7b;% NFkB
        delta_dimer(2,i)= delta_dimer(2,i) -flux2+flux4+flux9;      % NFkBn
        
        delta_dimer(cyt,i)= delta_dimer(cyt,i) ... % IkB:NFkB
                            +flux1-flux3-flux5+flux6-flux7-flux7b-flux8-flux10;
        delta_dimer(nuc,i)= delta_dimer(nuc,i) ... % IkB:NFkBn
                            +flux2-flux4+flux5-flux6-flux9-flux11;
                                      
        delta_ikb(2,j)  = delta_ikb(2,j) -flux1+flux3+flux10; % IkB
        delta_ikb(3,j)  = delta_ikb(3,j) -flux2+flux4+flux11; % IkBn 
        
        FLUX_TRIMER((i-1)*v.IKBS + j,:) = [flux1 flux2 flux3 flux4 flux5 flux6 flux10 flux11 flux8 flux9 flux7 flux7b];
    end
end

%--- Free IKB Reactions
%{
    The first section calculates the const and inducible txn for each IkB
    The second section uses those results and calculates the reaction fluxes
        Matrix math is used for this 2nd section to speed things up.
%}
% Calculate constitutive and inducible txn
ikbtxn = zeros(1,v.IKBS);
for i = 1:v.IKBS
    if ~v.IKB_KO(i) %KO
        continue;
    end
    
    if sum(v.pdi(7,:,i)) == 0 % not induced by any dimer
        ikbtxn(1,i) = v.pi(1,i);
        continue;
    end
    
    % Big ugly synthesis reaction
    %{
     Constitutive and NFkB-induced transcription
       const rate * ( 1 + sum( w*([d]/Kd)^hill)) / ...
                    ( 1 + sum( ([d]/Kd)^hill))  
     p50-Hom is assumed not to bind
    %}
     itxn = v.pi(1,i) * ...
      ((1+sum((v.pdi(7,IKB_DIMS{i},i)) .* (dimer(2,IKB_DIMS{i}) / v.pi(10,i)).^v.pi(9,i)))...
      / (1+sum((dimer(2,IKB_DIMS{i}) / v.pi(10,i)).^v.pi(9,i))));  

    if (v.PHASE == 1) || (v.pi(8,i) == 0) %skip if phase1 or no delay time
        ikbtxn(1,i) = itxn;   
    else  % push itxn to cache array and , if t>delay, calculate previous itxn value  
        if ~any(DELAY_IKB{i,2} == t,1) 
            DELAY_IKB{i,1}(DELAY_IKB{i,3}) = itxn; 
            DELAY_IKB{i,2}(DELAY_IKB{i,3}) = t; 
            DELAY_IKB{i,3} = DELAY_IKB{i,3} + 1;
        end
        % Calculate the concentration of delayed transcription 
        if ( t <= v.pi(8,i) )  % time <= delay so txn is at basal level
            ikbtxn(1,i) =DELAY_IKB{i,1}(1);
        else   % time > delay so interpolate for [mRNA] made at t - delay    
            %{
            start = find last index where t <= (t-delay)
            f_y = DELAY_IKB{i,1}(start);    (start value)
            f_x = DELAY_IKB{i,2}(start);    (start time)
            n_y = DELAY_IKB{i,1}(start+1);  (stop value)
            n_x = DELAY_IKB{i,2}(start+1);  (stop time)
            m = (n_y - f_y) / (n_x - f_x);  (slope)
            b = f_y;                        (start value = y-intercept)
            x = (t-v.pi(8,i)) - f_x;        (point on the line)
            ikbtxn(1,i) =  m*x + b;
            
            Old: ikbtxn(1,i) = interp1(DELAY_IKB{i,2}(1:DELAY_IKB{i,3}), ...
                 DELAY_IKB{i,1}(1:DELAY_IKB{i,3}),t-v.pi(8,i));
            %} 
            start = find(DELAY_IKB{i,2}(1:DELAY_IKB{i,3}-1) ...
                    <= (t-v.pi(8,i)),1,'last');
            
            ikbtxn(1,i)=((DELAY_IKB{i,1}(start+1) - DELAY_IKB{i,1}(start)) ...
                       /(DELAY_IKB{i,2}(start+1) - DELAY_IKB{i,2}(start))) ...
                       * ((t-v.pi(8,i)) - DELAY_IKB{i,2}(start)) ...
                       + DELAY_IKB{i,1}(start);  
        end
    end
end

% Calculate fluxes for all IkBs
fluxes3 = (v.pi(2,:) .* ikb(1,:));   % transcript deg
fluxes4 = (v.pi(3,:) .* ikb(1,:)); % translation
fluxes5 = (v.pi(4,:) .* ikb(2,:)); % cyt protein deg
fluxes6 = (v.pi(11,:).* ikb(2,:) * ikk_activity); % IKK2 deg
fluxes7 = (v.pi(5,:) .* ikb(3,:)); % nuc protein deg
fluxes8 = (v.pi(6,:) .* ikb(2,:)); % nuclear import
fluxes9 = (v.pi(7,:) .* ikb(3,:)); % nuclear export
%fluxes10 = (v.pi(15,:).* (ikb(2,:)*other(1))./(v.ps(6)+ikb(2,:))); %NIK deg
fluxes10 = (v.pi(15,:).* (ikb(2,:)*other(1))./(v.ps(6)*(1+monomer(2,4)/v.ps(3))+ikb(2,:))); %NIK deg

    delta_ikb(1,:) = delta_ikb(1,:) ... % IkB mRNA
                    + ikbtxn(1,:) ...            % mRNA synthesis
                    - fluxes3;   % transcript deg

    delta_ikb(2,:) = delta_ikb(2,:) ... % IkB free in cytoplasm
                    + fluxes4 ...% translation
                    - fluxes5 ...% cyt protein deg
                    - fluxes8 ...% nuclear import
                    + fluxes9 ...% nuclear export
                    - fluxes6... % IKK2 deg 
                     -fluxes10; % nik deg
                        % NIK deg :  rate * (IkBd*nik)/(Km + IkBd)
                        %- (v.pi(15,:).* (ikb(2,:)*other(1))); % nik deg

    delta_ikb(3,:) = delta_ikb(3,:) ... % IkB free in nucleus
                    - fluxes7 ...% nuc protein deg
                    + fluxes8 ...% nuclear import
                    - fluxes9;   % nuclear export
FLUX_IKB(:,:) = [v.pi(1,:); ikbtxn(1,:); fluxes3; fluxes4; fluxes5; fluxes6; fluxes7; fluxes8; fluxes9; fluxes10; ]';
 
%--- Free NFkB Monomers
% Calculate constitutive and inducible txn, NIK processing
   montxn = zeros(1,v.MONOMERS);
   FLUX_NIK_deg = zeros(1,v.MONOMERS); 
   for i = 1:v.MONOMERS
     if ~v.MONOMER_KO(i) || (v.pm(1,i) == 0) %KO or no const txn 
        continue;
     end
     
     %--- NIK-mediated processing of free monomer in the cytoplasm
     %sum of delta species = dimer(2,4)+dimer(3,5)
     if v.pm(10,i)
        %flux =  v.pm(9,i)* (monomer(2,i)*other(1))/(v.ps(3)+monomer(2,i));
        flux = v.pm(9,i)*monomer(2,i)*other(1) / (v.ps(3)*(1+sum(dimer(cyt,:))/v.ps(6)) + monomer(2,i));
        FLUX_NIK_deg(i) = flux;
        
        delta_monomer(2,i) = delta_monomer(2,i)- flux;
        delta_monomer(2,v.pm(10,i)) = delta_monomer(2,v.pm(10,i)) + flux;
     end
     
     %--- Synthesis
        % only if const txn rate > 0 and is induced by a dimer
     if (v.pm(1,i)==0) || (sum(v.pdm(1,:,i))== 0)  
         montxn(1,i) = v.pm(1,i);
         continue;
     end

     
     % modified txn rates
     %%%%wm=v.pdm(1,7,4)/50;    
      if i==4
           const2=v.pm(1,i).*(v.pdm(1,MON_DIMS{i},i)-1);
           itxn= v.pm(1,i) + sum(const2.*((dimer(2,MON_DIMS{i}) / v.pm(8,i)).^v.pm(7,i)))./(1+sum(dimer(2,MON_DIMS{i}) / v.pm(8,i)).^v.pm(7,i));
      else

     itxn = v.pm(1,i) + ...
     ((sum((v.pdm(1,MON_DIMS{i},i)) .* ((dimer(2,MON_DIMS{i})/ v.pm(8,i)).^v.pm(7,i))))...
     / ((1 +sum((dimer(2,MON_DIMS{i}) / v.pm(8,i)).^v.pm(7,i)))))./100000;
       end
     
     
     % no delay in Phase 1 or if delay== 0
     if (v.PHASE == 1) || (v.pm(6,i) == 0)
        montxn(1,i) = itxn; 
        
     else  % Phase 2 and delay_time > 0
        % Push time and  the sum of all mRNA synthesis onto cache array
        if ~any(DELAY_MON{i,2} == t,1)
            DELAY_MON{i,1}(DELAY_MON{i,3}) = itxn; 
            DELAY_MON{i,2}(DELAY_MON{i,3}) = t; 
            DELAY_MON{i,3} = DELAY_MON{i,3} + 1;
        end

        % Calculate the concentration of delayed transcription 
        %%%if ( wm<1  || t <= v.pm(6,i))  % time <= delay so txn is at basal level
        if (t <= v.pm(6,i))  % time <= delay so txn is at basal level
            montxn(1,i) = DELAY_MON{i,1}(1);
        else  % time > delay so interpolate for [mRNA] made at t - delay
            start = find(DELAY_MON{i,2}(1:DELAY_MON{i,3}-1) ...
                    <= (t-v.pm(6,i)),1,'last');
            
            montxn(1,i)=((DELAY_MON{i,1}(start+1) - DELAY_MON{i,1}(start)) ...
                       /(DELAY_MON{i,2}(start+1) - DELAY_MON{i,2}(start))) ...
                       * ((t-v.pm(6,i)) - DELAY_MON{i,2}(start)) ...
                       + DELAY_MON{i,1}(start);
        end
     end
     
   end
  
% Calculate fluxes for all monomers
fluxes2 = (v.pm(2,:) .* monomer(1,:));  % transcript deg
fluxes3 = (v.pm(3,:) .* monomer(1,:)); % translation
fluxes4 = (v.pm(4,:) .* monomer(2,:));   % protein deg cyt
fluxes5 = (v.pm(5,:) .* monomer(3,:));   % protein deg nuc

  delta_monomer(1,:) = delta_monomer(1,:) ... % monomer mRNA
                    + montxn(1,:) ...                % mRNA synthesis
                    - fluxes2;   % transcript deg
  
  delta_monomer(2,:) = delta_monomer(2,:) ... % monomer in cyt
                    + fluxes3 ...% translation
                    - fluxes4;   % protein deg cyt
  
  delta_monomer(3,:) = delta_monomer(3,:) ... % monomer in nuc
                    - fluxes5;   % protein deg nuc

FLUX_MON(:,:) = [v.pm(1,:); montxn(1,:); fluxes2; fluxes3; fluxes4; fluxes5; FLUX_NIK_deg]';
                
%--- Other reactions
  NIKDegTotal=max((7.7e-4*60./5)./12,v.ps(2)*nik_deg_mod);
  
  %delta_other(1) = delta_other(1) +v.ps(1)...          % NIK protein (cyt)
  %                  - (v.ps(2) * other(1) * nik_deg_mod);% NIK deg
  delta_other(1) = delta_other(1) +v.ps(1)...          % NIK protein (cyt)
                    - (NIKDegTotal * other(1));% NIK deg
  % hard-wired reaction: make sure you update this if you change the model
  %  p100 homodimerization in the cyt/nuc to form IkBd
  flux1 = v.ps(4)*monomer(2,4)*monomer(2,4); % p100 + p100 -> Ikbd
  flux2 = v.ps(5)*ikb(2,4); % IkBd -> p100 + p100
  
  flux3 = v.ps(4)*monomer(3,4)*monomer(3,4); % p100n + p100n -> Ikbdn
  flux4 = v.ps(5)*ikb(3,4); % IkBdn -> p100n + p100n
  
  delta_ikb(2,4) = delta_ikb(2,4) + flux1 - flux2;
  delta_monomer(2,4) = delta_monomer(2,4) - 2*flux1 + 2*flux2;
  
  delta_ikb(3,4) = delta_ikb(3,4) + flux3 - flux4;
  delta_monomer(3,4) = delta_monomer(3,4) - 2*flux3 + 2*flux4;
  
  FLUX_SPECIAL = [flux1 flux2 flux3 flux4]; %p100/Ikbd Cyt Dimerization, Cyt dissoc; Nuc dimerization, Nuc dissoc
                
%--- Display certain fluxes (debug mode only)
%if v.PHASE == 2 && t ==0
%disp(ikb(1,1));
%end

%--- Save concentrations for the next time step
delta  = [delta_ikb(:);delta_monomer(:);delta_dimer(:);delta_other(:)];

% save fluxes
%FLUX_CELL{size(FLUX_CELL,2)+1} = {FLUX_IKB, FLUX_MON, FLUX_DIMER, FLUX_TRIMER};

% Instead of adding a new row to FLUX_Vec each time, I add defaultSize(1000) at
% a time, then get rid of the extra 0's later.

if(size(FLUX_VEC,1) < timeStepCount)
    FLUX_VEC(timeStepCount + defaultSize,numel(FLUX_IKB) + numel(FLUX_MON) + numel(FLUX_DIMER) + numel(FLUX_TRIMER) + numel(FLUX_SPECIAL)) = 0;
    FLUX_T(timeStepCount + defaultSize) = 0;
    disp('extend')
end
FLUX_VEC(timeStepCount,:) = [reshape(FLUX_IKB',1,numel(FLUX_IKB)) , ...
    reshape(FLUX_MON',1,numel(FLUX_MON)) , ...
    reshape(FLUX_DIMER',1,numel(FLUX_DIMER)) , ...
    reshape(FLUX_TRIMER',1,numel(FLUX_TRIMER)), ...
    reshape(FLUX_SPECIAL',1,numel(FLUX_SPECIAL))];
FLUX_T(timeStepCount) = t;
 if any(~arrayfun(@isreal,FLUX_VEC(timeStepCount,:)))
     1;
 end
end