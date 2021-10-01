function equilstate = gradSignalingeffects_equil(paramFile)

    a = setDefaults('varyNIK'); 
    %species = {zeros(121,1)};
    %species_mod = species; 
    %species_mod = {[zeros(120,1); 1]}; 
    rate_mod={'NIK_multiplier',1};
    a(size(a,1)+1,:) = rate_mod;
    %a(size(a,1)+1,:) = {'SPECIES', species_mod};
    
%     a = setDefaults('med'); 
%     a(size(a,1)+1,:) = {'IKK2multipli', IKKmulti(i)}; 
%     a(size(a,1)+1,:) = {'BASAL_IKK2', IKKmulti(i)/100};     

    a = nfkb_run_62_passParams1_mod(a, paramFile);
    equilstate = a.OUTPUT1;
end