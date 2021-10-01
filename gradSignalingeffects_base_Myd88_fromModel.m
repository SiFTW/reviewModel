%trying to gradiently change the levels of NIK abundance and IKK activity over a gradient 

function v = gradSignalingeffects_base_Myd88(paramFile)

%IKKmulti = [100:-1:0]; 
%NIKab = [100:-1:0]; 

IKKmulti = [13.58,1]; 
%NIKab = [100,90,80,70,60,50,40,30,20,10,1]; 
NIKab = [1]; 

%species = {zeros(121,1)}; 

for i = 1:length(NIKab)   % so NIK = 1 nM 
    a = setDefaults('varyNIK'); 
    
    %species_mod = species; 
    %species_mod = {[zeros(120,1);NIKab(i)]}; 
    %a(size(a,1)+1,:) = {'SPECIES', species_mod};
    rate_mod={'NIK_multiplier',NIKab(i)};
    a(size(a,1)+1,:) = rate_mod;
    
    a = nfkb_run_62_passParams1_mod(a, paramFile);
    
    v.OUTPUT1(i,:) = a.OUTPUT1;
    v.A50(i,:) = sum(a.OUTPUT2(:,[32]),2);
    v.A52(i,:) = sum(a.OUTPUT2(:,92),2);
    v.A50basal(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,32),2);
    v.A52basal(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,92),2);
    v.B50(i,:) = sum(a.OUTPUT2(:,[42]),2);
    v.B52(i,:) = sum(a.OUTPUT2(:,52),2);
    v.B50basal(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,42),2);
    v.B52basal(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,52),2);
    v.C50(i,:) = sum(a.OUTPUT2(:,[72]),2);
    v.C52(i,:) = sum(a.OUTPUT2(:,82),2);
    v.C50basal(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,72),2);
    v.C52basal(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,82),2);
    v.C100(i,:) = sum(a.OUTPUT2(:,[62]),2);
    v.p5050(i,:) = sum(a.OUTPUT2(:,[102]),2);
    v.p5252(i,:) = sum(a.OUTPUT2(:,[112]),2);
    v.IKK(i,:) = a.IKK2;
    v.IKKbasal(i,:) = ones(1,1000)*a.BASAL_IKK2*a.IKK2multipli;
    v.NIK(i,:) = sum(a.OUTPUT2(:,[121]),2);
    v.NIKbasal(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,121),2);
    v.p100(i,:) = sum(a.OUTPUT2(:,[11 12 23 24 39 40 49 50 59 60 61:70 79 80 89 90 99 100]),2);
    v.p100basal(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,[11 12 23 24 39 40 49 50 59 60 61:70 79 80 89 90 99 100]),2);
    v.p52(i,:) = (sum(a.OUTPUT2(:,[26 27 51:60 81:90 91:100 111:120]),2));
    v.p52basal(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,[26 27 51:60 81:90 91:100 111:120]),2);
    v.IkBa(i,:) = (sum(a.OUTPUT2(:,[2 3 33 34 43 44 53 54 63 64 73 74 83 84 93 94]),2));
    v.IkBb(i,:) = (sum(a.OUTPUT2(:,[5 6 35 36 45 46 55 56 65 66 75 76 85 86 95 96]),2));
	v.Ikbe(i,:) = (sum(a.OUTPUT2(:,[8 9 37 38 47 48 57 58 67 68 77 78 87 88 97 98]),2));
	v.IkBd(i,:) = (sum(a.OUTPUT2(:,[11 12 39 40 49 50 59 60 69 70 79 80 89 90 99 100]),2));
    v.totalOutput(i,:,:) = a.OUTPUT2(:,:); 
    v.IKKCurve=a.IKK2;
end

for i = 1:length(IKKmulti)   % so IKK activity = 1 percent
    b = setDefaults('med'); 
    
    b(size(b,1)+1,:) = {'IKK2multipli', IKKmulti(i)}; 
    b(size(b,1)+1,:) = {'BASAL_IKK2', IKKmulti(i)/100}; 
    
    b = nfkb_run_62_passParams1_mod(b, paramFile);
    
    v.OUTPUT1_I(i,:) = b.OUTPUT1;
    v.A50_I(i,:) = sum(b.OUTPUT2(:,[32]),2);
    v.A52_I(i,:) = sum(b.OUTPUT2(:,92),2);
    v.A50basal_I(i,:) = ones(1,1000)*sum(b.OUTPUT2(1,32),2);
    v.A52basal_I(i,:) = ones(1,1000)*sum(b.OUTPUT2(1,92),2);
    v.B50_I(i,:) = sum(b.OUTPUT2(:,[42]),2);
    v.B52_I(i,:) = sum(b.OUTPUT2(:,52),2);
    v.B50basal_I(i,:) = ones(1,1000)*sum(b.OUTPUT2(1,42),2);
    v.B52basal_I(i,:) = ones(1,1000)*sum(b.OUTPUT2(1,52),2);
    v.C50_I(i,:) = sum(b.OUTPUT2(:,[72]),2);
    v.C52_I(i,:) = sum(b.OUTPUT2(:,82),2);
    v.C50basal_I(i,:) = ones(1,1000)*sum(b.OUTPUT2(1,72),2);
    v.C52basal_I(i,:) = ones(1,1000)*sum(b.OUTPUT2(1,82),2);
    v.C100_I(i,:) = sum(b.OUTPUT2(:,[62]),2);
    v.p5050_I(i,:) = sum(b.OUTPUT2(:,[102]),2);
    v.p5252_I(i,:) = sum(b.OUTPUT2(:,[112]),2);
    v.IKK_I(i,:) = b.IKK2;
    v.IKKbasal_I(i,:) = ones(1,1000)*b.BASAL_IKK2;
    v.NIK_I(i,:) = sum(b.OUTPUT2(:,[121]),2);
    v.NIKbasal_I(i,:) = ones(1,1000)*sum(b.OUTPUT2(1,121),2);
    v.p100_I(i,:) = sum(b.OUTPUT2(:,[11 12 23 24 39 40 49 50 59 60 61:70 79 80 89 90 99 100]),2);
    v.p100basal_I(i,:) = ones(1,1000)*sum(b.OUTPUT2(1,[11 12 23 24 39 40 49 50 59 60 61:70 79 80 89 90 99 100]),2);
    v.p52_I(i,:) = (sum(b.OUTPUT2(:,[26 27 51:60 81:90 91:100 111:120]),2));
    v.p52basal_I(i,:) = ones(1,1000)*sum(b.OUTPUT2(1,[26 27 51:60 81:90 91:100 111:120]),2);
    v.nfkb2mRNA_I(i,:) = sum(b.OUTPUT2(:,[22]),2);
    v.IkBa_I(i,:) = (sum(b.OUTPUT2(:,[2 3 33 34 43 44 53 54 63 64 73 74 83 84 93 94]),2));
    v.IkBb_I(i,:) = (sum(b.OUTPUT2(:,[5 6 35 36 45 46 55 56 65 66 75 76 85 86 95 96]),2));
	v.Ikbe_I(i,:) = (sum(b.OUTPUT2(:,[8 9 37 38 47 48 57 58 67 68 77 78 87 88 97 98]),2));
	v.IkBd_I(i,:) = (sum(b.OUTPUT2(:,[11 12 39 40 49 50 59 60 69 70 79 80 89 90 99 100]),2));
    v.totalOutput_I(i,:,:) = b.OUTPUT2(:,:);
end


for i = length(NIKab)   % so NIK = 1 nM 
    a = setDefaults('varyNIK'); 
    a(size(a,1)+1,:) = {'NIK_multiplier',NIKab(i)};
    a(size(a,1)+1,:) = {'IKK2_TYPE', 'IGM_mod'};
    
    a = nfkb_run_62_passParams1_mod(a, paramFile);
    
    %Total IkBd alpha instead of just nucleus
    v.OUTPUT1_A(i,:) = a.OUTPUT1;
    v.A50_A(i,:) = sum(a.OUTPUT2(:,[31 32]),2);
    v.A52_A(i,:) = sum(a.OUTPUT2(:,[91 92]),2);
    v.A50basal_A(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,32),2);
    v.A52basal_A(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,92),2);
    v.B50_A(i,:) = sum(a.OUTPUT2(:,[41 42]),2);
    v.B52_A(i,:) = sum(a.OUTPUT2(:,[51 52]),2);
    v.B50basal_A(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,42),2);
    v.B52basal_A(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,52),2);
    v.C50_A(i,:) = sum(a.OUTPUT2(:,[71 72]),2);
    v.C52_A(i,:) = sum(a.OUTPUT2(:,[81 82]),2);
    v.C50basal_A(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,72),2);
    v.C52basal_A(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,82),2);
    v.C100_A(i,:) = sum(a.OUTPUT2(:,[61 62]),2);
    v.p5050_A(i,:) = sum(a.OUTPUT2(:,[101 102]),2);
    v.p5252_A(i,:) = sum(a.OUTPUT2(:,[111 112]),2);
    v.IKK_A(i,:) = a.IKK2;
    v.IKKbasal_A(i,:) = ones(1,1000)*a.BASAL_IKK2*a.IKK2multipli;
    v.NIK_A(i,:) = sum(a.OUTPUT2(:,[121]),2);
    v.NIKbasal_A(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,121),2);
    v.p100_A(i,:) = sum(a.OUTPUT2(:,[11 12 23 24 39 40 49 50 59 60 61:70 79 80 89 90 99 100]),2);
    v.p100basal_A(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,[11 12 23 24 39 40 49 50 59 60 61:70 79 80 89 90 99 100]),2);
    v.p52_A(i,:) = (sum(a.OUTPUT2(:,[26 27 51:60 81:90 91:100 111:120]),2));
    v.p52basal_A(i,:) = ones(1,1000)*sum(a.OUTPUT2(1,[26 27 51:60 81:90 91:100 111:120]),2);
    v.IkBa_A(i,:) = (sum(a.OUTPUT2(:,[2 3 33 34 43 44 53 54 63 64 73 74 83 84 93 94]),2));
    v.IkBb_A(i,:) = (sum(a.OUTPUT2(:,[5 6 35 36 45 46 55 56 65 66 75 76 85 86 95 96]),2));
	v.Ikbe_A(i,:) = (sum(a.OUTPUT2(:,[8 9 37 38 47 48 57 58 67 68 77 78 87 88 97 98]),2));
	v.IkBd_A(i,:) = (sum(a.OUTPUT2(:,[11 12 39 40 49 50 59 60 69 70 79 80 89 90 99 100]),2));
    v.totalOutput_A(i,:,:) = a.OUTPUT2(:,:); 
end





end
