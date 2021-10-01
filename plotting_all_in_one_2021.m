%%v.totalOutput_I = IKK SCAN
%%v.totalOutput = NIK SCAN
%%
%setup
clear all;
%load('lastRun.mat');

%uncomment to run model
 speciesNames={ ...  % Component # 1-12 (ikb) , 13-30 (NFkB mon) , 31-40 (A50), 
                               %41-50 (B50), 51-60 (B52), 61-70 (C100), 71-80 (C50), 81-90 (C52), 91-100 (A52), 
                               % 101-110 (5050), 111-120 (5252), 121 (NIK)
            'tikba','ikba', 'ikban','tikbb','ikbb','ikbbn','tikbe','ikbe','ikben','tikbd','ikbd','ikbdn',...
            'tRelA', 'RelA','RelAn','tp50','p50','p50n', 'tRelB', 'RelB','RelBn','tp100','p1002','p100n','tp52','p52','p52n', 'tcRel', 'cRel', 'cReln'...
            'A50','A50n','ikbaA50','ikbaA50n','ikbbA50','ikbbA50n','ikbeA50','ikbeA50n','ikbdA50','ikbdA50n',...
            'B50','B50n','ikbaB50','ikbaB50n','ikbbB50','ikbbB50n','ikbeB50','ikbeB50n','ikbdB50','ikbdB50n',...
            'B52','B52n','ikbaB52','ikbaB52n','ikbbB52','ikbbB52n','ikbeB52','ikbeB52n','ikbdB52','ikbdB52n',...
            'C100','C100n','ikbaC100','ikbaC100n','ikbbC100','ikbbC100n','ikbeC100','ikbeC100n','ikbdC100','ikbdC100n',...
            'C50','C50n','ikbaC50','ikbaC50n','ikbbC50','ikbbC50n','ikbeC50','ikbeC50n','ikbdC50','ikbdC50n',...
            'C52','C52n','ikbaC52','ikbaC52n','ikbbC52','ikbbC52n','ikbeC52','ikbeC52n','ikbdC52','ikbdC52n',...
            'A52','A52n','ikbaA52','ikbaA52n','ikbbA52','ikbbA52n','ikbeA52','ikbeA52n','ikbdA52','ikbdA52n',...
            '5050', '5050n', 'ikba5050', 'ikba5050n', 'ikbb5050', 'ikbb5050n', 'ikbe5050', 'ikbe5050n', 'ikbd5050', 'ikbd5050n',...
            '5252', '5252n', 'ikba5252', 'ikba5252n', 'ikbb5252', 'ikbb5252n', 'ikbe5252', 'ikbe5252n', 'ikbd5252', 'ikbd5252n',...
            'NIK' ...
            }';
celltypes = {'mef', 'bcell', 'tcell', 'dc'};



for i = 1:numel(celltypes)
    filename = strcat('passParams030614_', celltypes{i});
    data(i) = gradSignalingeffects_base(filename);
    estate(i,:) = gradSignalingeffects_equil(filename);
end
%uncomment to run model

%%
%plotting the NIK curves
figure;
subplot(1,2,1);
plot(data(1).IKKCurve,'linewidth',3);
box off;
set(gca,'xtick',[0:60:9*60]);
set(gca,'xticklabels',[0:9]);
xlim([0,9*60]);
title('ikk curve');

subplot(1,2,2);
plot(data(1).NIK_I(end,:),'linewidth',3);
box off;
set(gca,'xtick',[0:60*8:8*40*60]);
set(gca,'xticklabels',[0:8:40]);
xlim([0,40*60]);
set(gcf,'color','w');
title('NIK curve');

%dimers = {'A50','B50','B52','C50','C52','A52','5050','5252'};
dimers = {'A50','B50','B52','C50'};
dimer_ind = [31,41,51,71,81,91,101,111];

for j = 1:numel(dimer_ind)
    sind = dimer_ind(j);
    for i = 1:numel(celltypes)
        res(j,i,:) = [estate(i,sind+0)+estate(i,sind+1), estate(i,sind+2)+estate(i,sind+3), ...
            estate(i,sind+4)+estate(i,sind+5), estate(i,sind+6)+estate(i,sind+7), ...
            estate(i,sind+8)+estate(i,sind+9)];
    end
end

%%
%plotting which inhibitors are bound to each NFkB dimer
figure;
imgi= 1;
cellcat = categorical(celltypes);
%unbound, alpha, beta, eps, delta
inhibitorColors=[0.25,0.25,1;1,1,0;0.75,0.75,0;1,0.75,0;1,0.5,0];

for idx = 1:numel(dimers)
    subplot(1,4,imgi);
    imgi = imgi + 1;
    H = bar(squeeze(res(idx,:,:)),'stacked');
    H(1).FaceColor = 'flat';
    H(1).CData = inhibitorColors(1,:);
    H(2).FaceColor = 'flat';
    H(2).CData = inhibitorColors(2,:);
    H(3).FaceColor = 'flat';
    H(3).CData = inhibitorColors(3,:);
    H(4).FaceColor = 'flat';
    H(4).CData = inhibitorColors(4,:);
    H(5).FaceColor = 'flat';
    H(5).CData = inhibitorColors(5,:);    
    set(gca,'xticklabel',celltypes);
    title(dimers{idx}, 'FontSize', 15);
    hLegend = legend(H, {'unbound','ikba','ikbb','ikbe','ikbd'});
    set(hLegend, 'Location','Best','position',[0.01 0.5 0.07 0.1]);
    ylim([0,450]);
    box off;
end
set(gcf,'color','w');

%%
%lets plot monomer abundances
proteinsToPlot={'A','50','ikba','ikbb','ikbe','C','B'};
indexesToPlot={};
for i=1:length(proteinsToPlot)
    theseIndexes=find(~cellfun(@isempty,strfind(speciesNames,proteinsToPlot{i})));    
    theseIndexes=theseIndexes(2:end);
    indexesToPlot{i}=theseIndexes;
end

figure;
title('monomers')
for cell = 1:numel(celltypes)
    for idx = 1:numel(proteinsToPlot)        
        subplot(1,numel(proteinsToPlot),idx);
        hold on;
        title(proteinsToPlot{idx});
        indexesThisBar=indexesToPlot{idx};
        thisBarVal=sum(estate(cell,indexesThisBar));
        H = bar(cell,thisBarVal,'FaceColor',[0.5,0.5,0.5]);
        set(gca,'xticklabel',celltypes);
        set(gca,'xTick',[1:length(celltypes)]);
        box off;
    end
end
set(gcf,'color','w');

%%
%let's plot the activity of the dimers in each cell type in response to IKK
%these are line graphs
names = {'A50n','C50n','B50n','B52n'};
indexesToPlot={};
for i=1:length(names)
    theseIndexes=find(strcmp(speciesNames,names{i}));    
    indexesToPlot{i}=theseIndexes;
end
figure;
imgi= 1;
title('dimers');
basalDimers=zeros(numel(celltypes),numel(indexesToPlot));
colorArray={[0,0.75,0],[1,0,0],[0.5,0,0.5],[0,0,1]};
for cell = 1:numel(celltypes)
    subplot(1,4,cell);
    for idx = 1:numel(names)
        toPlot=(data(cell).totalOutput_A(end,:,indexesToPlot{idx}));
        basalDimers(cell,idx)=toPlot(1);
        plot((1:2881)/(2881/48), toPlot, 'LineWidth', 2,'color',colorArray{idx});
        hold on;
        title(strcat(celltypes{cell}, ', ', names{idx},',fixed NIK=1nM'), 'FontSize', 10);
        ylim([0,250]);
        xlabel('time (hours)');
        ylabel(names{idx}, 'FontSize', 10);
        xlim([0,9]);
    end
end
set(gcf,'color','w');


%%
% let's plot dimers in each cell type in response to NIK
figure;
imgi= 1;
title('dimers');
for cell = 1:numel(celltypes)
    subplot(1,4,cell);
    for idx = 1:numel(names)
        toPlot=(data(cell).totalOutput_I(end,:,indexesToPlot{idx}));
        hold on;
        imgi = imgi + 1;
        plot((1:2881)/(2881/48), toPlot, 'LineWidth', 2,'color',colorArray{idx});
        title(strcat(celltypes{cell}, ', ', names{idx},',fixed IKK=1%'), 'FontSize', 10);
        ylim([0,200]);
        xlabel('time (hours)');
        ylabel(names{idx}, 'FontSize', 10);
        xlim([0,48]);
    end
end
set(gcf,'color','w');

%%
%NFkB dimer composition stacked bar graph
figure;
h= bar(basalDimers,'stacked');
colorsForUnbound={[21, 67, 96]./255;[26, 188, 156]./255;[93, 173, 226]./255;[165, 105, 189]./255};
set(h,{'FaceColor'},colorsForUnbound);
set(gca,'xticklabels',celltypes);
hLegend = legend(names);
set(hLegend, 'Location','Best','position',[0.01 0.5 0.07 0.1]);
set(gcf,'color','w');
title('unbound dimer composition');
box off;
%%
%heatmaps of NFkB responses with different levels of NIK
figure;
imgi=1;
for cell = 1:numel(celltypes)
    for idx = 1:numel(names)
        
        toPlot=[data(cell).totalOutput(:,:,indexesToPlot{idx})];
        toPlot=[toPlot,toPlot(:,end)];
        toPlot=[toPlot;toPlot(end,:)];
        subplot(4,4,imgi);
        imgi = imgi + 1;
        toPlot=flipud(toPlot);
        surf(toPlot);
        view(2);
        shading(gca,'interp');
        title(strcat(celltypes{cell}, ', ', names{idx}), 'FontSize', 10);
        zlim([0,550]);
        xlabel('time (mins)');
        zlabel(names{idx}, 'FontSize', 10);
        ylabel('NIK');        
        xlim([0,60*8]);
        caxis([0,320]);
        colormap('jet');
        box off;
        set(gca,'TickDir','out')
        ylim([1,size(toPlot,1)]);
    end
end
set(gcf,'color','w');
%%
%plot A50 highlight
cell=1;
idx=1;
figure;
 toPlot=[data(cell).totalOutput(:,:,indexesToPlot{idx})];
toPlot=[toPlot,toPlot(:,end)];
toPlot=[toPlot;toPlot(end,:)];
toPlot=flipud(toPlot);
surf(toPlot,'linestyle','none','FaceColor','interp');
hold on;
for i=1:size(toPlot,1)
    plot3([1:length(toPlot(i,:))]',repmat(i,1,2882)',toPlot(i,:)','k','linewidth',0.5);
end
%view(2);
view(26,19);
title(strcat(celltypes{cell}, ', ', names{idx}), 'FontSize', 10);
zlim([0,550]);
xlabel('time (mins)');
zlabel(names{idx}, 'FontSize', 10);
ylabel('NIK');
xlim([0,60*8]);
caxis([0,320]);
colormap('jet');
box off;
set(gcf,'color','w');

%%
%plot dimer highlights with increasing basal IKK
figure;
imgi=1;
for cell = 1:numel(celltypes)
    for idx = 1:numel(names)
        %nums = eval(strcat('data(cell).', names{idx}, '(100,:)'));
        
        toPlot=[data(cell).totalOutput_I(:,:,indexesToPlot{idx})];
        toPlot=[toPlot,toPlot(:,end)];
        toPlot=[toPlot;toPlot(end,:)];
        %totalToPlot=sum(toPlot,2);
        subplot(4,4,imgi);
        imgi = imgi + 1;
        %surf(toPlot,'linestyle','none');
        toPlot=flipud(toPlot);
        surf(toPlot);
        view(2);
        %view(26,19);
        shading(gca,'interp');
        title(strcat(celltypes{cell}, ', ', names{idx}), 'FontSize', 10);
        zlim([0,550]);
        xlabel('time (mins)');
        zlabel(names{idx}, 'FontSize', 10);
        ylabel('IKK');        
        xlim([0,60*40]);
        
        caxis([0,60]);

        colormap('jet');
        box off;
        set(gca,'TickDir','out')
        ylim([1,size(toPlot,1)]);
    end
end
set(gcf,'color','w');

names = {'p52','RelB','NIK'};
indexesToPlot={};
for i=1:length(names)
    theseIndexes=find(strcmp(speciesNames,names{i}));    
    indexesToPlot{i}=theseIndexes;
end

%%
%plot non-canonical players with increasing basal IKK
figure;
imgi=1;
for cell = 1:numel(celltypes)
    for idx = 1:numel(names)
        
        toPlot=[data(cell).totalOutput_I(:,:,indexesToPlot{idx})];
        toPlot=[toPlot,toPlot(:,end)];
        toPlot=[toPlot;toPlot(end,:)];
        toPlot=flipud(toPlot);
        subplot(4,3,imgi);
        imgi = imgi + 1;
        toPlot=flipud(toPlot);
        surf(toPlot);
        view(2);
        shading(gca,'interp');
        title(strcat(celltypes{cell}, ', ', names{idx}), 'FontSize', 10);
        zlim([0,550]);
        xlabel('time (mins)');
        zlabel(names{idx}, 'FontSize', 10);
        ylabel('IKK');        
        xlim([0,60*40]);
        
        caxis([0,60]);

        colormap('jet');
        box off;
        set(gca,'TickDir','out')
        ylim([1,size(toPlot,1)]);
    end
end
set(gcf,'color','w');

%%
%Rel B p52 highlight
names = {'A50n','C50n','B50n','B52n'};
indexesToPlot={};
for i=1:length(names)
    theseIndexes=find(strcmp(speciesNames,names{i}));    
    indexesToPlot{i}=theseIndexes;
end

cell=1;
idx=4;
figure;
 toPlot=[data(cell).totalOutput_I(:,:,indexesToPlot{idx})];
toPlot=[toPlot,toPlot(:,end)];
toPlot=[toPlot;toPlot(end,:)];
toPlot=flipud(toPlot);
surf(toPlot,'linestyle','none','FaceColor','interp');
hold on;
for i=1:size(toPlot,1)
    plot3([1:length(toPlot(i,:))]',repmat(i,1,2882)',toPlot(i,:)','k','linewidth',0.5);
end
view(26,19);
title(strcat(celltypes{cell}, ', ', names{idx}), 'FontSize', 10);
zlim([0,150]);
xlabel('time (mins)');
zlabel(names{idx}, 'FontSize', 10);
ylabel('IKK');
xlim([0,60*40]);
caxis([0,60]);
colormap('jet');
box off;
set(gcf,'color','w');




%%
%run myd88 mutant sims
%myd88 mutation = high basal NFkB (13.5nM) then look at NIK responsive dimers.
%celltypes = {'bcell'};
filename = strcat('passParams030614_bcell');
names={'A50_I','C50_I','B50_I','B52_I','IkBd_I'};
myd88Data = gradSignalingeffects_base_Myd88_fromModel(filename);

%%
%plot Myd88 mutation data
normalize=0;
plottedLines={4,numel(names)};
for idx = 1:numel(names)
    figure;    
    %wt
    toPlot=squeeze(myd88Data.(names{idx})(2,:));
    %toPlot=sum(dataToPlot,2);
    basalDimers(1,idx)=toPlot(1);
    %totalToPlot=sum(toPlot,2);
    if(normalize) toPlot=toPlot./toPlot(1); end
    plot((1:2881)/(2881/48), toPlot, 'LineWidth', 2,'color',colorArray{1});
    plottedLines{1,idx}=toPlot;
    hold on;
    

    %high basal IKK
    toPlot=squeeze(myd88Data.(names{idx})(1,:));
    %toPlot=sum(dataToPlot,2);
    if(normalize) toPlot=toPlot./toPlot(1); end
    plot((1:2881)/(2881/48), toPlot, 'LineWidth', 2,'color',colorArray{3});
    plottedLines{3,idx}=toPlot;

 

    title(strcat('bcell', ', ', names{idx},',varied Basal IKK'), 'FontSize', 10);
    %ylim([0,250]);
    ylim('auto');
    xlabel('time (hours)');
    ylabel(names{idx}, 'FontSize', 10);
    %xlim([0,9]);
    xlim([0,24]);
    legend({'WT','Myd88'});
    set(gcf,'color','w');
    box off;
    set(gca,'TickDir','out')
end

%%
%plot the input curves for Myd88 mutant
figure;plot((1:2881)/(2881/48),myd88Data.NIK_I(:,:)','LineWidth', 2,'color',[0,0,1]);
set(gcf,'color','w');
box off;
xlim([0,24]);
ylim([0,50]);
title('NIK');
figure;
plot((1:2881)/(2881/48),myd88Data.IKK_I(2,:)','LineWidth', 2,'color',colorArray{1});
hold on;
plot((1:2881)/(2881/48),myd88Data.IKK_I(1,:)','LineWidth', 2,'color',colorArray{3});
set(gcf,'color','w');
box off;
xlim([0,24]);
ylim([0,0.75]);
title('IKK');

%% IkBalpha knockout simulation . 
% let's check canonical stimuli with a couple of different NIK levels.
filename = strcat('passParams030614_bcell');
names={'A50','C50','B50','B52'};
IkBaData = gradSignalingeffects_base_Ikba(filename);
IkBaDataWT=gradSignalingeffects_base_IkbaWT(filename);
%% plt IkBa KO data
normalize=0;
plottedLines={4,numel(names)};
for idx = 1:numel(names)
    figure;
    %nums = eval(strcat('data(cell).', names{idx}, '(100,:)'));
    %this works uncomment to restore last figures
    %toPlot=(data(cell).totalOutput_A(100,:,indexesToPlot{idx}));
    
    %wt
    toPlot=squeeze(IkBaData.(names{idx})(1,:));
    %toPlot=sum(dataToPlot,2);
    basalDimers(1,idx)=toPlot(1);
    %totalToPlot=sum(toPlot,2);
    if(normalize) toPlot=toPlot./toPlot(1); end
    plot((1:2881)/(2881/48), toPlot, '-','LineWidth', 3,'color',colorArray{1});
    hold on;
    plot(([-500,0])./(2881/48), [toPlot(1),toPlot(1)], '-','LineWidth', 3,'color',colorArray{1});
    plottedLines{1,idx}=toPlot;
    hold on;
%     
%     %med basal NIK
%     toPlot=squeeze(IkBaData.(names{idx})(2,:));
%     %toPlot=sum(dataToPlot,2);
%     if(normalize) toPlot=toPlot./toPlot(1); end
%     plot((1:2881)/(2881/48), toPlot, '-','LineWidth', 3,'color',colorArray{2});
%     hold on;
%     plot(([-500,0])./(2881/48), [toPlot(1),toPlot(1)], '-','LineWidth', 3,'color',colorArray{2});
%     plottedLines{3,idx}=toPlot;    
% 
%     %high basal NIK
%     toPlot=squeeze(IkBaData.(names{idx})(3,:));
%     %toPlot=sum(dataToPlot,2);
%     if(normalize) toPlot=toPlot./toPlot(1); end
%     plot((1:2881)/(2881/48), toPlot, '-','LineWidth', 3,'color',colorArray{3});
%     hold on;
%     plot(([-500,0])./(2881/48), [toPlot(1),toPlot(1)], '-','LineWidth', 3,'color',colorArray{3});
%     plottedLines{3,idx}=toPlot;

 

    title(strcat('IkBa KObcell', ', ', names{idx},',fixed NIK'), 'FontSize', 10);
    ylim([0,50]);
    %ylim('auto');
    xlabel('time (hours)');
    ylabel(names{idx}, 'FontSize', 10);
    %xlim([0,9]);
    xlim([0,8]);
    %legend({'WT','med NIK','high NIK'});
    set(gcf,'color','w');
    box off;
    set(gca,'TickDir','out')
    
    % plot WT data
    %figure;
    %nums = eval(strcat('data(cell).', names{idx}, '(100,:)'));
    %this works uncomment to restore last figures
    %toPlot=(data(cell).totalOutput_A(100,:,indexesToPlot{idx}));
    
    %wt
    toPlot=squeeze(IkBaDataWT.(names{idx})(1,:));
    %toPlot=sum(dataToPlot,2);
    basalDimers(1,idx)=toPlot(1);
    %totalToPlot=sum(toPlot,2);
    if(normalize) toPlot=toPlot./toPlot(1); end
    plot((1:2881)/(2881/48), toPlot, '--','LineWidth', 2,'color',[colorArray{1},0.75]);
    hold on;
    plot(([-500,0])./(2881/48), [toPlot(1),toPlot(1)], '--','LineWidth', 2,'color',[colorArray{1},0.5]);
    plottedLines{1,idx}=toPlot;
    hold on;
    
%     %med basal NIK
%     toPlot=squeeze(IkBaDataWT.(names{idx})(2,:));
%     %toPlot=sum(dataToPlot,2);
%     if(normalize) toPlot=toPlot./toPlot(1); end
%     plot((1:2881)/(2881/48), toPlot, '--','LineWidth', 2,'color',[colorArray{2},0.25]);
%     hold on;
%     plot(([-500,0])./(2881/48), [toPlot(1),toPlot(1)], '--','LineWidth', 2,'color',[colorArray{2},0.25]);
%     plottedLines{3,idx}=toPlot;    
% 
%     %high basal NIK
%     toPlot=squeeze(IkBaDataWT.(names{idx})(3,:));
%     hold on;
%     plot(([-500,0])./(2881/48), [toPlot(1),toPlot(1)], '--','LineWidth', 2,'color',[colorArray{3},0.25]);
%     %toPlot=sum(dataToPlot,2);
%     if(normalize) toPlot=toPlot./toPlot(1); end
%     plot((1:2881)/(2881/48), toPlot, '--','LineWidth', 2,'color',[colorArray{3},0.25]);
%     plottedLines{3,idx}=toPlot;

 

    title(strcat('WT KObcell', ', ', names{idx},',fixed NIK'), 'FontSize', 10);
    %ylim([0,250]);
    %ylim('auto');
    xlabel('time (hours)');
    ylabel(names{idx}, 'FontSize', 10);
    %xlim([0,9]);
    xlim([-4,24]);
    %legend({'KO low NIK','KO med NIK','KO high NIK','WT low NIK','WT med NIK','WT high NIK'});
    %legend({'KO','WT'});
    set(gcf,'color','w');
    box off;
    set(gca,'TickDir','out')
end

%% plot the input curves for IkBa mutant
figure;plot((1:2881)/(2881/48),IkBaData.IKK(:,:)','LineWidth', 2,'color',[0,0,1]);
set(gcf,'color','w');
box off;
xlim([0,24]);
%ylim([0,50]);
title('IKK');
figure;
plot((1:2881)/(2881/48),IkBaData.NIK(1,:)','LineWidth', 2,'color',colorArray{1});
hold on;
% plot((1:2881)/(2881/48),IkBaData.NIK(2,:)','LineWidth', 2,'color',colorArray{2});
% plot((1:2881)/(2881/48),IkBaData.NIK(3,:)','LineWidth', 2,'color',colorArray{3});
set(gcf,'color','w');
box off;
xlim([0,24]);
%ylim([0,0.75]);
title('NIK');

%% nfkb2 constitutive processing mutation

%% IkBalpha knockout simulation . 
% let's check canonical stimuli with a couple of different NIK levels.
filename = strcat('passParams030614_dc_nfkb2');
nfkb2Data = gradSignalingeffects_base_IkbaWT(filename);
filename = strcat('passParams030614_dc');
nfkb2DataWT = gradSignalingeffects_base_IkbaWT(filename);
%% plt IkBa KO data
names={'A50','C50','B50','B52','p5252'};
normalize=0;
plottedLines={4,numel(names)};
for idx = 1:numel(names)
    figure;
    %nums = eval(strcat('data(cell).', names{idx}, '(100,:)'));
    %this works uncomment to restore last figures
    %toPlot=(data(cell).totalOutput_A(100,:,indexesToPlot{idx}));
    
    %wt
    toPlot=squeeze(nfkb2Data.(names{idx})(1,:));
    %toPlot=sum(dataToPlot,2);
    basalDimers(1,idx)=toPlot(1);
    %totalToPlot=sum(toPlot,2);
    if(normalize) toPlot=toPlot./toPlot(1); end
    plot((1:2881)/(2881/48), toPlot, '-','LineWidth', 3,'color',colorArray{4});
    hold on;
    plot(([-500,0])./(2881/48), [toPlot(1),toPlot(1)], '-','LineWidth', 3,'color',colorArray{4});
    plottedLines{1,idx}=toPlot;
    hold on;
%     
%     %med basal NIK
%     toPlot=squeeze(IkBaData.(names{idx})(2,:));
%     %toPlot=sum(dataToPlot,2);
%     if(normalize) toPlot=toPlot./toPlot(1); end
%     plot((1:2881)/(2881/48), toPlot, '-','LineWidth', 3,'color',colorArray{2});
%     hold on;
%     plot(([-500,0])./(2881/48), [toPlot(1),toPlot(1)], '-','LineWidth', 3,'color',colorArray{2});
%     plottedLines{3,idx}=toPlot;    
% 
%     %high basal NIK
%     toPlot=squeeze(IkBaData.(names{idx})(3,:));
%     %toPlot=sum(dataToPlot,2);
%     if(normalize) toPlot=toPlot./toPlot(1); end
%     plot((1:2881)/(2881/48), toPlot, '-','LineWidth', 3,'color',colorArray{3});
%     hold on;
%     plot(([-500,0])./(2881/48), [toPlot(1),toPlot(1)], '-','LineWidth', 3,'color',colorArray{3});
%     plottedLines{3,idx}=toPlot;

 

    title(strcat('nfkb2 mut dc', ', ', names{idx},',fixed NIK'), 'FontSize', 10);
    %ylim([0,250]);
    ylim('auto');
    xlabel('time (hours)');
    ylabel(names{idx}, 'FontSize', 10);
    %xlim([0,9]);
    xlim([0,8]);
    legend({'WT','med NIK','high NIK'});
    set(gcf,'color','w');
    box off;
    set(gca,'TickDir','out')
    
    % plot WT data
    %figure;
    %nums = eval(strcat('data(cell).', names{idx}, '(100,:)'));
    %this works uncomment to restore last figures
    %toPlot=(data(cell).totalOutput_A(100,:,indexesToPlot{idx}));
    
    %wt
    toPlot=squeeze(nfkb2DataWT.(names{idx})(1,:));
    %toPlot=sum(dataToPlot,2);
    basalDimers(1,idx)=toPlot(1);
    %totalToPlot=sum(toPlot,2);
    if(normalize) toPlot=toPlot./toPlot(1); end
    plot((1:2881)/(2881/48), toPlot, '--','LineWidth', 2,'color',[colorArray{4},0.5]);
    hold on;
    plot(([-500,0])./(2881/48), [toPlot(1),toPlot(1)], '--','LineWidth', 2,'color',[colorArray{4},0.5]);
    plottedLines{1,idx}=toPlot;
    hold on;
    
%     %med basal NIK
%     toPlot=squeeze(IkBaDataWT.(names{idx})(2,:));
%     %toPlot=sum(dataToPlot,2);
%     if(normalize) toPlot=toPlot./toPlot(1); end
%     plot((1:2881)/(2881/48), toPlot, '--','LineWidth', 2,'color',[colorArray{2},0.25]);
%     hold on;
%     plot(([-500,0])./(2881/48), [toPlot(1),toPlot(1)], '--','LineWidth', 2,'color',[colorArray{2},0.25]);
%     plottedLines{3,idx}=toPlot;    
% 
%     %high basal NIK
%     toPlot=squeeze(IkBaDataWT.(names{idx})(3,:));
%     hold on;
%     plot(([-500,0])./(2881/48), [toPlot(1),toPlot(1)], '--','LineWidth', 2,'color',[colorArray{3},0.25]);
%     %toPlot=sum(dataToPlot,2);
%     if(normalize) toPlot=toPlot./toPlot(1); end
%     plot((1:2881)/(2881/48), toPlot, '--','LineWidth', 2,'color',[colorArray{3},0.25]);
%     plottedLines{3,idx}=toPlot;

 

    title(strcat('WT dc', ', ', names{idx},',fixed NIK'), 'FontSize', 10);
    ylim([0,250]);
    %ylim('auto');
    xlabel('time (hours)');
    ylabel(names{idx}, 'FontSize', 10);
    %xlim([0,9]);
    xlim([-4,24]);
    %legend({'KO low NIK','KO med NIK','KO high NIK','WT low NIK','WT med NIK','WT high NIK'});
    legend({'KO','WT'});
    set(gcf,'color','w');
    box off;
    set(gca,'TickDir','out')
end

%% plot the input curves for IkBa mutant
figure;plot((1:2881)/(2881/48),nfkb2DataWT.IKK(:,:)','LineWidth', 2,'color',[0,0,1]);
set(gcf,'color','w');
box off;
xlim([0,24]);
%ylim([0,50]);
title('IKK');
figure;
plot((1:2881)/(2881/48),nfkb2DataWT.NIK(1,:)','LineWidth', 2,'color',colorArray{1});
hold on;
% plot((1:2881)/(2881/48),IkBaData.NIK(2,:)','LineWidth', 2,'color',colorArray{2});
% plot((1:2881)/(2881/48),IkBaData.NIK(3,:)','LineWidth', 2,'color',colorArray{3});
set(gcf,'color','w');
box off;
xlim([0,24]);
%ylim([0,0.75]);
title('NIK');

%% save data
save 'lastRun.mat'

