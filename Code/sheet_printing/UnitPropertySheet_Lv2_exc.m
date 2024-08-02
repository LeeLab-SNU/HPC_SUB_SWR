% ------------------ visual display
%

load([ROOT.Unit '\' clusterID '.mat']);
thisMap.d = d;
%% this field map 1D _ session
thisField = table;
temp=[];
temp = Spk.t_spk(~Spk.area_spk(:,5) & Spk.correctness_spk(:,1));
temp = sortrows(temp,1);

thisField.ts = temp(:,1);
thisMap.thisFieldMap1D = getFieldMaps(clusterID,thisField,'session',ROOT.Raw.Mother,ROOT.Info);
thisMap.thisFieldMap1D_trial = getFieldMaps(clusterID,thisField,'trial',ROOT.Raw.Mother,ROOT.Info);
FRMap = [];
for i=1:numel(thisMap.thisFieldMap1D.skaggsMap1D)
    FRMap(i,:) = thisMap.thisFieldMap1D.skaggsMap1D{i};
end

FRMap(6,:) = thisMap.thisFieldMap1D.skaggsMap_left1D;
FRMap(7,:) = thisMap.thisFieldMap1D.skaggsMap_right1D;
Dv =  (max(Pos.y)-diverging_point)/thisFRMapSCALE;

thisFieldMap1D_trial = thisMap.thisFieldMap1D_trial;
%%

% sheetPOS = [0 0 29.7 21.0];

szFONT = 10;
txtFONT = 14;
Px2cm=.23;

FigPos.Background = [100 100 1400 800];
FigPos.LinFRMaps = [.05 .65 .2 .2];
FigPos.TextBox = [.05 .9 .4 .1];
FigPos.ReactProp = [.05 .25 .2 .2];
FigPos.Raster = [.05 .06 .2 .4];
FigPos.RDI_bar = [.43 .92 .1 .08];
FigPos.NormReact = [.54 .56 .18 .04];
FigPos.ReactProbScatter = [.61 .48 .14 .14];

clist = {'#E72416','#4FBFD6','#711419','#4F63AD','#EF972C','#00506A'};
clist = {'#E72416','#4FBFD6','#711419','#4F63AD','#E72416','#4FBFD6'};
sheet = figure('position',FigPos.Background,'color','w');

switch sessiontype(5:end)
    case 'fc'
        CxtList = {'Forest','','City',''};
    case 'zpbm'
        CxtList = {'Zebra','Pebbles','Bamboo','Mountains'};
    case 'ds'
        CxtList = {'Dot','Square','',''};
    case 'dszp'
        CxtList = {'Dot','Square','Zebra','Pebbles'};
end


%% clusterID & printed date
subplot(5,6,1);
set(gca,'units','normalized','position',[0.03 0.97 0.94 0.02]); axis off;
text(0.02, 0.5, ['Rat' thisRID '-S' thisSID '-TT' thisTTID '-cluster' thisCLID],'FontSize',txtFONT);
text(0.02, -2, ['RDI:   Lscene=' jjnum2str(d(1),2) ', Rscene= ' jjnum2str(d(2),2) ', Choice= ' jjnum2str(d(3),2)],'FontSize',txtFONT);
text(0.85, 0.5, ['printed on ' date],'FontSize',szFONT);

%%
subplot('position',FigPos.LinFRMaps)
hold on
temp_trialFR = thisMap.thisFieldMap1D_trial;

DrawLinFR(FRMap(1,:),nanstd(temp_trialFR,[],2)/sqrt(size(temp_trialFR,2)),'k')
SetBasicPropFR(FRMap,Dv)
title(['overall'],'fontweight','b')

% if ~isnan(d(1))
subplot('position',FigPos.LinFRMaps+[.23 0 0 .05])
hold on
DrawLinFR(FRMap(2,:),nanstd(temp_trialFR(:,Pos.trial_context==1),[],2)/sqrt(size(temp_trialFR(:,Pos.trial_context==1),2)),hex2rgb(clist(1)))
l1 = plot(FRMap(2,:),'color',clist{1},'linewidth',1);
DrawLinFR(FRMap(3,:),nanstd(temp_trialFR(:,Pos.trial_context==3),[],2)/sqrt(size(temp_trialFR(:,Pos.trial_context==3),2)),hex2rgb(clist(3)))
l2 = plot(FRMap(3,:),'color',clist{3},'linewidth',1);
SetBasicPropFR(FRMap,Dv)
title(['left choice (Lscene RDI = ' num2str(round(d(1),3)) ')'],'fontweight','b')
legend([l1 l2],{CxtList{1},CxtList{3}},'location','northoutside','orientation','horizontal')
% end

% if ~isnan(d(2))
subplot('position',FigPos.LinFRMaps+[.46 0 0 .05])
hold on
DrawLinFR(FRMap(4,:),nanstd(temp_trialFR(:,Pos.trial_context==2),[],2)/sqrt(size(temp_trialFR(:,Pos.trial_context==2),2)),hex2rgb(clist(2)))
l1 = plot(FRMap(4,:),'color',clist{2},'linewidth',1);
DrawLinFR(FRMap(5,:),nanstd(temp_trialFR(:,Pos.trial_context==4),[],2)/sqrt(size(temp_trialFR(:,Pos.trial_context==4),2)),hex2rgb(clist(4)))
l2 = plot(FRMap(5,:),'color',clist{4},'linewidth',1);
SetBasicPropFR(FRMap,Dv)
title(['right choice (Rscene RDI = ' num2str(round(d(2),3)) ')'],'fontweight','b')
legend([l1 l2],{CxtList{2},CxtList{4}},'location','northoutside','orientation','horizontal')
% end

% if ~isnan(d(3))
subplot('position',FigPos.LinFRMaps+[.69 0 0 .05])
hold on
DrawLinFR(FRMap(6,:),nanstd(temp_trialFR(:,Pos.trial_side==1),[],2)/sqrt(size(temp_trialFR(:,Pos.trial_side==1),2)),hex2rgb(clist(5)))
l1 = plot(FRMap(6,:),'color',hex2rgb(clist(5)),'linewidth',1);
DrawLinFR(FRMap(7,:),nanstd(temp_trialFR(:,Pos.trial_side==2),[],2)/sqrt(size(temp_trialFR(:,Pos.trial_side==2),2)),hex2rgb(clist(6)))
l2 = plot(FRMap(7,:),'color',hex2rgb(clist(6)),'linewidth',1);
SetBasicPropFR(FRMap,Dv)
title(['left vs. right choice (choice RDI = ' num2str(round(d(3),3)) ')'],'fontweight','b')
legend([l1 l2],{'Left','Right'},'location','northoutside','orientation','horizontal')
% end


%%
FRMap(isnan(FRMap))=0;
colormap gray

subplot('position',FigPos.ReactProp+[0 .3 0 -.17])
imagesc(max(max(FRMap))-FRMap(1,:)); line([Dv Dv], [.5 1.5],'color','r','linewidth',1.5); yticks([]); xticks([]); caxis([0 max(max(FRMap))])
title('overall')

% if ~isnan(d(1))
subplot('position',FigPos.ReactProp+[.23 .3 0 -.17])
imagesc(max(max(FRMap))-FRMap(2,:)); line([Dv Dv], [.5 1.5],'color','r','linewidth',1.5); yticks([]); xticks([]); caxis([0 max(max(FRMap))])
title(CxtList{1})

subplot('position',FigPos.ReactProp+[.23 .23 0 -.17])
imagesc(max(max(FRMap))-FRMap(3,:)); line([Dv Dv], [.5 1.5],'color','r','linewidth',1.5); yticks([]); xticks([]); caxis([0 max(max(FRMap))])
title(CxtList{3})
% end

% if ~isnan(d(2))
subplot('position',FigPos.ReactProp+[.46 .3 0 -.17])
imagesc(max(max(FRMap))-FRMap(4,:)); line([Dv Dv], [.5 1.5],'color','r','linewidth',1.5); yticks([]); xticks([]); caxis([0 max(max(FRMap))])
title(CxtList{2})

subplot('position',FigPos.ReactProp+[.46 .23 0 -.17])
imagesc(max(max(FRMap))-FRMap(5,:)); line([Dv Dv], [.5 1.5],'color','r','linewidth',1.5); yticks([]); xticks([]); caxis([0 max(max(FRMap))])
title(CxtList{4})
% end

% if ~isnan(d(3))
subplot('position',FigPos.ReactProp+[.69 .3 0 -.17])
imagesc(max(max(FRMap))-FRMap(6,:)); line([Dv Dv], [.5 1.5],'color','r','linewidth',1.5); yticks([]); xticks([]); caxis([0 max(max(FRMap))])
title(['Left (' CxtList{1} ' & ' CxtList{3} ')'])

subplot('position',FigPos.ReactProp+[.69 .23 0 -.17])
imagesc(max(max(FRMap))-FRMap(7,:)); line([Dv Dv], [.5 1.5],'color','r','linewidth',1.5); yticks([]); xticks([]); caxis([0 max(max(FRMap))])
title(['Right (' CxtList{2} ' & ' CxtList{4} ')'])
% end

saveas(gca, [ROOT.Save2 '\rat' clusterID '_fr.svg']);
%%
subplot('position',FigPos.Raster+[.03 .0 -.08 -.23])

b = bar(d);
b.FaceColor = 'flat';
if d(1)>=0, b.CData(1,:) = hex2rgb(clist(1)); else, b.CData(1,:) = hex2rgb(clist(3)); end
if d(2)>=0, b.CData(2,:) = hex2rgb(clist(2)); else, b.CData(2,:) = hex2rgb(clist(4)); end
if d(3)>=0, b.CData(3,:) = hex2rgb(clist(5)); else, b.CData(3,:) = hex2rgb(clist(6)); end
set(gca,'ylim',[-1 1],'xticklabels',{'L','R','Ch'})
ylabel(sprintf(['Bamboo           RDI           Zebra\n'...
    'Mountains                          Pebbles\n'...
    'Right                         Left']))
% ylabel('RDI')
% axis ij

%%
y_linearized = get_linearized_position_JM(clusterID,x,y,Pos.side,diverging_point);
Spk.y_linearized = interp1(Pos.t,y_linearized,Spk.t_spk);
thisSpk = Spk.y_linearized(inspk);
[~,trial_temp] = max(Spk.trial_spk,[],2);
switch exper, case 'LSM', xlims=[-5 120]; case 'JS', xlims=[0 300]; otherwise, xlims=[0 100]; end

subplot('position',FigPos.Raster+[0 .27 0 -.2])

hold on
r=rasterplot_JM(trial_temp(inspk),Spk.y_linearized(inspk),'k',1);
SetBasicPropRaster(xlims,diverging_point,trial_temp)

% if ~isnan(d(1))
subplot('position',FigPos.Raster+[.23 0 0 0])
hold on
t1 = find(Pos.trial_context==1); t2=find(Pos.trial_context==3);
[x1,r1] = sort_unique(trial_temp(inspk& Spk.cont_spk(:,1)),t1);
r=rasterplot_JM(x1,Spk.y_linearized(inspk & Spk.cont_spk(:,1)),hex2rgb(clist(1)),1);
[x2,r2] = sort_unique(trial_temp(inspk& Spk.cont_spk(:,3)),t2);
x2=x2+length(t1)+2;
r=rasterplot_JM(x2,Spk.y_linearized(inspk & Spk.cont_spk(:,3)),hex2rgb(clist(3)),1);
line([xlims],[length(t1)+.5 length(t1)+.5],'color','k','linestyle','--')
SetBasicPropRaster(xlims,diverging_point,trial_temp)
set(gca, 'ylim',[0 length(find(Pos.trial_side==1))],'ytick',[unique(x1);unique(x2)]-.5,'yticklabels', [r1;r2])
set(gca, 'ylim',[0 length(find(mod(Pos.trial_context,2)==1))+2],'ytick',[1;length(t1);length(t1)+2;length([t1;t2])+2]-.5,...
    'yticklabels', [1;length(t1);1;length(t2)])
% end
% if ~isnan(d(2))
subplot('position',FigPos.Raster+[.46 0 0 0])
hold on
t1 = find(Pos.trial_context==2); t2=find(Pos.trial_context==4);
[x1,r1] = sort_unique(trial_temp(inspk& Spk.cont_spk(:,2)),t1);
r=rasterplot_JM(x1,Spk.y_linearized(inspk & Spk.cont_spk(:,2)),hex2rgb(clist(2)),1);
[x2,r2] = sort_unique(trial_temp(inspk& Spk.cont_spk(:,4)),t2);
x2=x2+length(t1)+2;
r=rasterplot_JM(x2,Spk.y_linearized(inspk & Spk.cont_spk(:,4)),hex2rgb(clist(4)),1);
line([xlims],[length(t1)+.5 length(t1)+.5],'color','k','linestyle','--')
SetBasicPropRaster(xlims,diverging_point,trial_temp)
set(gca, 'ylim',[0 length(find(mod(Pos.trial_context,2)==0))+2],'ytick',[1;length(t1);length(t1)+2;length([t1;t2])+2]-.5,...
    'yticklabels', [1;length(t1);1;length(t2)])
% end
% if ~isnan(d(3))
subplot('position',FigPos.Raster+[.69 0 0 0])
hold on
t1 = find(mod(Pos.trial_context,2)==1); t2=find(mod(Pos.trial_context,2)==0);
[x1,r1] = sort_unique(trial_temp(inspk& Spk.side_spk(:,1)),t1);
r=rasterplot_JM(x1,Spk.y_linearized(inspk & Spk.side_spk(:,1)),hex2rgb(clist(5)),1);
[x2,r2] = sort_unique(trial_temp(inspk& Spk.side_spk(:,2)),t2);
x2=x2+length(t1)+4;
r=rasterplot_JM(x2,Spk.y_linearized(inspk & Spk.side_spk(:,2)),hex2rgb(clist(6)),1);
line(xlims,[length(t1)+1.5 length(t1)+1.5],'color','k','linestyle','--')
SetBasicPropRaster(xlims,diverging_point,trial_temp)
set(gca, 'ylim',[0 length(find(Pos.trial_side~=0))+4],'ytick',[1;length(t1);length(t1)+4;length([t1;t2])+4]-.5,...
    'yticklabels', [1;length(t1);1;length(t2)])
% end
%%


%%

    cd(ROOT.Save2)
    saveas(gca, [ROOT.Save2 '\rat' clusterID '.svg']);
close all
%%
function DrawRawSpkMap(Pos,Spk,Title,num,xlim,ylim,in,inspk,exper)
if ~isempty(Title{num})
    
    scatter(Pos.x(in) - xlim(1), Pos.y(in)  - ylim(1)+1,10,[.8 .8 .8],'filled')
    hold on
    plot(Spk.x_spk(inspk&Spk.cont_spk(:,num)) - xlim(1), Spk.y_spk(inspk&Spk.cont_spk(:,num)) - ylim(1)+1, '.', 'MarkerSize', 2, 'color', 'k');
    title(sprintf([Title{num} '\n(Nspks=' num2str(length(Spk.x_spk(inspk&Spk.cont_spk(:,num)))) ')']),'FontWeight','Bold','FontSize',10)
    
    if strcmp(exper,'JS')
        set(gca,'XLim', [0 xlim(2)-xlim(1)], 'YLim', [0 ylim(2)-ylim(1)],'xtick',[],'ytick',[]);
    else
        set(gca, 'YDir', 'rev','XLim', [0 xlim(2)-xlim(1)], 'YLim', [0 ylim(2)-ylim(1)],'xtick',[],'ytick',[]);
    end
end
axis off;

end
function DrawLinFR(m,error,c)
m(isnan(m))=0;
error(isnan(error))=0;
plot(m,'color',c,'linewidth',2)
    x_vector = [[1:length(m)]'; flip([1:length(m)]')];
    patch = fill(x_vector, [m'+error; flip(m'-error)], c);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', .3);
    s1=ScatMaxFR(m);
end
function s = ScatMaxFR(FRMap)
[n,m] = max(FRMap);
s=scatter(m,n,20,'r','filled');
text(m*1.1,n,[jjnum2str(n,2) 'Hz'])
end

function SetBasicPropFR(FRMap,Dv)
line([Dv Dv],[0 max(max(FRMap))],'color','r','linewidth',1.5)
set(gca,'ylim',[0 max(max(FRMap))*1.1])
if size(FRMap,2)-Dv<3
    xticks([1 size(FRMap,2)]); xticklabels({'Startbox','End'});
else
xticks([1 Dv size(FRMap,2)]); xticklabels({'Startbox','Dv','Fdwell'});
end
xlim([1 size(FRMap,2)])
ylabel('FR (Hz)')
end

function SetBasicPropRaster(xlims,diverging_point,trial_temp)
line([diverging_point,diverging_point],[0 max(trial_temp)],'color','r','linewidth',1.5)
xlabel('linearized y (cm)')
ylabel('trial')
set(gca,'XDir','reverse','xlim',xlims)
end