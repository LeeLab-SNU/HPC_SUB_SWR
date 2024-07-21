Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Save '\ripples_mat\R0'];
ROOT.Rip = [ROOT.Save '\ripples_mat\R3'];
ROOT.Rip4 = [ROOT.Save '\ripples_mat\R4_10ms'];
ROOT.Fig3 = [ROOT.Save '\ripples_mat\ProfilingSheet\R11_AI_hot'];
ROOT.Units = [ROOT.Save '\units_mat\U1'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];

dir = '';
ROOT.Fig = [ROOT.Fig3];

if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


thisRegion = 'CA1';
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis_RDI.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '.xlsx']);
UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion '.xlsx']);
UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable = UnitsTable_A;
TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

SessionTable  = readtable([ROOT.Save '\SessionCountTable_' thisRegion '_forAnalysis_RDI.xlsx']);

%%
close all
unit_sort = 'CorrectTrials';
ripple_sort = 'STindex';
Units_all=[];
for sid = 1:size(SessionTable,1)
    thisRID = SessionTable.rat(sid);
    thisSID = SessionTable.session(sid);
% Ripples = RipplesTable;
% Units = UnitsTable;
    Ripples = RipplesTable(RipplesTable.rat==thisRID & RipplesTable.session==thisSID,:);

Ripples = Ripples(Ripples.StartTime_fromTrialEnd<10,:);
Units = UnitsTable(UnitsTable.rat==thisRID & UnitsTable.session==thisSID,:);
Units.abs_RDI_C = abs(Units.RDI_LR);
%%
Units = sortrows(Units,unit_sort);
Ripples = sortrows(Ripples, ripple_sort);

RC = Ripples(Ripples.correctness==1,:);
RW = Ripples(Ripples.correctness==2,:);
%%
% figure;
% histogram(RC.StartTime_fromTrialEnd,'normalization','probability','binwidth',.5)
% hold on
% histogram(RW.StartTime_fromTrialEnd,'normalization','probability','binwidth',.5)
% xlim([0 25])
%%
URC=zeros(size(Units,1),size(RC,1));

for r=1:size(RC,1)
    thisRipUnits = ReactTable.UnitID(find(strcmp(RC.ID(r),ReactTable.RippleID)));
    for u=1:size(Units,1)
        if(ismember(Units.ID(u),thisRipUnits))
            URC(u,r)=1;
        end
    end
end

URW=zeros(size(Units,1),size(RW,1));
for r=1:size(RW,1)
    thisRipUnits = ReactTable.UnitID(find(strcmp(RW.ID(r),ReactTable.RippleID)));
    for u=1:size(Units,1)
        if(ismember(Units.ID(u),thisRipUnits))
            URW(u,r)=1;
        end
    end
end

%%
UR = nan(size(Ripples,1),size(Units,1));

for r=1:size(Ripples,1)
    thisRipUnits = ReactTable.UnitID(find(strcmp(Ripples.ID(r),ReactTable.RippleID)));
    n = Ripples.correctness(r);
    s = mod(Ripples.context(r),2);
    for u=1:size(Units,1)
        if(ismember(Units.ID(u),thisRipUnits))
            UR(r,u)=n;
        end
    end
end
%%
for u=1:size(Units,1)
    Units.nRipples(u) = sum(UR(:,u)>0);
    Units.CorrectTrials(u) = 1-(nanmean(UR(:,u))-1);
    Units.CorrectTrials_norm(u) = Units.CorrectTrials(u)/(1-mean(Ripples.correctness-1));
end
Units_all = [Units_all;Units];
    %%

figure('position',[2043,50,498,940])

imagesc(UR)
hold on
for r = 1:size(Ripples,1)-1
if Ripples.tr_temp(r)<Ripples.tr_temp(r+1)
    line([0 u],[r+.5 r+.5],'color','k')
end
end
u = find(isnan(Units.(unit_sort)),1,"first");
if isempty(u), u=size(Units,1)+1; end
% line([u-.5 u-.5],[0 size(Ripples,1)],'color','r')

u = find(Units.(unit_sort)>(1-mean(Ripples.correctness-1)),1,"first");
if isempty(u), u=size(Units,1)+1; end
 line([u-.5 u-.5],[0 size(Ripples,1)],'color','k')

clim([0 3])
title([num2str(thisRID) '-' jmnum2str(thisSID,2)])
xlabel(['Unit (' unit_sort ' sort)'],'Interpreter','none')
ylabel(['Ripple (' ripple_sort ' sort)'],'Interpreter','none')
colormap([1 1 1;0 0 0;1 0 0; 0 0 1])
% figure;
% imagesc([URC URW])

% saveas(gca,['D:\HPC-SWR project\temp\' num2str(thisRID) '-' num2str(thisSID) '.png'])
close all
end
%%

Units = Units_all;
Units = Units(Units.nRipples>4,:);
Units = Units(~isnan(Units.CorrectTrials),:);
Units = Units(~isnan(Units.RDI_LR),:);

figure;
histogram(Units.CorrectTrials_norm)

figure;
subplot(1,2,1)
scatter(Units.RDI_LR, Units.CorrectTrials,20,'k','filled')
 set(gca,'fontsize',12,'fontweight','b')
subplot(1,2,2)
scatter(abs(Units.RDI_LR), Units.CorrectTrials,20,'k','filled')
 set(gca,'fontsize',12,'fontweight','b')
corrcoef(abs(Units.RDI_LR),Units.CorrectTrials)

y = Units.CorrectTrials_norm;
X1 = abs(Units.RDI_LR);
x1 = ones(size(X1,1),1);
X = [x1 X1];

[~,~,~,~,stats] = regress(y,X)