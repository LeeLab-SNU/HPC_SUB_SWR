%% Unit Pair_Field_analysis
Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat\R2'];
ROOT.Rip4 = [ROOT.Mother '\Processed Data\ripples_mat\R4_SUB_refCA1'];
ROOT.Rip5 = [ROOT.Mother '\Processed Data\ripples_mat\R5_cell_AllPopul'];
ROOT.Fig = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R25_ca1'];
ROOT.Unit1 = [ROOT.Mother '\Processed Data\units_mat\U1'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U2'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];


CList = [ [207 8 23]/255;[23 84 181]/255];


RegionList = {'SUB','CA1'};

RipplesTable.SUB = readtable([ROOT.Save '\RipplesTable_SUB_refCA1_forAnalysis' '.xlsx']);
ReactTable.SUB = readtable([ROOT.Save '\ReactTable_SUB_SUB.xlsx']);
UnitsTable.SUB = readtable([ROOT.Save '\UnitsTable_SUB_forAnalysis_TP.xlsx']);
UnitsTable_field.SUB = readtable([ROOT.Save '\UnitsTable_SUB_field_forAnalysis.xlsx']);

RipplesTable.CA1 = readtable([ROOT.Save '\RipplesTable_CA1_forAnalysis' '.xlsx']);
ReactTable.CA1 = readtable([ROOT.Save '\ReactTable_CA1_CA1.xlsx']);
UnitsTable.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_forAnalysis_TP.xlsx']);
UnitsTable_field.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_field_forAnalysis.xlsx']);


%%
% UnitPairF.CA1.Nsp_L(UnitPairF.CA1.Nsp_L>0.9999)=nan;
% UnitPairF.CA1.Nsp_R(UnitPairF.CA1.Nsp_R>0.9999)=nan;
% UnitPairF.CA1.Nsp_C(UnitPairF.CA1.Nsp_C>0.9999)=nan;
% 
% UnitPairF.SUB.Nsp_L(UnitPairF.SUB.Nsp_L>0.9999)=nan;
% UnitPairF.SUB.Nsp_R(UnitPairF.SUB.Nsp_R>0.9999)=nan;
% UnitPairF.SUB.Nsp_C(UnitPairF.SUB.Nsp_C>0.9999)=nan;
% 
%     writetable(UnitPairF.CA1,[ROOT.Save '\UnitPair_CA1_field.xlsx'],'writemode','replacefile');
%  writetable(UnitPairF.SUB,[ROOT.Save '\UnitPair_SUB_field.xlsx'],'writemode','replacefile');

UnitPairF.CA1 = readtable([ROOT.Save '\UnitPair_CA1_field.xlsx']);
UnitPairF.SUB = readtable([ROOT.Save '\UnitPair_SUB_field.xlsx']);
    UFca1 = UnitPairF.CA1;
    UFsub = UnitPairF.SUB;
    UFsub.Nsp_max = nanmax([UFsub.Nsp_L,UFsub.Nsp_R,UFsub.Nsp_C],[],2);
    UFca1.Nsp_max = nanmax([UFca1.Nsp_L,UFca1.Nsp_R,UFca1.Nsp_C],[],2);
        UFsub.Nsp_min = nanmin([UFsub.Nsp_L,UFsub.Nsp_R,UFsub.Nsp_C],[],2);
    UFca1.Nsp_min = nanmin([UFca1.Nsp_L,UFca1.Nsp_R,UFca1.Nsp_C],[],2);
        UFsub.Nsp_avg = nanmean([UFsub.Nsp_L,UFsub.Nsp_R,UFsub.Nsp_C],2);
    UFca1.Nsp_avg = nanmean([UFca1.Nsp_L,UFca1.Nsp_R,UFca1.Nsp_C],2);
%% scatterplot
figure;
s = 'C'; lincl='r';
subplot(1,2,1); 
x=UFsub.Sp; y=UFsub.(['Nsp_' s]);c=UFsub.SameUnit;
sc = scatterhistogram(x,y,'groupdata',c,'HistogramDisplayStyle','smooth', 'LineStyle','-','MarkerFilled','off','MarkerStyle','x');
sc.Color = {'Blue','Black'};
sc.BinWidths = .05;
mdl = fitlm(x,y);
% h=plot(mdl);
% delete(h(1))
% set(h(2),'LineWidth',3,'color',lincl)
% set(h(3),'LineWidth',2,'color',lincl)
% set(h(4),'LineWidth',2,'color',lincl)
x(isnan(y))=[]; y(isnan(y))=[];
y(isnan(x))=[]; x(isnan(x))=[];
[p1,p2] = corr(x,y);
% text(0, .8,['corr = ' jjnum2str(p1,3)],'color',lincl)
% text(0, .7,['p = ' num2str(p2)],'color',lincl)
xlim([-1 1]); ylim([-1 1]); xlabel('Spatial correlation'); ylabel(['Non-spatial correlation (' s ')']); title('SUB')
set(gca,'fontsize',12,'LegendTitle','Same unit')

subplot(1,2,2); 
x=UFca1.Sp; y=UFca1.(['Nsp_' s]);c=UFca1.SameUnit;
sc = scatterhistogram(x,y,'groupdata',c,'HistogramDisplayStyle','smooth', 'LineStyle','-','MarkerFilled','off','MarkerStyle','x');
sc.Color = {'Blue','Black'};
sc.BinWidths = .05;
mdl = fitlm(x,y);
% h=plot(mdl);
% delete(h(1))
% set(h(2),'LineWidth',3,'color',lincl)
% set(h(3),'LineWidth',2,'color',lincl)
% set(h(4),'LineWidth',2,'color',lincl)
x(isnan(y))=[]; y(isnan(y))=[];
y(isnan(x))=[]; x(isnan(x))=[];
[p1,p2] = corr(x,y);
% text(0, .8,['corr = ' jjnum2str(p1,3)],'color',lincl)
% text(0, .7,['p = ' num2str(p2)],'color',lincl)
xlim([-1 1]); ylim([-1 1]); xlabel('Spatial correlation'); ylabel(['Non-spatial correlation (' s ')']); title('CA1')
set(gca,'fontsize',12,'LegendTitle','Same unit')


%% SVM

rng(1); % 재현성을 위해

% SVM 분류기 학습
CVSVMModel = fitcsvm(X,Y,'Holdout',0.15,'ClassNames',{'b','g'},'Standardize',true);
CompactSVMModel = CVSVMModel.Trained{1}; % 학습된, 컴팩트한 분류기 추출

% 테스트 샘플의 레이블 예측
[label,score] = predict(CompactSVMModel,XTest);

% 결과 테이블 생성
table(YTest(1:10),label(1:10),score(1:10,2),'VariableNames',{'TrueLabel','PredictedLabel','Score'})


%%
thisRegion = 'CA1';

UPca1 = UnitPair.CA1; 
UPca1.crp = UPca1.pq ./ UPca1.un; UPca1.crp(UPca1.p==0 | UPca1.q==0) = nan;

d = 'Max';
figure; hold on

scatter3(UPca1.(['SpCorr_' d]),UPca1.(['NSpCorr_L_' d]),UPca1.crp,30,CList(2,:),'filled')
xlabel('Spatial corr'); ylabel('Non-spatial corr (L Scene)');  zlabel('Co-react. Prob.');



UPsub = UnitPair.SUB; 
UPsub.crp = UPsub.pq ./ UPsub.un; UPsub.crp(UPsub.p==0 | UPsub.q==0) = nan;
scatter3(UPsub.(['SpCorr_' d]),UPsub.(['NSpCorr_L_' d]),UPsub.crp,30,CList(1,:),'filled')
xlabel('Spatial corr'); ylabel('Non-spatial corr (L Scene)');  zlabel('Co-react. Prob.');


title([d ' corr'])
legend({'CA1','SUB'})

[h,p] = ttest2(UPca1.(['SpCorr_' d]),UPsub.(['SpCorr_' d]))

model = fitlm(UPsub, ['crp ~ SpCorr_' d ' + NSpCorr_L_' d ' + NSpCorr_R_' d ' + NSpCorr_C_' d]);
disp(model)
model.Rsquared.Adjusted

model = fitlm(UPca1, ['crp ~ SpCorr_' d ' + NSpCorr_L_' d ' + NSpCorr_R_' d ' + NSpCorr_C_' d]);
disp(model)
model.Rsquared.Adjusted