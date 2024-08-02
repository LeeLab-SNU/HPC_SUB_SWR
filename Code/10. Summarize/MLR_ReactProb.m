prob = RipplesTable.SUB.DecodingP_all;
Prob2 = RipplesTable.CA1.DecodingP_all;

Sp1 = RipplesTable.SUB.DecodingP_all;
Sp2 = RipplesTable.CA1.DecodingP_all;


Nsp1 = nanmin([RipplesTable.SUB.pRDI_L_UV,RipplesTable.SUB.pRDI_R_UV,RipplesTable.SUB.pRDI_C_UV],[],2);
Nsp2 =  nanmin([RipplesTable.CA1.pRDI_L_UV,RipplesTable.CA1.pRDI_R_UV,RipplesTable.CA1.pRDI_C_UV],[],2);

figure;
% subplot(1,2,1)
% Nsp1(isnan(Nsp1))=1; Nsp2(isnan(Nsp2))=1;
% cdfplot(Sp1); hold on; cdfplot(Sp2); 
% 
% 
% subplot(1,2,2)
Nsp1(isnan(Nsp1))=rand([sum(isnan(Nsp1)),1])/2+0.5; Nsp2(isnan(Nsp2))=rand([sum(isnan(Nsp2)),1])/2+0.5;
cdfplot(Nsp1); hold on; cdfplot(Nsp2); 
xlim([0 .5])


[h,p] = kstest2(Sp1,Sp2)
[h,p] = kstest2(Nsp1,Nsp2)


[h,p] = ttest2(Sp1,Sp2)

[h,p] = kstest2(Nsp1,Nsp2)

nanmean(Sp1)


%%
figure
subplot(1,2,1)
% pie([sum(Sp1<0.05 & ~(Nsp1<0.05)) sum(Sp1<0.05 & Nsp1<0.05) sum(~(Sp1<0.05) & Nsp1<0.05) sum(~(Sp1<0.05) & ~(Nsp1<0.05))])
pie([sum(~(Sp1<0.05) & Nsp1<0.05) sum((Sp1<0.05) | ~(Nsp1<0.05))])
subplot(1,2,2)
% pie([sum(Sp2<0.05 & ~(Nsp2<0.05)) sum(Sp2<0.05 & Nsp2<0.05) sum(~(Sp2<0.05) & Nsp2<0.05) sum(~(Sp2<0.05) & ~(Nsp2<0.05))])
pie([sum(~(Sp2<0.05) & Nsp2<0.05) sum((Sp2<0.05) | ~(Nsp2<0.05))])

%%
UnitPairT=struct;
UnitPairT.CA1 = load([ROOT.Save '\UnitPair_CA1.mat']); UnitPairT.CA1= UnitPairT.CA1.UnitPair;
UnitPairT.SUB = load([ROOT.Save '\UnitPair_SUB.mat']); UnitPairT.SUB= UnitPairT.SUB.UnitPair;

UnitPairT.SUB(UnitPairT.SUB.p==0 | UnitPairT.SUB.q==0,:)=[];

%%
% 데이터를 불러옵니다.

UP = UnitPair.CA1; 
UP.crp = UP.pq ./ UP.un; UP.crp(UP.p==0 | UP.q==0) = nan;

data_o = table;
data_o.prob = UP.crp; 
data_o.scorr = UP.SpCorr_Max; 
data_o.ncorr = nanmax([UP.NSpCorr_L_Max,UP.NSpCorr_R_Max,UP.NSpCorr_C_Max],[],2);

% bins = linspace(min(data_o.scorr),max(data_o.scorr),200);
bins = linspace(-1,1,200);
movewin = 0.2;
acc_bins = zeros(length(bins),3);
corr_bins = zeros(length(bins),2);
for b = 1:length(bins)
% 연속변수 prob의 예측 모델 생성 (임의의 예측 모델)
% data = data_o(find(data_o.scorr>bins(b)-movewin & data_o.scorr<=bins(b)+movewin),:);
data = data_o(find(data_o.scorr>bins(b)),:); 

if isempty(data), continue; end
data(isnan(data.prob),:)=[];
prob = data.prob; scorr=data.scorr; ncorr = data.ncorr;

if size(data,1)<10, acc_bins(b) = nan; continue; end
% training set과 test set 분리

% 반복 횟수
nIter = 100;

% 정확도 초기화
accuracy = zeros(nIter,3);

parfor i = 1:nIter
    % 데이터 분할
    cv = cvpartition(height(data), 'HoldOut', 0.2);
    idx = cv.test;
    dataTrain = data(~idx,:);
    dataTest  = data(idx,:);

    % 선형 회귀 모델 생성
    mdl1 = fitlm(dataTrain, 'prob ~ scorr');
    mdl2 = fitlm(dataTrain, 'prob ~ ncorr');
    mdl3 = fitlm(dataTrain, 'prob ~ ncorr + scorr');
%     model  = fitglm(ncorr,prob);
%     [~,dev_noconstant] = glmfit(ones(length(prob),1),prob,'normal','Constant','off');
%     D = dev_noconstant - dev;
% prediction = predict(model,ncorr);

% [ConfusionMat,order] = confusionmat(prob,prediction)

    % 테스트 세트에 대한 예측
    yPred1 = predict(mdl1, dataTest); 
    yPred2 = predict(mdl2, dataTest);
    yPred3 = predict(mdl3, dataTest);
     accuracy(i,:) = [sqrt(sum((yPred1 - dataTest.prob).^2)/length(yPred1)),...
         sqrt(sum((yPred2 - dataTest.prob).^2)/length(yPred2)),...
         sqrt(sum((yPred3 - dataTest.prob).^2)/length(yPred3))];

% accuracy(i) = (mdl.Rsquared.Adjusted);
% 
% corrcoef(i) = corr(dataTrain.prob,dataTrain.scorr);

end

% 평균 예측 정확도 출력
acc_bins(b,:) = nanmean(accuracy);
corr_bins(b,:) = [corr(prob,scorr) corr(prob,ncorr)];
end
%%
% acc_ca1 = acc_bins; corr_ca1 = corr_bins;
figure
% acc_bins(isnan(acc_bins)) = [];
plot(smooth(acc_sub(:,2)));
hold on
plot(smooth(acc_ca1(:,2)));
xticks([1:10:200])
xticklabels(string([-1:0.1:1]))
% title('CA1')
xlim([00 200]);
% ylim([-0.01 .15])
disp(['Average test set accuracy: ', num2str(mean(accuracy))]);

%%
% 데이터를 로드합니다.
data = readtable([ROOT.Save '\UnitPair_All.xlsx']);

data.Nsp_Max = nanmax([data.Nsp_L,data.Nsp_R,data.Nsp_C],[],2);
data.Nsp_min = nanmin([data.Nsp_L,data.Nsp_R,data.Nsp_C],[],2);

% a의 범위를 정의합니다. 여기서는 a가 0부터 100까지 10단위로 변한다고 가정하였습니다.
Sp_values = 0:0.5:1;

% 각 변수에 대한 R-squared 값을 저장할 구조체를 초기화합니다.
coef_sub = []; coef_ca1 = [];
 
for r=1:2
    coef=[];
for i = 1:length(Sp_values)-1
    % a가 특정 값 이상인 데이터만 선택합니다.
    data_subset = data(data.Sp >= Sp_values(i) & data.Sp < Sp_values(i+1) & data.region==r, :);

    dat = [ones(size(data_subset,1),1),data_subset.Sp,...
        data_subset.Nsp_L,data_subset.Nsp_R,data_subset.Nsp_C,data_subset.Nsp_Max,data_subset.Nsp_min];
%     if isempty(data_subset), continue; end

        % GLM 모델을 만듭니다.
        b = regress(data_subset.crp,dat);
%         model = fitlm(data_subset, 'crp ~ Sp + Nsp_L + Nsp_R + Nsp_C + Nsp_Max + Nsp_min');

        % R-squared 값을 저장합니다.
        coef = [coef;b'];

end
if r==1, coef_ca1=coef; else, coef_sub=coef; end
end

%% 각 변수에 대한 R-squared 값의 변화를 line plot으로 그립니다.
figure
hold on
plot(Sp_values, coef_sub(:,2), 'DisplayName', 'Sp','linewidth',2)
plot(Sp_values, coef_sub(:,3), 'DisplayName', 'Nsp-L','linewidth',2)
plot(Sp_values, coef_sub(:,4), 'DisplayName', 'Nsp-R','linewidth',2)
plot(Sp_values, coef_sub(:,5), 'DisplayName', 'Nsp-C','linewidth',2)
plot(Sp_values, coef_sub(:,6), 'DisplayName', 'Nsp-Max','linewidth',2)
% plot(Sp_values, coef_sub(:,7), 'DisplayName', 'Nsp-min','linewidth',2)


xlabel('Sp value')
ylabel('coefficient')
title('Change in MLR coefficient in SUB')
legend('location','northwest')
line([0 .9],[0 0],'color','k')
hold off
set(gca,'fontsize',12,'fontweight','b')
xlim([0 .9]); ylim([-1 1])
% ylim([0 1])
figure
hold on
plot(Sp_values, coef_ca1(:,2), 'DisplayName', 'Sp','linewidth',2)
plot(Sp_values, coef_ca1(:,3), 'DisplayName', 'Nsp-L','linewidth',2)
plot(Sp_values, coef_ca1(:,4), 'DisplayName', 'Nsp-R','linewidth',2)
plot(Sp_values, coef_ca1(:,5), 'DisplayName', 'Nsp-C','linewidth',2)
plot(Sp_values, coef_ca1(:,6), 'DisplayName', 'Nsp-Max','linewidth',2)
% plot(Sp_values, coef_ca1(:,7), 'DisplayName', 'Nsp-min','linewidth',2)


xlabel('Sp value')
ylabel('coefficient')
title('Change in MLR coefficient in CA1')
legend('location','northwest')
line([0 .9],[0 0],'color','k')
hold off
set(gca,'fontsize',12,'fontweight','b')
xlim([-0 .9]); ylim([-1 1])
%%
figure;
b = bar([coef_sub(2,2:6);coef_ca1(2,2:6)]');
b(1).FaceColor = CList(1,:);
b(2).FaceColor = CList(2,:);

xticks([1:5]); ylim([-.15 .15])
xticklabels({'Spatial', 'LScene', 'RScene','Choice','NonSpatial-Max'})
ylabel('coefficient')
legend({'SUB','CA1'},'location','eastoutside')
set(gca,'fontsize',12,'fontweight','b')
%% UP data 전처리
n=10;  thisRegion = 'CA1';
UP = UnitPair.(thisRegion);
UP.crp = UP.pq ./ UP.un; UP.crp(UP.p==0 | UP.q==0) = nan;
UP.NSpCorr_M_Max = nanmax([UP.NSpCorr_L_Max,UP.NSpCorr_R_Max,UP.NSpCorr_C_Max],[],2);
UP =  UP(UP.crp>0 & UP.p>n & UP.q>n,:); UP(isnan(UP.crp),:)=[];
UnitPair.(thisRegion) = UP;

Utemp=table;
Utemp.UID = [UP.UID1;UP.UID2]; Utemp.p = [UP.p;UP.q];
[x,ia,ic] = unique(Utemp.UID);

Utemp = Utemp(ia,:);

figure;
histogram(Utemp.p,'binwidth',5)
xlabel('Ripple 참여 횟수'); ylabel('Unit 수'); ylim([0 50])
 sum(Utemp.p>=10)
%% scatterplot with ANCOVA
% [X,ia,ic] = unique(UP.UID1);
% UP = UP(ia,:);
% sum(UP.p<10)
% figure; histogram(UP.p,'binwidth',5)
% xlabel('참여한 ripple의 수')

lincl = 'k'; n=10; s=2;
xvar = 'SpCorr_Max'; yvar = 'crp';
ImgName = ['Sp_' 'In_SFSF'];
if strcmp(xvar, 'SpCorr_Max'), xlab = 'Spatial Correlation (max)'; 
elseif strcmp(xvar, 'NSpCorr_M_Max'), xlab = 'Non-Spatial Correlation (max)'; 
end

figure('position',[263,165,1167,548]); 
subplot(1,2,1)
hold on
UP = UnitPair.SUB;
UP =  UP(UP.crp>0 & UP.p>n & UP.q>n,:); UP(isnan(UP.crp),:)=[];
sid1 = UP.numFields_1==1;
sid2 = UP.numFields_2==1;
sid = sid1+sid2;
UP.(yvar)(sid~=s)=nan;
x=UP.(xvar); y=UP.(yvar);
% scatter(x,y,40,CList(1,:),'filled')
scatter(x(~sid1 & ~sid2),y(~sid1 & ~sid2),40,CList(1,:),'filled')
scatter(x(logical(abs(sid1 - sid2))),y(logical(abs(sid1 - sid2))),40,CList(1,:))
scatter(x(sid1 & sid2),y(sid1 & sid2),40,'k','x')
mdl = fitlm(x,y);
h=plot(mdl);
delete(h(1))
set(h(2),'LineWidth',3,'color',lincl)
set(h(3),'LineWidth',2,'color',lincl)
set(h(4),'LineWidth',2,'color',lincl)
x(isnan(y))=[]; y(isnan(y))=[];
y(isnan(x))=[]; x(isnan(x))=[];
[p1,p2] = corr(x,y);
text(-.5, .8,['corr = ' jjnum2str(p1,3)],'color',lincl)
text(-.5, .7,['p = ' num2str(p2)],'color',lincl)
title('SUB'); xlabel(xlab); ylabel('Co-React. Prob.');
xlim([-1 1]); ylim([0 1])
set(gca,'fontsize',12,'fontweight','b'); legend off
corr(x,y)
legend({'MF-MF','SF-MF','SF-SF'})
UP = UP(~(sid~=s),:);
UP_SUB=UP;

subplot(1,2,2)
hold on
UP = UnitPair.CA1;
UP =  UP(UP.crp>0 & UP.p>n & UP.q>n,:); UP(isnan(UP.crp),:)=[];
sid1 = UP.numFields_1==1;
sid2 = UP.numFields_2==1;
sid = sid1+sid2;
UP.(yvar)(sid~=s)=nan;
x=UP.(xvar); y=UP.(yvar);
% scatter(x,y,40,CList(2,:),'filled')
scatter(x(~sid1 & ~sid2),y(~sid1 & ~sid2),40,CList(2,:),'filled')
scatter(x(logical(abs(sid1 - sid2))),y(logical(abs(sid1 - sid2))),40,CList(2,:))
scatter(x(sid1 & sid2),y(sid1 & sid2),40,'k','x')
h=plot(mdl);
delete(h(1))
set(h(2),'LineWidth',3,'color',lincl)
set(h(3),'LineWidth',2,'color',lincl)
set(h(4),'LineWidth',2,'color',lincl)
x(isnan(y))=[]; y(isnan(y))=[];
y(isnan(x))=[]; x(isnan(x))=[];
[p1,p2] = corr(x,y);
text(-.5, .8,['corr = ' jjnum2str(p1,3)],'color',lincl)
text(-.5, .7,['p = ' num2str(p2)],'color',lincl)
title('CA1'); xlabel(xlab); ylabel('Co-React. Prob.'); 
xlim([-1 1]); ylim([0 1])
set(gca,'fontsize',12,'fontweight','b'); legend off
corr(x,y)
legend({'MF-MF','SF-MF','SF-SF'})
UP = UP(~(sid~=s),:);
UP_CA1=UP;

% ANCOVA 
data1 = UP_SUB;
data2 = UP_CA1;

% 데이터 생성
x1 = data1.(xvar);
y1 = data1.(yvar);
x2 = data2.(xvar);
y2 = data2.(yvar);

% 두 데이터 세트 결합
x = [x1; x2];
y = [y1; y2];

% 그룹 변수 생성
group = [ones(size(x1)); 2*ones(size(x2))];

% ANCOVA 수행
[p,~,~] = anovan(y,{group,x},'varnames',{'Region','Spatial Corr.'},'continuous',2, 'model','interaction');
if p(3)<0.001, pstr='p<0.001'; else, pstr = ['p=' jjnum2str(p(3),3)]; end
sgtitle([pstr '; ANCOVA'],'fontsize',12,'fontweight','b')

% 결과 출력
% fprintf('p-value for the interaction between group and x: %.4f\n', p(3));


    saveas(gca,['D:\HPC-SWR project\Manuscript figures\fig7\UnitPair_nFields\' ImgName '.png'])
        saveas(gca,['D:\HPC-SWR project\Manuscript figures\fig7\UnitPair_nFields\' ImgName '.svg'])
%%

%%
lincl = 'k'; n=10; 
figure('position',[263,165,1167,548]); 
subplot(1,2,1); hold on

UP = UnitPair.SUB;
b = MkBar_UP_nFields(UP,'NSpCorr_M_Max',n);
b.FaceColor = 'flat';
b.FaceColor = CList(1,:);
title('SUB')

subplot(1,2,2); hold on
UP = UnitPair.CA1;
b = MkBar_UP_nFields(UP,'NSpCorr_M_Max',n);
b.FaceColor = 'flat';
b.FaceColor = CList(2,:);
title('CA1')


%% field correlation scatter & line plot
n=10; y= 'Sp'; thisRegion = 'CA1'; y2 = 'SpCorr_Max';
UP = UnitPair.(thisRegion); UPF = UnitPair_field.(thisRegion); 
UP.crp = UP.pq ./ UP.un; UP.crp(UP.p==0 | UP.q==0) = nan;
UP.NSpCorr_M_Max = nanmax([UP.NSpCorr_L_Max,UP.NSpCorr_R_Max,UP.NSpCorr_C_Max],[],2);
UP =  UP(UP.crp>0 & UP.p>n & UP.q>n,:); UP(isnan(UP.crp),:)=[];
UP.minNumFields = min([UP.numFields_1,UP.numFields_2],[],2);
UP = sortrows(UP,{'minNumFields',y2});

figure; hold on
for uid = 1:size(UP,1)
    fid = find(strncmp(UPF.UID1,UP.UID1{uid},12) & strncmp(UPF.UID2,UP.UID2{uid},12));
    UP_thisFields = UPF(fid,:);
    scatter(UP_thisFields.(y),uid*ones(size(UP_thisFields,1)),20,'k','.')
    if size(UP_thisFields,1)>1
        line([nanmin(UP_thisFields.(y)) nanmax(UP_thisFields.(y))],[uid uid],'color','k')
    end
end

%%
rng('default');  % 결과를 재현할 수 있도록 난수 생성기 초기화
x1 = rand(100,1);  % data1의 x1 생성
y1 = 3*x1 + randn(100,1);  % data1의 y1 생성 (회귀계수는 3)
x2 = rand(100,1);  % data2의 x2 생성
y2 = 2*x2 + randn(100,1);  % data2의 y2 생성 (회귀계수는 2)

% 두 데이터 세트 결합
x = [x1; x2];
y = [y1; y2];

% 그룹 변수 생성
group = [ones(size(x1)); 2*ones(size(x2))];

% ANCOVA 수행
[p,tbl,stats] = anovan(y,{group,x},'varnames',{'group','x'},'continuous',2);

% 결과 출력
fprintf('p-value for the interaction between group and x: %.4f\n', p(3));
%%
function b = MkBar_UP_nFields(UP,y,n)
if strncmp(y,'SpCorr',6), ylab = 'Spatial correlation';
elseif strncmp(y,'NSpCorr',7), ylab = 'In-field scene/choice selectivity correlation (max)';
elseif strcmp(y,'crp'), ylab = 'Co-reactivation probability'; end
UP.crp = UP.pq ./ UP.un; UP.crp(UP.p==0 | UP.q==0) = nan;
UP.NSpCorr_M_Max = nanmax([UP.NSpCorr_L_Max,UP.NSpCorr_R_Max,UP.NSpCorr_C_Max],[],2);
UP =  UP(UP.crp>0 & UP.p>n & UP.q>n,:); UP(isnan(UP.crp),:)=[];
sid1 = UP.numFields_1==1;
sid2 = UP.numFields_2==1;
x = UP.(y); 

x3 = x(~sid1 & ~sid2);
x2 = x(logical(abs(sid1 - sid2)));
x1 = x(sid1 & sid2);

b = bar([mean(x1) mean(x2) mean(x3)]');

err = [nanstd(x1)/sqrt(length(x1)) nanstd(x2)/sqrt(length(x2)) nanstd(x3)/sqrt(length(x3))];
    errorbar([1 2 3], [mean(x1) mean(x2) mean(x3)], err, 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);


text(1,mean(x1)+0.05,['n=' num2str(length(x1))]);...
    text(2,mean(x2)+0.05,['n=' num2str(length(x2))]);text(3,mean(x3)+0.05,['n=' num2str(length(x3))])
[~,p1] = ttest2(x2,x3); [~,p2] = ttest2(x2,x1); 
text(2.5,.95,['p=' jjnum2str(p1,2)]); text(1.5,.95,['p=' jjnum2str(p2,2)])
ylim([0 1]); 
ylabel(ylab)
xticks([1 2 3]); xticklabels({'SF-SF','SF-MF','MF-MF'})
end
