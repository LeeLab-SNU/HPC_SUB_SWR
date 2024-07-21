Initial_SWRFilter_common

thisRegion = 'CA1';
thisRegion2 = 'CA1';

ROOT.Save = [ROOT.Processed];
ROOT.Unit = [ROOT.Save '\units_mat\U2'];
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis.xlsx']);
% RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_field_RDIs_UV.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion2 '.xlsx']);
ROOT.Rip = [ROOT.Save '\ripples_mat\R4'];

if ~exist(ROOT.Rip), mkdir(ROOT.Rip); end


rng('shuffle');
% gcp = parpool(4);

Params.lineFittingMethod = 'linear regression'; % line finding

%% decoding parameters
Params.tbinDuration = 0.01;
Params.N_threshold = 3;
Params.fitting_threshold = 1.5;
a_test = 0.05;

%% variables for display
max_position = [1.1 1.1 0.1 0.1];
imagePosition = [1600 700 800 300];

ADbituVolts = 0.000000076296274187370727 * 10^6;
% Replay_Bayesian_decoding_batch

%% Ripple info
for RipNum = 1:size(RipplesTable,1)
    try
        [nCells,posterior,p_test,R_actual,R_shuffled,v,c] = ...
            Replay_Bayesian_decoding(ROOT,Params,RipplesTable(RipNum,:),ReactTable,5);
        save([ROOT.Rip '\' RipplesTable.ID{RipNum} '.mat'],"c","v",'R_shuffled','R_actual','posterior','p_test')
%         RipplesTable.nPCs(RipNum) = nCells;
        RipplesTable.DecodingP_all(RipNum) = p_test(1);
        if length(p_test)>1
        RipplesTable.DecodingP_Z(RipNum) = p_test(2);
        RipplesTable.DecodingP_P(RipNum) = p_test(3);
        RipplesTable.DecodingP_B(RipNum) = p_test(4);
        RipplesTable.DecodingP_M(RipNum) = p_test(5);
        end
    catch
        RipplesTable.DecodingP_all(RipNum) = nan;
        RipplesTable.DecodingP_Z(RipNum) = nan;
        RipplesTable.DecodingP_P(RipNum) = nan;
        RipplesTable.DecodingP_B(RipNum) = nan;
        RipplesTable.DecodingP_M(RipNum) = nan;
        disp([RipplesTable.ID{RipNum} ' error'])
    end
end

writetable(RipplesTable,[ROOT.Rip '\RipplesTable_' thisRegion '_ReplayP.xlsx'],'writemode','replacefile')
RipplesTable(isnan(RipplesTable.DecodingP_all),:)=[];
writetable(RipplesTable,[ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis.xlsx'],'writemode','replacefile')
