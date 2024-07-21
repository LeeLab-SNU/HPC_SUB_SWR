Initial_SWRFilter_common;
warning off
ROOT.Old = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Save = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Fig = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\L0'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U0'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
if ~exist(ROOT.Save), mkdir(ROOT.Save); end
if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);
Cluster_List = readtable([ROOT.Units '\ClusterList_SWR_CA1_filtered.xlsx']);
thisRegion = 'CA1';
Experimenter = {'LSM'};
RipplesTable_all = table;

for sid=1:size(SessionList,1)
    if SessionList.include(sid) & ismember(SessionList.experimenter{sid},Experimenter)
        thisRID = SessionList.rat(sid);
        thisSID = SessionList.session(sid);
        thisRSID = [jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2)];
        load([ROOT.Behav '\' thisRSID '.mat']);
          BehavTable = readtable([ROOT.Behav '\' thisRSID '.xlsx']);
         
          Behav.t2 = [Behav.t(1):0.033:Behav.t(end)]';
          Behav.v2 = Behav.velocity.speed(knnsearch(Behav.t,Behav.t2));
          
          
          
        Recording_region_TT = Recording_region({thisRSID},:);
        TargetTT = GetTargetTT(ROOT,thisRSID,thisRegion,Params,1);
        if ~exist([ROOT.Save '\' thisRSID '_' thisRegion '.csv']), continue; end
        ripples = readtable([ROOT.Save '\' thisRSID '_' thisRegion '.csv']);
        EEG = LoadEEGData(ROOT, thisRSID, TargetTT,Params,Params_Ripple);
            Spike=LoadSpikeData(ROOT, thisRSID, [1:24], Params.cellfindn);
            
        clusters = Cluster_List(Cluster_List.rat==thisRID & Cluster_List.session==thisSID,:);
       RipplePropertySheet_Lv0_exc(ROOT,SessionList(sid,:),ripples,EEG,TargetTT,Spike,clusters,Behav,BehavTable,Params_Ripple);
       
       disp([thisRSID ' is finished!'])
    end
end
        