% UnitPropertySheet_Lv2


Initial_SWRFilter_common;
ROOT.Save0 = [ROOT.Save '\units_mat\ProfilingSheet\U0'];
ROOT.Save1 = [ROOT.Save '\units_mat\ProfilingSheet\U1'];
ROOT.Save2 = [ROOT.Save '\units_mat\ProfilingSheet\U2'];
if ~exist(ROOT.Save2), mkdir(ROOT.Save2); end
ROOT.Unit = [ROOT.Save '\units_mat\U1'];

Session_List = readtable([ROOT.Info '\SessionList_SWR.xlsx']);
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv']);

Cluster_List = readtable([ROOT.Info '\ClusterList_SWR_CA1_filtered.xlsx']);
exper = {'LSM','JS','SEB'};
exper = {'LSM'};
TargRegion = 'CA1';
 
for cid = 1:size(Cluster_List,1)
if ismember(Cluster_List.experimenter{cid},exper) && Cluster_List.onMazeAvgFR(cid)>=0.5
        clusterID = Cluster_List.ID{cid};
        sid = find(Session_List.rat==Cluster_List.rat(cid) & Session_List.session==Cluster_List.session(cid));
        try
            

[thisMap,fig,nspks,onmazeAvgFR] = JMGetClusterQuals_4sheet(clusterID, ROOT,Session_List.type{sid},Cluster_List.experimenter{cid},2);
        catch
            disp([clusterID ' is failed!'])
        end
end
end