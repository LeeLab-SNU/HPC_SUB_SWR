
Initial_SWRFilter_common;

% mother_root = 'H:\CA1&SUB_SCSM\ephys_analysis';
mother_root = ['Y:\EPhysRawData\HPC-SWR\EPhysRawData\RawData'];

Cluster_List = readtable([ROOT.Info '\ClusterList_SWR_SUB.xlsx']);


%
cd(mother_root);
% [~,inputCSV] = xlsread('D:\HPC-LFP project\Information Sheet\ClusterList.csv');

stCellRun = 1;

for cellRUN = 1 : size(Cluster_List,1)
%     if ~strcmp(Cluster_List.experimenter{cellRUN},'JS') 

    
    thisCLUSTER = Cluster_List.ID{cellRUN};
           
    createParsedSpike(mother_root, thisCLUSTER);
    disp([thisCLUSTER ' has processed.']);

end