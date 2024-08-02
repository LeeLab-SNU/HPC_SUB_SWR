Initial_SWRFilter_common;
warning off
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
ROOT.React = [ROOT.Mother '\Processed Data\react_mat'];
ROOT.Save = [ROOT.Mother '\Processed Data'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);

thisRegion = 'CA1';
filt_time = 0;

if filt_time==0, suff = ''; else, suff = ['_' num2str(filt_time) 's']; end

RipplesTable = readtable([ROOT.Save '\RipplesTable_Ensemble_' thisRegion  suff '.xlsx']);
UnitsTable = readtable([ROOT.Save '\UnitsTable_filtered_' thisRegion '.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion suff '.xlsx']);


thisRipples=table;
RipCountTable = table;
%%
for clRip = 1:size(RipplesTable,1)
    RipID = cell2mat(RipplesTable.ID(clRip));
    thisSID = RipID(1:6);
    
    if clRip~=size(RipplesTable,1)
    RipID = cell2mat(RipplesTable.ID(clRip+1));
    thisSID_p = RipID(1:6);
    end
    
    if strcmp(thisSID,thisSID_p) && clRip~=size(RipplesTable,1)
        thisRipples = [thisRipples;RipplesTable(clRip,:)];
    else
        thisRipples = [thisRipples;RipplesTable(clRip,:)];
        RipCountT=table;
        RipCountT.SID = {thisSID};
        RipCountT.Zebra = sum(thisRipples.context==1);
        RipCountT.Pebbles = sum(thisRipples.context==2);
        RipCountT.Bamboo = sum(thisRipples.context==3);
        RipCountT.Mountains = sum(thisRipples.context==4);
        RipCountTable = [RipCountTable;RipCountT];
        thisRipples=table;
    end
end
writetable(RipCountTable,[ROOT.Save '\RipplesCountTable_' thisRegion suff '.xlsx'])