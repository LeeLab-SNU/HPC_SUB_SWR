Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip = [ROOT.Save '\ripples_mat'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];
ROOT.React = [ROOT.Save '\react_mat'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);

thisRegion = 'SUB';
thisRegion2 = 'SUB';
filt_time = 0;

if filt_time==0, suff = ''; else, suff = ['_' num2str(filt_time) 's']; end

RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis' '.xlsx']);
UnitsTable = readtable([ROOT.Save '\UnitsTable_' thisRegion '_forAnalysis' '.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion2 suff '.xlsx']);



for clRip = 1:size(RipplesTable,1)
    RipID= cell2mat(RipplesTable.ID(clRip));
    id = find(cellfun(Params.cellfind(RipID),(ReactTable.RippleID)));
    thisReactTable = ReactTable(id,:);
    
    unit=[];
    thisUnitTable =[];
    for clUnit=1:size(thisReactTable,1)
        uid = cell2mat(thisReactTable.UnitID(clUnit));
                id = find(cellfun(Params.cellfind(uid),(UnitsTable.ID)));
        if isempty(id), continue; end

        unit(clUnit) = str2double(uid(8:9))*10^2+str2double(uid(11:12));
        thisUnitTable = [thisUnitTable; UnitsTable(id(1),:)];
        thisReactTable.SI(clUnit) = UnitsTable.SI(id(1));
        thisReactTable.PeakBin(clUnit) = UnitsTable.PeakBin(id(1));
    end
    unit(unit==0)=[];
    if ~isempty(unit)
        [unq_unit,ia,ic] = unique(unit','rows');
        RipplesTable.Spikes(clRip) = length(unit);
        RipplesTable.Ensemble_all(clRip) = length(unq_unit);
        RipplesTable.Ensemble_PC(clRip) = sum(thisUnitTable.SI(ia)>=0.5);
        RipplesTable.Ensemble_OnMazePC(clRip) = sum((thisUnitTable.SI(ia)>=0.5) & (thisUnitTable.PeakBin(ia)>1));
        
        ReactTable.Ensemble(id) = RipplesTable.Ensemble_all(clRip);
    end
    
    
end


RipplesTable_Trial = RipplesTable(RipplesTable.context>0,:);

RipplesTable_Ensemble = RipplesTable(RipplesTable.Ensemble_all>=3,:);

writetable(RipplesTable_Ensemble,[ROOT.Save '\RipplesTable_Ensemble_' thisRegion suff '.xlsx'],'writemode','overwrite')
