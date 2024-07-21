Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Mother '\Processed Data\ripples_mat\R1'];

thisRegion0 ='SUB_refCA1'
RipplesTable_ori = readtable([ROOT.Rip0 '\RipplesTable_' thisRegion0 '.xlsx']);

RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion0 '_forAnalysis.xlsx']);
RipplesTable(isnan(RipplesTable.DecodingP_all),:)=[];
writetable(RipplesTable,[ROOT.Save '\RipplesTable_' thisRegion0 '_forAnalysis.xlsx'],'WriteMode','replacefile');
% RipplesTable(RipplesTable.DecodingP_all<0.05,:)

RipplesTable  = Rip_SUB2
for rid=1:size(RipplesTable,1)
    rid2 = find(strcmp(RipplesTable.ID(rid), RipplesTable_ori.ID));
    
    if ~isempty(rid2)
    RipplesTable.NormAmp_v3(rid) = RipplesTable_ori.NormAmp_v3(rid2);
    else
        RipplesTable.NormAmp_v3(rid)=nan;
    end
    
end
RipplesTable(isnan(RipplesTable.NormAmp_v3),:)=[];
Rip_SUB2 = RipplesTable

sum(RipplesTable_ori.NormFilteredPeak<1.2)

sum(RipplesTable.NormFilteredPeak<1.2)

sum(RipplesTable.NormAmp_v3<1)

histogram(RipplesTable.NormFilteredPeak)
hold on
histogram(RipplesTable.NormAmp_v3,'binwidth',0.1)

cdfplot(RipplesTable.NormAmp_v3)
set(gca,'fontsize',12,'FontWeight','b')
title('baseline = ITI non-ripple period amp. mean')
xlabel('normalized amplitude')
ylabel('# ripples')