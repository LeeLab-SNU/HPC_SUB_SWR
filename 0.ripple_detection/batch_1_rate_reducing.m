Initial_SWRFilter_common;

thisRID = 069; %551 5-11, 564 4-11 578 4-12 588 4-11
raw_dir = ['Y:\EPhysRawData\CA1,CA3,V2L Recording_VR_JS'];
save_dir = 'F:\EPhysRawData\RawData';
for thisSID = 1:1
RateReducer(thisRID, thisSID, raw_dir, save_dir,'AG')
end

