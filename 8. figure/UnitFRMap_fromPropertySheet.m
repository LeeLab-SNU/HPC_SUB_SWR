    Clist_f = {'#FF6A68','#00B455','#00A2FF','#FA7344','#00B7A1'};

figure; hold on
for s=1:1
    plot([1:size(FRMap_A,2)],smooth(FRMap_A(1,:,s)),'linewidth',2,'color',[.5 .5 .5])
plot([1:size(FRMap_A,2)],zeros(size(FRMap_A,2),1),'color','k');
    UID = Unitsa{s};
Fidx = find(strncmp(UID,Units(:,1),2));

for f=1:size(Fidx,1)
    fields = smooth(FRMap_B(1,:,Fidx(f)));
% fields(fields==0)=nan;
    plot([1:size(FRMap_A,2)],fields,'linewidth',2,'color',hex2rgb(Clist_f{f}))
    ylim([0 16])
end
end

figure;
imagesc(smooth(FRMap_A(1,:,1))')
colormap jet