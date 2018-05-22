close all
clear all
root = load_root(); %For now, this points to the tutorial directory
addpath([root 'mFiles\']);
cd ../Outputs
load('HexPlotCoords.mat')
AllData=[];


for s=[6 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31]
load(['sb' num2str(s) 'V1_V4ParticipantCheck_Avg.mat'])
AllData=cat(3,AllData,Avg);
end
Group=mean(AllData,3);
res=[171 171];
for c=1:37
ax(c)=subplot('position', [xCord(c)-.03 yCord(c) .1 .1]);
imagesc(reshape(Group(c,:),res(2),res(1)));
axis equal xy off;
end
match_clim(ax,[-.02 .1])
set(gcf, 'Position', [100, 100, 500, 500])
saveas(gcf,'MapCheck.png')