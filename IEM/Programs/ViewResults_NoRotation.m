clear all
% close all
% clc
cd ..
mydir=pwd;
cd ..
mydir2=pwd;
cd IEM/Programs
root = load_root(); %For now, this points to the tutorial directory
addpath([root 'mFiles\']);
figure

xx = linspace(-17/2,17/2,171);
yy = linspace(-17/2,17/2,171);
res=[171 171];

sbc=1;
cd([mydir '\Outputs']);
for bin=3:3
    c=1;
    clear Hold Equal Shift AvgHold AvgEqual AvgShift
    for s=[6 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31]
        load(['sb' num2str(s) 'V1_V4_' num2str(bin) '_Reconstructions_NoRotation.mat']);
        aC1L(c,:)=mean(C1L,1);
        aC1R(c,:)=mean(C1R,1);
        aC2L(c,:)=mean(C2L,1);
        aC2R(c,:)=mean(C2R,1);
        aC3L(c,:)=mean(C3L,1);
        aC3R(c,:)=mean(C3R,1);
        clearvars -except xx yy s c aC1L aC1R aC2L aC2R aC3L aC3R res sbc bin ax B1_Output B2_Output B3_Output
        c=c+1;
    end
    
    AvgC1L=mean(aC1L,1);
    AvgC1R=mean(aC1R,1);
    AvgC2L=mean(aC2L,1);
    AvgC2R=mean(aC2R,1);
    AvgC3L=mean(aC3L,1);
    AvgC3R=mean(aC3R,1);


    
    
    ax(sbc)=subplot(3,2,sbc);
    imagesc(xx,yy,reshape(AvgC3L,res(2),res(1)));
    axis equal xy off;
    sbc=sbc+1;
    ax(sbc)=subplot(3,2,sbc);
    imagesc(xx,yy,reshape(AvgC1R,res(2),res(1)));
    axis equal xy off;
    sbc=sbc+1;
    ax(sbc)=subplot(3,2,sbc);
    imagesc(xx,yy,reshape(AvgC2L,res(2),res(1)));
    axis equal xy off;
    sbc=sbc+1;
    ax(sbc)=subplot(3,2,sbc);
    imagesc(xx,yy,reshape(AvgC2R,res(2),res(1)));
    axis equal xy off;
    sbc=sbc+1;
    ax(sbc)=subplot(3,2,sbc);
    imagesc(xx,yy,reshape(AvgC1L,res(2),res(1)));
    axis equal xy off;
    sbc=sbc+1;
    ax(sbc)=subplot(3,2,sbc);
    imagesc(xx,yy,reshape(AvgC3R,res(2),res(1)));
    axis equal xy off;
    sbc=sbc+1;

end
match_clim(ax,[-.02 .1])
set(gcf, 'Position', [100, 100, 500, 500])
saveas(gcf,'IEM_NoRotation.png')