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
for bin=1:3
    c=1;
    clear Hold Equal Shift AvgHold AvgEqual AvgShift
    for s=[6 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31]
        load(['sb' num2str(s) 'V1_V4_' num2str(bin) '_Reconstructions.mat']);
        Hold(c,:)=mean(HoldRecon,1);
        Equal(c,:)=mean(EqualRecon,1);
        Shift(c,:)=mean(ShiftRecon,1);
        BoxHold=reshape(mean(HoldRecon,1),res(2),res(1));
        BoxEqual=reshape(mean(EqualRecon,1),res(2),res(1));
        BoxShift=reshape(mean(ShiftRecon,1),res(2),res(1));
        Hold_Left=mean(mean(BoxHold(76:96,46:66))); %81:91,51:61
        Hold_Right=mean(mean(BoxHold(76:96,106:126)));
     
        Equal_Left=mean(mean(BoxEqual(76:96,46:66)));
        Equal_Right=mean(mean(BoxEqual(76:96,106:126)));
     
        Shift_Left=mean(mean(BoxShift(76:96,46:66)));
        Shift_Right=mean(mean(BoxShift(76:96,106:126)));
        
        if Bin==1
            B1_Output(c,:)=[Hold_Left Equal_Left Shift_Left Hold_Right Equal_Right Shift_Right];
        elseif Bin==2
            B2_Output(c,:)=[Hold_Left Equal_Left Shift_Left Hold_Right Equal_Right Shift_Right];
        elseif Bin==3
            B3_Output(c,:)=[Hold_Right Equal_Right Shift_Right Hold_Left Equal_Left Shift_Left];
        end
        

        clearvars -except xx yy s c Hold Equal Shift res sbc bin ax B1_Output B2_Output B3_Output mydir mydir2
        c=c+1;
    end
    
    AvgHold=mean(Hold,1);
    AvgEqual=mean(Equal,1);
    AvgShift=mean(Shift,1);


    
    
    ax(sbc)=subplot(3,3,sbc);
    imagesc(xx,yy,reshape(AvgHold,res(2),res(1)));
    axis equal xy off;
    sbc=sbc+1;
    ax(sbc)=subplot(3,3,sbc);
    imagesc(xx,yy,reshape(AvgEqual,res(2),res(1)));
    axis equal xy off;
    sbc=sbc+1;
    ax(sbc)=subplot(3,3,sbc);
    imagesc(xx,yy,reshape(AvgShift,res(2),res(1)));
    axis equal xy off;
    sbc=sbc+1;

end
match_clim(ax,[-.02 .09])
set(gcf, 'Position', [100, 100, 500, 500])
saveas(gcf,'IEM_Main.png')

%% Save the averages for statistics
cd([mydir2 '\OutputFiles']);
fid = fopen('IEM_Pretrial.txt','w+');
fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t \n', 'MHLeft', 'EQLeft', 'MSLeft', 'MHRight', 'EQRight', 'MSRight');
fclose(fid);
dlmwrite('IEM_Pretrial.txt', B1_Output, 'delimiter', '\t', '-append')

fid = fopen('IEM_Post.txt','w+');
fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t \n', 'MHLeftShift', 'EQLeftShift', 'MSLeftShift', 'MHRightShift', 'EQRightShift', 'MSRightShift', 'MHRightHold', 'EQRightHold', 'MSRightHold', 'MHLeftHold', 'EQLeftHold', 'MSLeftHold' );
fclose(fid);
dlmwrite('IEM_Post.txt', [B2_Output B3_Output], 'delimiter', '\t', '-append')