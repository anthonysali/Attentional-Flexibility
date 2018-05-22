% clear all
sc=1;
close all
clc
cd ../BehavioralData
AllData=dlmread(['BehavioralData.txt'],'\t',1,0);
for s=[2 6 8 9 10 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31]
    clearvars -except s sc Output Devs Ps ver1 Coeffs AllData
    rlist=1:4;
    sbjData=zeros(60*4,3);
    
    blocks=[];
    cat_targ=[];
    cat_cues=[];
    for r =rlist
        data=AllData((AllData(:,1)==s & AllData(:,2)==r),3:13);
        
        
        StableContext=data(3,10);
        EqualContext=data(4,10);
        FlexContext=data(5,10);
        for t=[1,21,41]
            if StableContext==data(t,11)
                CurrentContext=1;
            elseif EqualContext==data(t,11)
                CurrentContext=2;
            elseif FlexContext==data(t,11)
                CurrentContext=3;
            end
            blocks=[blocks ones(1,20)*CurrentContext];
        end
        
        sbjData((1:60)+(60*(r-1)),1)=data(:,4);
        sbjData((1:60)+(60*(r-1)),2)=data(:,7);
        sbjData((1:60)+(60*(r-1)),3)=data(:,9)-1;
    end
    
    ctx=blocks'; 
    bl_1=zeros(240,1);
    bl_2=zeros(240,1);
    bl_1(blocks==1)=1;
    bl_2(blocks==3)=1;
    
    
    blocks(blocks==3)=.75;
    blocks(blocks==2)=.5;
    blocks(blocks==1)=.25;
    
    RT = sbjData(:,1);
    acc = sbjData(:,2);
    sw = sbjData(:,3);  % 0 = hold, 1 = shift attention
    
    badTrials=(acc==0);
    badTrials(find(acc(1:end-1)==0)+1)=1;
    idx = ~(badTrials);  % trial index of all valid trials
    idxSw1 = sw==1 & ctx==1 & idx;
    idxRep1= sw==0 & ctx==1 & idx;
    idxSw2 = sw==1 & ctx==3 & idx;
    idxRep2= sw==0 & ctx==3 & idx;
    idxSw3 = sw==1 & ctx==2 & idx;
    idxRep3= sw==0 & ctx==2 & idx;
    
    trNum = size(sbjData,1);
    estRT = zeros(trNum,1);
    pSw   = ones(trNum,1)*.5; % initialize predicted switch at.5
    
    %% exhaustively try all possible alpha
    iteration = 0;
    for alpha = 0.01:0.01:0.99
        
        iteration = iteration+1;
        for n=1:trNum-1
            pSw(n+1) = pSw(n) + alpha*(sw(n)-pSw(n));
        end
        
        tmp_PE = abs(pSw-sw);
        dmSw1  = [tmp_PE(idxSw1)  ones(sum(idxSw1),1) ones(sum(idxSw1),1) zeros(sum(idxSw1),1) ones(sum(idxSw1),1)];
        dmRep1 = [tmp_PE(idxRep1) zeros(sum(idxRep1),1) ones(sum(idxRep1),1) zeros(sum(idxRep1),1) ones(sum(idxRep1),1)];
        dmSw2  = [tmp_PE(idxSw2)  ones(sum(idxSw2),1) zeros(sum(idxSw2),1) ones(sum(idxSw2),1) ones(sum(idxSw2),1)];
        dmRep2 = [tmp_PE(idxRep2) zeros(sum(idxRep2),1) zeros(sum(idxRep2),1) ones(sum(idxRep2),1) ones(sum(idxRep2),1)];
        dmSw3  = [tmp_PE(idxSw3)  ones(sum(idxSw3),1) zeros(sum(idxSw3),1) zeros(sum(idxSw3),1) ones(sum(idxSw3),1)];
        dmRep3 = [tmp_PE(idxRep3) zeros(sum(idxRep3),1) zeros(sum(idxRep3),1) zeros(sum(idxRep3),1) ones(sum(idxRep3),1)];
        y     = [RT(idxSw1);RT(idxRep1);RT(idxSw2);RT(idxRep2);RT(idxSw3);RT(idxRep3)];
        dm    = [dmSw1;dmRep1;dmSw2;dmRep2;dmSw3;dmRep3];
        
        x = pinv(dm)*y;
        yhat = dm * x;
        diff = y-yhat;
        SSE= sum(diff.^2);
        searchLog(iteration,1:2)=[SSE alpha];
        
    end
    bestAlpha = searchLog(searchLog(:,1)==min(searchLog(:,1)),2);
    bestSSE = searchLog(searchLog(:,1)==min(searchLog(:,1)),1);
    
    %% use the best Alpha to obtain the trial-by-trial prediction and PE
    
    for n=1:trNum-1
        pSw(n+1) = pSw(n) + bestAlpha*(sw(n)-pSw(n));
    end
    unsignPE = abs(pSw-sw); % Absolute PE, what you'd enter for subsequent fMRI analysis
    signPE = pSw-sw;
    
    
    RT_s=RT(acc==1);
    PE_s=unsignPE(acc==1);
    sw_s=sw(acc==1);
    bl_1=bl_1(acc==1);
    bl_2=bl_2(acc==1);
    
    %% For across subject correlation between some kind of RT interaction effect.. and Model fit
    X1 = [PE_s sw_s bl_1 bl_2 ];       %% can also enter low-level trial features.. e.g. exact stim repetition?
    y  = RT_s;
    [B,DEV,STATS] = glmfit(X1,y, 'normal');
    sbjModelFit = STATS.t(2); % t-stat... for across sbj Correlation
    
    Output(sc,1:2)=[bestAlpha sbjModelFit];
    sc=sc+1;
    save(['../fMRI_Files/PEs/s' num2str(s) '_PE_final.mat'],'unsignPE')
end
s=[2 6 8 9 10 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31];
cd ../OutputFiles
fid = fopen('Alphas.txt','w+');
fprintf(fid,'%s\t %s\t %s\t \n', 'ID', 'Alpha', 'Fit');
fclose(fid);
dlmwrite('Alphas.txt', [s' Output(:,1:2) ], 'delimiter', '\t', '-append')
