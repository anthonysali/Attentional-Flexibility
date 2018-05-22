clear all
clc
dir=pwd;
cd ..
dir2=pwd;
cd Programs
for s=[2 6 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31]
    AllXYDeg=[];
    AllTimes=[];

    cd(['../../BehavioralData/MappingTask'])
    for r=1:4
        onsettime=[];
        XYDeg=[];
        load([num2str(s) '_IEM_ModelTrain_sess01_run0' num2str(r) '.mat']);
        tc=1;
        clear XYDeg onsettime
        for t=1:47
            if p.targPresent(t)==0 %First, make sure we exclude the target present trials
                XYDeg(tc,1:2) = p.stimLocsDeg(t,:);
                onsettime(tc) = p.stimStart(t)-p.startExp-8;
                tc=tc+1;
            end
        end
        %If it isn't the first run, add the appropriate offset
        AllXYDeg=[AllXYDeg; XYDeg];
        AllTimes=[AllTimes ceil(ceil(onsettime)/2)+4+(r-1)*172];
    end
    
    %% Change Times into TRs and write output files
    AllTimes_TR2=AllTimes+1;
    cd([dir2 '\TimingFiles']);
    dlmwrite(['S' num2str(s) 'SampleOne.1D'],AllTimes);
    dlmwrite(['S' num2str(s) 'SampleTwo.1D'],AllTimes_TR2);
    %% Save the stimulus locations (in degrees)
    save(['S' num2str(s) 'All_trn_conds.mat'],'AllXYDeg')
    cd(dir);
end