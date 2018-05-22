clear all
clc
dir=pwd;
cd ..
dir2=pwd;
cd Programs
for s=[2 6 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31]%[2 6 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26]
    Locations=[];
    CueTypes=[];
    AllTimes=[];

    cd(['../../BehavioralData'])
    AllData=dlmread(['BehavioralData.txt'],'\t',1,0);
    for r=1:4
        data=AllData((AllData(:,1)==s & AllData(:,2)==r),3:13);
        tc=1;
        clear tmpLoc context tmpCue tmpTimes
        tmpTimes=[];
        tmpLoc=[];
        tmpCue=[];

        for t=1:60
            if rem(t-1,20)==0
                context=data(t,11);
            end
            if data(t,7)==1 && data(t,9)==2 
                tmpTimes=[tmpTimes round(data(t,2)-8)];
                tmpLoc=[tmpLoc context];
                if data(t,5)==1
                    tmpCue=[tmpCue 3];
                elseif data(t,5)==2
                    tmpCue=[tmpCue 4];
                end
            end
        end
        %If it isn't the first run, add the appropriate offset
        
        Locations=[Locations; tmpLoc'];
        CueTypes=[CueTypes; tmpCue'];
        
        AllTimes=[AllTimes ceil(ceil(tmpTimes)/2)+4+(r-1)*274];
        
    end
    
    %% Change Times into TRs and write output files
    %this is the TR of cue onset.
    AllTimes_TR2=AllTimes+1;
    cd([dir2 '\TimingFiles']);
    dlmwrite(['S' num2str(s) 'Shift_SampleOne.1D'],AllTimes);
    dlmwrite(['S' num2str(s) 'Shift_SampleTwo.1D'],AllTimes_TR2);
    %% Save the stimulus locations (in degrees)
    save(['S' num2str(s) 'Shift_test_conds.mat'])
    cd(dir);
end