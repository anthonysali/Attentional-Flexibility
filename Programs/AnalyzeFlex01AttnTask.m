%Analyze RSVP Data for Flex01 Study

clear all
clc
sc=1;
trim=0;
total=0;
cd('../BehavioralData')
AllData=dlmread(['BehavioralData.txt'],'\t',1,0);
for s=[2 6 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31]
    
    StableHold=[];
    StableShift=[];
    EqualHold=[];
    EqualShift=[];
    FlexHold=[];
    FlexShift=[];
    SHAcc=0;
    SSAcc=0;
    EHAcc=0;
    ESAcc=0;
    FHAcc=0;
    FSAcc=0;
    
    TotSH=0;
    TotSS=0;
    TotEH=0;
    TotES=0;
    TotFH=0;
    TotFS=0;
    
    sums=0;
    for r=1:4
        data=AllData((AllData(:,1)==s & AllData(:,2)==r),3:13);
        StableContext=data(3,10);
        EqualContext=data(4,10);
        FlexContext=data(5,10);
        for t=1:60
            if rem(t-1,20)==0
                if StableContext==data(t,11)
                    CurrentContext=1;
                elseif EqualContext==data(t,11)
                    CurrentContext=2;
                elseif FlexContext==data(t,11)
                    CurrentContext=3;
                end
            end
            if data(t,7)==1
                if CurrentContext==1
                    if data(t,9)==1
                        StableHold=[StableHold data(t,4)];
                        SHAcc=SHAcc+1;
                        TotSH=TotSH+1;
                    elseif data(t,9)==2
                        StableShift=[StableShift data(t,4)];
                        SSAcc=SSAcc+1;
                        TotSS=TotSS+1;
                    end
                elseif CurrentContext==2
                    if data(t,9)==1
                        EqualHold=[EqualHold data(t,4)];
                        EHAcc=EHAcc+1;
                        TotEH=TotEH+1;
                    elseif data(t,9)==2
                        EqualShift=[EqualShift data(t,4)];
                        ESAcc=ESAcc+1;
                        TotES=TotES+1;
                    end
                elseif CurrentContext==3
                    if data(t,9)==1
                        FlexHold=[FlexHold data(t,4)];
                        FHAcc=FHAcc+1;
                        TotFH=TotFH+1;
                    elseif data(t,9)==2
                        FlexShift=[FlexShift data(t,4)];
                        FSAcc=FSAcc+1;
                        TotFS=TotFS+1;
                    end
                end
                
            else %inaccurate
                if CurrentContext==1
                    if data(t,9)==1
                        TotSH=TotSH+1;
                    elseif data(t,9)==2
                        TotSS=TotSS+1;
                    end
                elseif CurrentContext==2
                    if data(t,9)==1
                        TotEH=TotEH+1;
                    elseif data(t,9)==2
                        TotES=TotES+1;
                    end
                elseif CurrentContext==3
                    if data(t,9)==1
                        TotFH=TotFH+1;
                    elseif data(t,9)==2
                        TotFS=TotFS+1;
                    end
                end
            end
        end
        sums=sums+sum(data(:,7));
    end
    
    [tStableHold, tc1]= RtTrim(StableHold);
    [tStableShift, tc2]= RtTrim(StableShift);
    [tEqualHold, tc3]= RtTrim(EqualHold);
    [tEqualShift, tc4]= RtTrim(EqualShift);
    [tFlexHold, tc5]= RtTrim(FlexHold);
    [tFlexShift, tc6]= RtTrim(FlexShift);
    
    trim=trim+tc1+tc2+tc3+tc4+tc5+tc6;
    total=total+length(StableHold)+length(StableShift)+length(EqualHold)+length(EqualShift)+length(FlexHold)+length(FlexShift);
    
    Output(sc,1:7)=[s mean(tStableHold)*1000 mean(tStableShift)*1000 mean(tEqualHold)*1000 mean(tEqualShift)*1000 mean(tFlexHold)*1000 mean(tFlexShift)*1000];
    Accuracy(sc,1:7)=[s (SHAcc/TotSH)*100 (SSAcc/TotSS)*100 (EHAcc/TotEH)*100 (ESAcc/TotES)*100 (FHAcc/TotFH)*100 (FSAcc/TotFS)*100];
    sc=sc+1;
end
cd ../OutputFiles
fid = fopen('RTs.txt','w+');
fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t \n', 'SubID', 'HoldMH', 'ShiftMH', 'HoldEQ', 'ShiftEQ', 'HoldMS', 'ShiftMS');
fclose(fid);
dlmwrite('RTs.txt', Output, 'delimiter', '\t', 'precision',10, '-append')
fid = fopen('Accuracies.txt','w+');
fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t \n', 'SubID', 'HoldMH', 'ShiftMH', 'HoldEQ', 'ShiftEQ', 'HoldMS', 'ShiftMS');
fclose(fid);
dlmwrite('Accuracies.txt', Accuracy, 'delimiter', '\t', 'precision',10, '-append')