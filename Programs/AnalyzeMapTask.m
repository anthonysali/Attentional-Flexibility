%Analyze RSVP Data for Flex01 Study

clear all
clc
sc=1;
for s=[6 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31]
    cd(['../BehavioralData/MappingTask'])
    Acc=[];
    FA=[];
    for r=1:4
       load([num2str(s) '_IEM_ModelTrain_sess01_run0' num2str(r)])
       Acc=[Acc p.accuracy/10];
       FA=[FA p.nFalseAlarms];
    end
    Output(sc,1)=sum(Acc)/40;
    Output(sc,2)=sum(FA);
    sc=sc+1;
    cd ..
end
cd ../OutputFiles
fid = fopen('MapTask.txt','w+');
fprintf(fid,'%s\t %s\t \n', 'Hits', 'FalseAlarm');
fclose(fid);
dlmwrite('MapTask.txt', Output, 'delimiter', '\t', '-append')