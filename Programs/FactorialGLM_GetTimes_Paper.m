%Generate Timing Scripts for Flex Attention Task
clear all
clc

Nums=[20495 20674 20973 20979 20983 20987 20996 20997 21005 21037 21070 21191 21192 21197 21229 21239 21241 21249 21356 21379 21394 21421 21424];
c=1;
cd ../fMRI_Files
if exist('FactorialGLM')~=7
    mkdir('FactorialGLM');
else
    rmdir('FactorialGLM', 's');
    mkdir('FactorialGLM');
end
cd ../BehavioralData
AllData=dlmread(['BehavioralData.txt'],'\t',1,0);
cd ../fMRI_Files/FactorialGLM
for s=[2 6 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31]
    clearvars -except s c Nums AllData
    for r=1:4
        data=AllData((AllData(:,1)==s & AllData(:,2)==r),3:13);
        Hold_MH=[];
        Shift_MH=[];
        Hold_EQ=[];
        Shift_EQ=[];
        Hold_MS=[];
        Shift_MS=[];
        Error=[];
        
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
            thetime=round(data(t,2)-8);
            if data(t,7)==1
                if data(t,9)==1 %hold trial
                    if CurrentContext==1
                        Hold_MH=[Hold_MH thetime];
                    elseif CurrentContext==2
                        Hold_EQ=[Hold_EQ thetime];
                    elseif CurrentContext==3
                        Hold_MS=[Hold_MS thetime];
                    end
                elseif data(t,9)==2 %shift trial
                    if CurrentContext==1
                        Shift_MH=[Shift_MH thetime];
                    elseif CurrentContext==2
                        Shift_EQ=[Shift_EQ thetime];
                    elseif CurrentContext==3
                        Shift_MS=[Shift_MS thetime];
                    end
                end
            else
                Error=[Error thetime];
            end
        end
        
        Hold_MH_out=[Hold_MH' ones(length(Hold_MH),1)*2.25 ones(length(Hold_MH),1)];
        Hold_EQ_out=[Hold_EQ' ones(length(Hold_EQ),1)*2.25 ones(length(Hold_EQ),1)];
        Hold_MS_out=[Hold_MS' ones(length(Hold_MS),1)*2.25 ones(length(Hold_MS),1)];
        
        Shift_MH_out=[Shift_MH' ones(length(Shift_MH),1)*2.25 ones(length(Shift_MH),1)];
        Shift_EQ_out=[Shift_EQ' ones(length(Shift_EQ),1)*2.25 ones(length(Shift_EQ),1)];
        Shift_MS_out=[Shift_MS' ones(length(Shift_MS),1)*2.25 ones(length(Shift_MS),1)];
        
        Error_out=[Error' ones(length(Error),1)*2.25 ones(length(Error),1)];
        
        if r==1
            str='1';
        elseif r==2
            str='2';
        elseif r==3
            str='3';
        elseif r==4
            str='4';
        end
        
        dlmwrite([num2str(Nums(c)) 'Hold_MH_r' str '.txt'],Hold_MH_out,'delimiter', '\t', '-append');
        dlmwrite([num2str(Nums(c)) 'Hold_EQ_r' str '.txt'],Hold_EQ_out,'delimiter', '\t', '-append');
        dlmwrite([num2str(Nums(c)) 'Hold_MS_r' str '.txt'],Hold_MS_out,'delimiter', '\t', '-append');
        
        dlmwrite([num2str(Nums(c)) 'Shift_MH_r' str '.txt'],Shift_MH_out,'delimiter', '\t', '-append');
        dlmwrite([num2str(Nums(c)) 'Shift_EQ_r' str '.txt'],Shift_EQ_out,'delimiter', '\t', '-append');
        dlmwrite([num2str(Nums(c)) 'Shift_MS_r' str '.txt'],Shift_MS_out,'delimiter', '\t', '-append');
        
        dlmwrite([num2str(Nums(c)) 'Error_r' str '.txt'],Error_out,'delimiter', '\t', '-append');
        
    end
    c=c+1;
end