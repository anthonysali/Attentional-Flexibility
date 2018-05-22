%Generate Timing Scripts for Flex Attention Task
clear all
clc
Nums=[20495 20674 20973 20979 20983 20987 20996 20997 21005 21037 21070 21191 21192 21197 21229 21239 21241 21249 21356 21379 21394 21421 21424];
c=1;
cd ../fMRI_Files
if exist('PE_GLM')~=7
    mkdir('PE_GLM');
else
    rmdir('PE_GLM', 's');
    mkdir('PE_GLM');
end
cd ../BehavioralData
AllData=dlmread(['BehavioralData.txt'],'\t',1,0);
cd ../fMRI_Files
for s=[2 6 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31]
    clearvars -except s c Nums AllData
    load(['PEs/s' num2str(s) '_PE_final.mat'])
    runnum=4;
    r_list=[1 2 3 4];
    cd ..
    for r=1:runnum
        data=AllData((AllData(:,1)==s & AllData(:,2)==r),3:13);
        times=[];
        Hold=[];
        Shift=[];
        Error=[];
        
       if r==1
           run_unsignPE=unsignPE(1:60);
       elseif r==2
           run_unsignPE=unsignPE(61:120);
       elseif r==3
           run_unsignPE=unsignPE(121:180);
       elseif r==4
           run_unsignPE=unsignPE(181:240);
       end
        
        tc=0;
        for t=1:60
            tc=tc+1;
            thetime=round(data(t,2)-8);
            if data(t,7)==1
                times=[times thetime];
                if data(t,9)==1 %hold trial
                    Hold=[Hold thetime];
                elseif data(t,9)==2 %shift trial
                    Shift=[Shift thetime];
                end
            elseif data(t,7)==0
                run_unsignPE(tc)=NaN;
                run_signPE(tc)=NaN;
                Error=[Error thetime];
            end
        end
        run_unsignPE(isnan(run_unsignPE))=[];
        run_signPE(isnan(run_signPE))=[];

        uns_modulator_out=[times' ones(length(times),1)*2.25 run_unsignPE];
        
        hold_out=[Hold' ones(length(Hold),1)*2.25 ones(length(Hold),1)*1];
        shift_out=[Shift' ones(length(Shift),1)*2.25 ones(length(Shift),1)*1];
        error_out=[Error' ones(length(Error),1)*2.25 ones(length(Error),1)*1];

        if r_list(r)==1
            str='1';
        elseif r_list(r)==2
            str='2';
        elseif r_list(r)==3
            str='3';
        elseif r_list(r)==4
            str='4';
        end
        cd fMRI_Files/PE_GLM

        dlmwrite([num2str(Nums(c)) 'Final_UnsignModulator_r' str '.txt'],uns_modulator_out,'delimiter', '\t', '-append');
        dlmwrite([num2str(Nums(c)) 'Final_Hold_r' str '.txt'],hold_out,'delimiter', '\t', '-append');
        dlmwrite([num2str(Nums(c)) 'Final_Shift_r' str '.txt'],shift_out,'delimiter', '\t', '-append');
        dlmwrite([num2str(Nums(c)) 'Final_Error_r' str '.txt'],error_out,'delimiter', '\t', '-append');
        cd ../..
        
    end
    c=c+1;
    cd fMRI_Files
end