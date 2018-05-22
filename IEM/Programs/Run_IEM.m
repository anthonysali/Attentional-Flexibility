%Run the IEM analysis for the V1->V4 ROI
% clc

%This script will test time points prior to cue onset and following cue
%onset in three bins (see labels below). Pre-cue is always T=-2,0. Post cue
%is always T = 6,8).

%This code is adapted from T. Sprague.

cd ..
mydir=pwd;
cd ..
mydir2=pwd;
cd BehavioralData
AllData=dlmread(['BehavioralData.txt'],'\t',1,0);
cd ../IEM/Programs

for Bin=1:3 %1=Pretrial, 2=Post Shift Trials, 3=Post Hold Trials
    if Bin==1
        Exclusion=0; %0=all, 1=Shift, 2=Hold
    elseif Bin==2
        Exclusion=1;
    elseif Bin==3
        Exclusion=2;
    end
    clearvars -except Bin Exclusion mydir mydir2 %Get rid of everything in memory except for the Bin variable
    si=0; %This is a counter variable used for the indices below. 
    %We'll run the analysis for the following participants:
    for s=[6 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31]
        si=si+1;
        %We need to specify the index of each participant according to the
        %scanner ID number. Each entry here corresponds to one in s.
        sub_indices=[20674 20973 20979 20983 20987 20996 20997 21005 21037 21070 21191 21192 21197 21229 21239 21241 21249 21356 21379 21394 21421 21424];
        
        %% Read in the data. The data have already been saved in text files using AFNI commands for easy reading.
        cd([mydir '\Data']);
        %First, we'll read the data from the checkerboard task. This will
        %be used to train the encoding model
        map1=dlmread(['s' num2str(sub_indices(si)) '.map1.V1_V4.txt'],' ');
        map2=dlmread(['s' num2str(sub_indices(si)) '.map2.V1_V4.txt'],' ');
        map3=dlmread(['s' num2str(sub_indices(si)) '.map3.V1_V4.txt'],' ');
        map4=dlmread(['s' num2str(sub_indices(si)) '.map4.V1_V4.txt'],' ');
        
        %Next, we'll read in the data from the attention task. This is the
        %set that will generate the reconstructions.
        attn1=dlmread(['s' num2str(sub_indices(si)) '.attn1.V1_V4.txt'],' ');
        attn2=dlmread(['s' num2str(sub_indices(si)) '.attn2.V1_V4.txt'],' ');
        attn3=dlmread(['s' num2str(sub_indices(si)) '.attn3.V1_V4.txt'],' ');
        attn4=dlmread(['s' num2str(sub_indices(si)) '.attn4.V1_V4.txt'],' ');
        
        %% Now load in the timing files and trial labels for the training data
        %The following file contains the variable AllXYDeg. This specifies
        %the location of the checkerboard stimulus for each of the training
        %trials.
        load([mydir '\TimingFiles\S' num2str(s) 'All_trn_conds.mat'])
        %Sample one and sample two specify the time of stimulus
        %presentation. These have been shifted so that they represent the
        %TRs aquired as close to 6 seconds and 8 seconds following the
        %actual stimulus onset.
        mb1=dlmread([mydir '\TimingFiles\S' num2str(s) 'SampleOne.1D']);
        mb2=dlmread([mydir '\TimingFiles\S' num2str(s) 'SampleTwo.1D']);
        %concatenate the data in time across the four runs.
        map=[map1'; map2'; map3'; map4'];
        
        %% Censor Motion: 
        %Code in this section will read in the motion censor 
        %files that were used in the GLM and remove any trials that have at 
        %least one time point of the two that are averaged that were marked as outliers
        %There is a row for each timepoint and a column for each volume
        %censored
        
        try
            mot1=dlmread([mydir2 '\fMRI_Files\MotionFiles\sb.' num2str(sub_indices(si)) '.map1.motion']);
        catch
            warning(['Motion File sb.' num2str(sub_indices(si)) '.map1.motion not detected.  Assuming no motion censoring.']);
            mot1=zeros(172,1);
        end
        try
            mot2=dlmread([mydir2 '\fMRI_Files\MotionFiles\sb.' num2str(sub_indices(si)) '.map2.motion']);
        catch
            warning(['Motion File sb.' num2str(sub_indices(si)) '.map2.motion not detected.  Assuming no motion censoring.']);
            mot2=zeros(172,1);
        end
        try
            mot3=dlmread([mydir2 '\fMRI_Files\MotionFiles\sb.' num2str(sub_indices(si)) '.map3.motion']);
        catch
            warning(['Motion File sb.' num2str(sub_indices(si)) '.map3.motion not detected.  Assuming no motion censoring.']);
            mot3=zeros(172,1);
        end
        try
            mot4=dlmread([mydir2 '\fMRI_Files\MotionFiles\sb.' num2str(sub_indices(si)) '.map4.motion']);
        catch
            warning(['Motion File sb.' num2str(sub_indices(si)) '.map4.motion not detected.  Assuming no motion censoring.']);
            mot4=zeros(172,1);
        end
        %Collapse across columns so we have a column vector where 1=censor,
        %0 = do not censor
        motion=[mean(mot1,2); mean(mot2,2); mean(mot3,2); mean(mot4,2)];
        %Get the indices of the volumes that should be censored.
        motI=find(motion);
        
        %Now, see what timepoints in mb1(2) are not to-be-censored
        first_diff = setdiff(mb1,motI);
        second_diff = setdiff(mb2,motI);
        mborig=mb1;
        clear mb1 mb2
        
        %The following code checks to see if both the first and second
        %to-be-averaged samples would survive censoring. If at least one
        %sample does not survive, throw the trial out.
        cc=1;
        for ww=1:length(first_diff)
            if sum(ismember(first_diff(ww)+1,second_diff))~=0
                mb1(cc)=first_diff(ww);
                cc=cc+1;
            end
        end
        cc=1;
        for ww=1:length(second_diff)
            if sum(ismember(second_diff(ww)-1,first_diff))~=0
                mb2(cc)=second_diff(ww);
                cc=cc+1;
            end
        end
        
        %Average the samples and then apply the censoring to the trial labels
        firstBrik=map(mb1,:);
        secondBrik=map(mb2,:);
        trn=(firstBrik+secondBrik)/2; %gets the average of the two time points.
        
        [~,Ind]=setdiff(mborig,mb1);
        trn_conds=AllXYDeg;
        trn_conds(Ind,:)=[];
        
        %% Now read in the test phase data timepoints
        if Exclusion==0
            ab1=dlmread([mydir '\TimingFiles\S' num2str(s) 'Pretrial_SampleOne.1D']);
            ab2=dlmread([mydir '\TimingFiles\S' num2str(s) 'Pretrial_SampleTwo.1D']);
        elseif Exclusion==1
            ab1=dlmread([mydir '\TimingFiles\S' num2str(s) 'Shift_SampleOne.1D']);
            ab2=dlmread([mydir '\TimingFiles\S' num2str(s) 'Shift_SampleTwo.1D']);
        elseif Exclusion==2
            ab1=dlmread([mydir '\TimingFiles\S' num2str(s) 'Hold_SampleOne.1D']);
            ab2=dlmread([mydir '\TimingFiles\S' num2str(s) 'Hold_SampleTwo.1D']);
        end
        attn=[attn1'; attn2'; attn3'; attn4'];
        
        %% Censor Motion (see description above)
        
        try
            mot1=dlmread([mydir2 '\fMRI_Files\MotionFiles\sb.' num2str(sub_indices(si)) '.test1.motion']);
        catch
            warning(['Motion File sb.' num2str(sub_indices(si)) '.test1.motion not detected.  Assuming no motion censoring.']);
            mot1=zeros(274,1);
        end
        try
            mot2=dlmread([mydir2 '\fMRI_Files\MotionFiles\sb.' num2str(sub_indices(si)) '.test2.motion']);
        catch
            warning(['Motion File sb.' num2str(sub_indices(si)) '.test2.motion not detected.  Assuming no motion censoring.']);
            mot2=zeros(274,1);
        end
        try
            mot3=dlmread([mydir2 '\fMRI_Files\MotionFiles\sb.' num2str(sub_indices(si)) '.test3.motion']);
        catch
            warning(['Motion File sb.' num2str(sub_indices(si)) '.test3.motion not detected.  Assuming no motion censoring.']);
            mot3=zeros(274,1);
        end
        try
            mot4=dlmread([mydir2 '\fMRI_Files\MotionFiles\sb.' num2str(sub_indices(si)) '.test4.motion']);
        catch
            warning(['Motion File sb.' num2str(sub_indices(si)) '.test4.motion not detected.  Assuming no motion censoring.']);
            mot4=zeros(274,1);
        end
        motion=[mean(mot1,2); mean(mot2,2); mean(mot3,2); mean(mot4,2)];
        motI=find(motion);
        
        first_diff = setdiff(ab1,motI);
        second_diff = setdiff(ab2,motI);
        aborig=ab1;
        clear ab1 ab2
        cc=1;
        for ww=1:length(first_diff)
            if sum(ismember(first_diff(ww)+1,second_diff))~=0
                ab1(cc)=first_diff(ww);
                cc=cc+1;
            end
        end
        cc=1;
        for ww=1:length(second_diff)
            if sum(ismember(second_diff(ww)-1,first_diff))~=0
                ab2(cc)=second_diff(ww);
                cc=cc+1;
            end
        end
        
        [~,Ind]=setdiff(aborig,ab1);
        
        % Average the samples:
        firstBrikA=attn(ab1,:);
        secondBrikA=attn(ab2,:);
        tst=(firstBrikA+secondBrikA)/2; %gets the average of the two time points.
        
        %% Load the test phase trial labels
        if Exclusion==0
            load([mydir '\TimingFiles\S' num2str(s) 'Pretrial_test_conds.mat'], 'CueTypes')
            load([mydir '\TimingFiles\S' num2str(s) 'Pretrial_test_conds.mat'], 'Locations')
        elseif Exclusion==1
            load([mydir '\TimingFiles\S' num2str(s) 'Shift_test_conds.mat'], 'CueTypes')
            load([mydir '\TimingFiles\S' num2str(s) 'Shift_test_conds.mat'], 'Locations')
        elseif Exclusion==2
            load([mydir '\TimingFiles\S' num2str(s) 'Hold_test_conds.mat'], 'CueTypes')
            load([mydir '\TimingFiles\S' num2str(s) 'Hold_test_conds.mat'], 'Locations')
        end
        
        %combine the labels into a single variable, tst_conds, and then
        %apply the censoring computed above to these labels.
        tst_conds=[Locations CueTypes];
        tst_conds(Ind,:)=[];

        %get rid of anything that is no longer needed.
        clearvars -except trn tst trn_conds tst_conds s Bin sub_indices si Exclusion mydir mydir2
        
        %% Run the IEM
        cd([mydir '\Programs']);
        root = load_root(); %The next two lines add the support functions to the matlab search path
        addpath([root 'mFiles\']);
        stim_size = 0.9; % This is the radius of the training stimulus in degrees of visual angle.
        
        
        res = [171 171]; %The resolution at which we will compute the reconstructions.
        %The aperature in the experiment was ~8.35 degrees. Here, we will
        %set the bounding box at 8.5 degrees.
        
        %In xx and yy, each row represents a pixel in the reconstruction.
        [xx, yy] = meshgrid(linspace(-17/2,17/2,res(1)),linspace(-17/2,17/2,res(2)));
        xx = reshape(xx,numel(xx),1);yy = reshape(yy,numel(yy),1);
        
        stim_mask = zeros(length(xx),size(trn_conds,1)); %mask size: rows: num of pixels; columns: number of trials
        
        for ii = 1:size(trn_conds,1) %trn_conds is the x and y position of the stimulus
            
            rr = sqrt((xx-trn_conds(ii,1)).^2+(yy-trn_conds(ii,2)).^2);
            stim_mask(rr <= stim_size,ii) = 1; % assume constant stimulus size here - There is a column for each trial here
            
        end
        
        %stim_mask is variable that is 1s for stimulus present at that
        %location and 0s for stimulus absent at that location.
        
        %Now we'll generate the basis functions of the model. Start by
        %defining a hexagonal grid as in the task file.
        rfPtsX = [-1.5 -0.5 0.5 1.5  -2 -1 0 1 2   -2.5 -1.5 -0.5 0.5 1.5 2.5  -3 -2 -1 0 1 2 3  -2.5 -1.5 -0.5 0.5 1.5 2.5  -2 -1 0 1 2  -1.5 -0.5 0.5 1.5];
        rfPtsY = [-3 -3 -3 -3       -2 -2 -2 -2 -2  -1 -1 -1 -1 -1 -1           0 0 0 0 0 0 0    1 1 1 1 1 1                  2 2 2 2 2         3 3 3 3]*(sqrt(3)/2);
        
        rfPtsX=rfPtsX*2; %We'll space the grid so that it is slightly larger than that used to acquire the data.
        rfPtsY=rfPtsY*2; 
        
        rfGridX = reshape(rfPtsX,numel(rfPtsX),1);rfGridY = reshape(rfPtsY,numel(rfPtsY),1);
        rfSize = 1.25*stim_size;   % This parameter determines the size of each basis function. (see Sprague & Serences, 2013)
        
        basis_set = nan(size(xx,1),size(rfGridX,1)); % initialize
        
        for bb = 1:size(basis_set,2)
            %As in earlier studies, we use gaussian-like basis functions
            %which reach zero at 2.5166x the FWHM of the function. The FWHM
            %will be equal to rfSize above, so the size parameter of our
            %function is 2.5166*rfSize. The function make2dcos generates
            %each individual basis function.
            
            basis_set(:,bb) = make2dcos(xx,yy,rfGridX(bb),rfGridY(bb),rfSize*2.5166,7);
        end
        
        %Now we will fit the model. We start by computing the overlap of
        %each stimulus mask onto the basis functions, which will serve as
        %the design matrix.
        
        trnX = stim_mask.' * basis_set;
        trnX = trnX./max(trnX(:));        % and normalize to 1
        
        %Next determine if the design matrix is full rank (this should be
        %37).
        
        fprintf('Rank of design matrix: %i\n',rank(trnX));
        
        % So we've established that the design matrix is full-rank (the rank is
        % equal to the number of columns). Let's estimate the channel weights.
        
        % this is just pinv(trnX) * trn; - but I wrote out the full math for
        % completness
        w = pinv(trnX) * trn;
        %This could also be written w = inv(trnX.'*trnX)*trnX.' * trn; or as: w = trnX\trn;
        
        % w is now n_channels x n_voxels. Next, invert the weight matrix 
        
        
        chan_resp = (inv(w*w')*w*tst').';
        
        
        %% Rotating the Reconstructions
        
        % We need to adjust the basis set so that it is reconstructed at the
        % horizontal locations
        
        %These variables are used below to separate out shift and hold
        %trials for the post-cue analyses.
        if Exclusion==0
            lookbelow=99;
            lookabove=0;
        elseif Exclusion==1
            lookbelow=99;
            lookabove=2;
        elseif Exclusion==2
            lookbelow=3;
            lookabove=0;
        end
        
        %Context 1 (/), Left:
        clear adj_set thisidx context_chan_resp
        adj_set = RotateBasis(rfPtsX,rfPtsY,-60,rfSize,xx,yy);
        thisidx = tst_conds(:,1)==1 & (tst_conds(:,2)==4 | tst_conds(:,2)==1) & tst_conds(:,2) > lookabove & tst_conds(:,2) < lookbelow;
        context_chan_resp = chan_resp(thisidx,:);
        stim_reconstructions_1_L = context_chan_resp * adj_set.';
        
        %Context 1 (/), Right:
        clear adj_set thisidx context_chan_resp
        adj_set = RotateBasis(rfPtsX,rfPtsY,120,rfSize,xx,yy);
        thisidx = tst_conds(:,1)==1 & (tst_conds(:,2)==3 | tst_conds(:,2)==2) & tst_conds(:,2) > lookabove & tst_conds(:,2) < lookbelow;;
        context_chan_resp = chan_resp(thisidx,:);
        stim_reconstructions_1_R = context_chan_resp * adj_set.';
        
        %Context 2(--), Left:
        clear adj_set thisidx context_chan_resp
        adj_set = RotateBasis(rfPtsX,rfPtsY,0,rfSize,xx,yy);
        thisidx = tst_conds(:,1)==2 & (tst_conds(:,2)==4 | tst_conds(:,2)==1) & tst_conds(:,2) > lookabove & tst_conds(:,2) < lookbelow;;
        context_chan_resp = chan_resp(thisidx,:);
        stim_reconstructions_2_L = context_chan_resp * adj_set.';
        
        %Context 2 (--), Right:
        clear adj_set thisidx context_chan_resp
        adj_set = RotateBasis(rfPtsX,rfPtsY,-180,rfSize,xx,yy);
        thisidx = tst_conds(:,1)==2 & (tst_conds(:,2)==3 | tst_conds(:,2)==2) & tst_conds(:,2) > lookabove & tst_conds(:,2) < lookbelow;;
        context_chan_resp = chan_resp(thisidx,:);
        stim_reconstructions_2_R = context_chan_resp * adj_set.';
        
        %Context 3(\), Left:
        clear adj_set thisidx context_chan_resp
        adj_set = RotateBasis(rfPtsX,rfPtsY,60,rfSize,xx,yy);
        thisidx = tst_conds(:,1)==3 & (tst_conds(:,2)==4 | tst_conds(:,2)==1) & tst_conds(:,2) > lookabove & tst_conds(:,2) < lookbelow;;
        context_chan_resp = chan_resp(thisidx,:);
        stim_reconstructions_3_L = context_chan_resp * adj_set.';
        
        %Context 3 (\), Right:
        clear adj_set thisidx context_chan_resp
        adj_set = RotateBasis(rfPtsX,rfPtsY,-120,rfSize,xx,yy);
        thisidx = tst_conds(:,1)==3 & (tst_conds(:,2)==3 | tst_conds(:,2)==2) & tst_conds(:,2) > lookabove & tst_conds(:,2) < lookbelow;;
        context_chan_resp = chan_resp(thisidx,:);
        stim_reconstructions_3_R = context_chan_resp * adj_set.';
        
        C1 = [stim_reconstructions_1_L; stim_reconstructions_1_R];
        C2 = [stim_reconstructions_2_L; stim_reconstructions_2_R];
        C3 = [stim_reconstructions_3_L; stim_reconstructions_3_R];
        
        %% Assign the Reconstructions
        
        %Read in the behavioral data to get the context labels.
        data=AllData((AllData(:,1)==s & AllData(:,2)==1),3:13);
        Hold=data(3,10);
        Equal=data(4,10);
        Shift=data(5,10);
        
        if Hold==1
            HoldRecon=C1;
        elseif Hold==2
            HoldRecon=C2;
        elseif Hold==3
            HoldRecon=C3;
        end
        
        if Equal==1
            EqualRecon=C1;
        elseif Equal==2
            EqualRecon=C2;
        elseif Equal==3
            EqualRecon=C3;
        end
        
        if Shift==1
            ShiftRecon=C1;
        elseif Shift==2
            ShiftRecon=C2;
        elseif Shift==3
            ShiftRecon=C3;
        end
        
        %% Save the data
        clearvars -except s HoldRecon EqualRecon ShiftRecon ROI Bin sub_indices si Exclusion seed ROI_list mydir mydir2
        cd([mydir '\Outputs']);
        save(['sb' num2str(s) 'V1_V4_' num2str(Bin) '_Reconstructions.mat']);
    end
end