clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%PARAMETERS%%%              %
%%%%%%%%%%%%%%%%              %
ROI='V1_V4';                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
mydir=pwd;
cd ..
mydir2=pwd;
cd IEM/Programs
clearvars -except ROI Bin Exclusion Bin_list mydir mydir2
si=0;
for s=[6 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31]
        si=si+1;
        sub_indices=[20674 20973 20979 20983 20987 20996 20997 21005 21037 21070 21191 21192 21197 21229 21239 21241 21249 21356 21379 21394 21421 21424];
        
        %% Read in the data
        cd ../Data
        map1=dlmread(['s' num2str(sub_indices(si)) '.map1.V1_V4.txt'],' ');
        map2=dlmread(['s' num2str(sub_indices(si)) '.map2.V1_V4.txt'],' ');
        map3=dlmread(['s' num2str(sub_indices(si)) '.map3.V1_V4.txt'],' ');
        map4=dlmread(['s' num2str(sub_indices(si)) '.map4.V1_V4.txt'],' ');
        
        
        attn1=dlmread(['s' num2str(sub_indices(si)) '.attn1.V1_V4.txt'],' ');
        attn2=dlmread(['s' num2str(sub_indices(si)) '.attn2.V1_V4.txt'],' ');
        attn3=dlmread(['s' num2str(sub_indices(si)) '.attn3.V1_V4.txt'],' ');
        attn4=dlmread(['s' num2str(sub_indices(si)) '.attn4.V1_V4.txt'],' ');
        
        load([mydir '\TimingFiles\S' num2str(s) 'All_trn_conds.mat'])
        mb1=dlmread([mydir '\TimingFiles\S' num2str(s) 'SampleOne.1D']);
        mb2=dlmread([mydir '\TimingFiles\S' num2str(s) 'SampleTwo.1D']);
        Alltrn_conds=AllXYDeg;
        map=[map1'; map2'; map3'; map4'];
        firstBrik=map(mb1,:);
        secondBrik=map(mb2,:);
    
    Alltrn=(firstBrik+secondBrik)/2; %gets the average of the two time points.
    RunNums=[ones(37,1); ones(37,1)+1; ones(37,1)+2; ones(37,1)+3];
    Alltrn_conds=[Alltrn_conds RunNums];
    clearvars -except i currentDir Opt Alltrn s Bin ROI sub_indices si Alltrn_conds Exclusion Bin_list Ind mydir mydir2
    
    Output=zeros(37,2);
    nc=1;
    stim_reconstructions=zeros(37*4,171*171);
    for iter=1:4
        trn=Alltrn(Alltrn_conds(:,3)~=iter,:);
        tst=Alltrn(Alltrn_conds(:,3)==iter,:);
        trn_parse=Alltrn_conds(:,3)~=iter;
        trn_conds=Alltrn_conds(trn_parse,1:2);
        tst_parse=Alltrn_conds(:,3)==iter;
        tst_conds=Alltrn_conds(tst_parse,1:2);
        %% Run the IEM
        cd([mydir '\Programs']);
        root = load_root(); 
        addpath([root 'mFiles\']);
        stim_size = 0.9; 
        
        
        res = [171 171];
        
        [xx, yy] = meshgrid(linspace(-17/2,17/2,res(1)),linspace(-17/2,17/2,res(2)));
        xx = reshape(xx,numel(xx),1);yy = reshape(yy,numel(yy),1);
        
        stim_mask = zeros(length(xx),size(trn_conds,1)); %mask size: rows: num of pixels; columns: number of trials
        
        for ii = 1:size(trn_conds,1) %trn_conds is the x and y position of the stimulus -- this is going from 1 to the number of trials
            
            rr = sqrt((xx-trn_conds(ii,1)).^2+(yy-trn_conds(ii,2)).^2);
            stim_mask(rr <= stim_size,ii) = 1; % assume constant stimulus size here - There is a column for each trial here
            
        end
        
        
        rfPtsX = [-1.5 -0.5 0.5 1.5  -2 -1 0 1 2   -2.5 -1.5 -0.5 0.5 1.5 2.5  -3 -2 -1 0 1 2 3  -2.5 -1.5 -0.5 0.5 1.5 2.5  -2 -1 0 1 2  -1.5 -0.5 0.5 1.5];
        rfPtsY = [-3 -3 -3 -3       -2 -2 -2 -2 -2  -1 -1 -1 -1 -1 -1           0 0 0 0 0 0 0    1 1 1 1 1 1                  2 2 2 2 2         3 3 3 3]*(sqrt(3)/2);
        
        rfPtsX=rfPtsX*2; %step size from running script plus extra to make larger
        rfPtsY=rfPtsY*2; %step size from running script plus extra to make larger
        
        rfGridX = reshape(rfPtsX,numel(rfPtsX),1);rfGridY = reshape(rfPtsY,numel(rfPtsY),1);
        rfSize = 1.25*stim_size;   % for the filter shape we use (see below), this size/spacing ratio works well (see Sprague & Serences, 2013)
        

        
        basis_set = nan(size(xx,1),size(rfGridX,1)); % initialize
        
        for bb = 1:size(basis_set,2)
            
            
            basis_set(:,bb) = make2dcos(xx,yy,rfGridX(bb),rfGridY(bb),rfSize*2.5166,7);
        end
        
        trnX = stim_mask.' * basis_set;
        trnX = trnX./max(trnX(:));        % and normalize to 1
 
        
        fprintf('Rank of design matrix: %i\n',rank(trnX));
        
        w = pinv(trnX) * trn;
        
        
        chan_resp = (inv(w*w')*w*tst').';
        
        %% Align the reconstruction for averaging
        %Each trial needs it's own, UNIQUE, basis set because of the
        %jittering
        orig_stim_reconstructions = chan_resp * basis_set.';
        returnDir=pwd;
        cd([mydir2 '\BehavioralData\MappingTask'])
        load([num2str(s) '_IEM_ModelTrain_sess01_run0' num2str(iter) '.mat'])
        [~,index]=sort(p.rndInd(p.targPresent~=1));
        unSort=chan_resp(index,:); %reorder the trials into their original position

        [xx, yy] = meshgrid(linspace(-17/2,17/2,res(1)),linspace(-17/2,17/2,res(2)));
        xx = reshape(xx,numel(xx),1);yy = reshape(yy,numel(yy),1);
        for jc=1:37
            unJitterXDeg=rfGridX-p.jitterDeg(jc,1);
            unJitterYDeg=rfGridY-p.jitterDeg(jc,2);
            [thgrid, rgrid] = cart2pol(unJitterXDeg,unJitterYDeg);
            thgrid = thgrid + deg2rad(p.angOffset);
            [xLoc, yLoc] = pol2cart(thgrid,rgrid);
            new_set = nan(size(xx,1),size(rfGridX,1)); % initialize
            for bb = 1:size(new_set,2)
                new_set(:,bb) = make2dcos(xx,yy,xLoc(bb),yLoc(bb),2*rfSize*2.5166,7);
            end
            stim_reconstructions(jc+37*(iter-1),:)=unSort(jc,:)*new_set.';
        end
        cd(returnDir)
    end
    %% Get the maximum to compare
    first=stim_reconstructions(1:37,:);
    second=stim_reconstructions(38:74,:);
    third=stim_reconstructions(75:111,:);
    fourth=stim_reconstructions(112:148,:);
    z=cat(3,first,second,third,fourth);
    Avg=mean(z,3);
    xx = linspace(-17/2,17/2,171);
    yy = linspace(-17/2,17/2,171);
    for n=1:size(Avg,1)
        [~,position]=max(Avg(n,:));
        getMax=zeros(1,size(Avg,2));
        getMax(position)=1;
        getMax=reshape(getMax,res(2),res(1));
        [y_coord,x_coord]=find(getMax);
        Output(nc,1)=xx(x_coord);
        Output(nc,2)=yy(y_coord);
        nc=nc+1;
    end
    
    %% Save the data
    clearvars -except s ROI sub_indices si Exclusion Bin_list Output Alltrn_conds Avg mydir mydir2
    cd([mydir '\Outputs']);
    save(['sb' num2str(s) ROI 'ParticipantCheck_Avg.mat']);
end
