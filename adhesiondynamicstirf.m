function [Data, ShiftSF, Stack]=adhesiondynamicstirf(imFile,regtf,thresrg, sizeFilt, frameFilt, maskType)

%This function computes the time delay (as well as other metrics) between
%the accumulation dynamics of two proteins at focal adhesions. Input data
%are
% 2-channels time-series of cultured cells acquired with a Nikon TIRF
% microscope, as described in Saurav et al, Cooper lab, Fred Hutch Cancer
% Center.

%INPUT          -imFile: string name of a .nd2 or exported .tif timelapse
%                movie file with 2 channels (R and G, ch#1 -mCherry- and
%                ch#2 -GFP-, respectively), corresponding to 2 proteins of
%                interest.
%
%               -regtf: boolean for applying image registration (1) or not
%               (0)
%
%               -thresrg: channel(s) used for adhesion cluster
%               thresholding:
%                       'r': red (first) channel only 'g': green (second)
%                       channel only 'rg': union of both channels
%
%               -sizeFilt: scalar defining the minimum size (in pixels) of
%               clusters (default 20)
%
%               -frameFilt: scalar defining the minimum number of frames a
%                           cluster should persist (default 3)
%
%               -Masktf: switch scalar defining mask behavior
%                            0 - inside mask excluded 1 -  outside mask
%                            excluded
%                            
%
%OUTPUT:        -Data: a structure array with 3 fields:
%                   --Data.FA a nCluster-by-14 table aggregating 14 metrics
%                   for each tracked adhesion cluster:
%
% FA_Mask        | x-by-y-by-nFrames binary mask of a given tracked cluster
% Starting Frame | index of frame where tracked cluster is first detected
% Intensities    | nFrames-by-2 array containing average mean intensities
%                  of channel 1 (col1) and channel 2 (col2) over time
% Shift          | time delay (in seconds) between ch1 and ch2 protein
%                  (intensities) accumulation
% R/Gwidth       | width of intensity timeseries peak (if any) in seconds
%                  for ch1/ch2
% aRateR/G       | rate of accumulation (sec-1) %dRateR/G      | rate of
% disappearance (sec-1) strengthR/G    | ratio of max mean intensity at
% cluster over max global
%                  intensity
% ShapeR/G       | overall shape of the timeseries (bell, upward, downward,
%                  u-shape)
%
%                   --Binary: a binary XYT stack containing the FA masks at
%                   each
%                     timepoint
%
%                   --Thresholds: a 2-by-1 vector containing the threshold
%               values of red and green (#1 and #2) channels.
%
%               -ShiftSF: n-by-2 matrix containing the time shift (first
%               column) and frame start (second column) for each tracked
%               protein cluster
%
%               -Stack: x-by-y-by-nFrames-by-nChannels containing the image
%               data
%
%
% Example:  Analyze 'DemoFATIRF.tif' after  registering the  stack, using
% both red and green channels to define clusters, size filtered of 20
% pixels and track length of at least 3  frames:
%
% [Data, ShiftSF, Stack]=adhesiondynamicstirf('DEMOFATIRF.tif',1,'rg',...
% 20, 3, 1);
%
% First a figure will pop up showing the last frame of the first channel,
% asking the user to define with freehand tools an ROI.  The outside of the
% ROI wqill be excluded from the analysis (last argument being '1'). To
% exclude the inside of the ROI, the last argument should be set to '0'
% (this option will show the first frame of the first 'R' channel).
% 
% After the analysis is completed, the script will create 2 figures: a
% montage with 2 panels, with the footprint of each tracked FA highlighted
% and colorcoded using time a scatter plot showing the time shift
% distribution
%
% Created by Julien Dubrulle, 2022.



%% load and read imagestack

IS = bfopen(imFile);
R = cat(3,IS{1,1}{1:2:end,1}); % first channel (traditionally mCherry, hence 'R')
G = cat(3,IS{1,1}{2:2:end,1}); % second channel (traditionally GFP, hence 'G')


%% Stack registration

if regtf
    Stack = cat(4,R,G);
    StackReg = ImRegStack(Stack); %extra function
    R = StackReg(:,:,:,1);
    G = StackReg(:,:,:,2);    
end

Stack = cat(4,R,G);

%% Channel preprocessing
nT = size(R,3);
framerate = 30 * 60 / nT; % each movie lasts 30 minutes

Rc = zeros(size(R),'like',R);
Gc = zeros(size(G),'like',G);

for i = 1:size(R,3)
    Rc(:,:,i) = imtophat(medfilt2(R(:,:,i)),strel('disk',5));
    Gc(:,:,i) = imtophat(medfilt2(G(:,:,i)),strel('disk',5));
end

if maskType == 0
figure, imshow(R(:,:,1),[]);
roi = drawfreehand;
MASKCELL = createMask(roi);
close;
elseif maskType == 1
    figure, imshow(R(:,:,end),[]);
    roi = drawfreehand;
    MASKCELL = createMask(roi);
    MASKCELL =~ MASKCELL;
    close;
end

%% Focal Adhesion segmentation using the 1/4 intensity values of
% the 50st brightest pixel (from each channel) as threshold

RWS = false(size(R));
RW = false(size(R));
GW = false(size(R));
TR = sort(Rc(:),'descend');
zr = TR(50);
kr = zr * 0.25;

TG = sort(Gc(:),'descend');
zg = TG(50);
kg = zg * 0.25;
tvalues = [kr kg];
for i = 1:nT
    Rw = Rc(:,:,i) > kr;
    Gw = Gc(:,:,i) > kg;
    RW(:,:,i) = bwareaopen(Rw,sizeFilt);
    GW(:,:,i) = bwareaopen(Gw,sizeFilt);
    if strcmp(thresrg,'rg')
        RWS(:,:,i) = RW(:,:,i) | GW(:,:,i);
    elseif strcmp(thresrg,'r')
        RWS(:,:,i) = RW(:,:,i);
    elseif strcmp(thresrg,'g')
        RWS(:,:,i) = GW(:,:,i);
    else
        error('wrong channels selected')
    end
    RWStemp = RWS(:,:,i);
    DW = imreconstruct(MASKCELL,RWStemp);
    RWStemp(DW) = false;
    RWS(:,:,i) = RWStemp;
end


%% FA tracking

%%% overlap object tracking 3D %%%%%%
% NB: we perform rough tracking by assuming that footprints of 
% individual moving protein clusters overlap from one frame to the  other.
% applying bwlabeln on a XYT stack create one label per  such cluster.
% One caveat is that split clusters  are defined as part of the same FA
% cluster. This caveat will be addressed in the next iteration.


progressText(0,'FA tracking and signal intensity extraction'); %extra function

[L, nl] = bwlabeln(RWS,6);

for i = 1:nl
    M = L == i;
    if sum(max(M,[],[1 2])) < frameFilt
        L(L  == i) = 0;
    end
end

nlu = unique(L);
nlu(1) = [];
nnt = numel(nlu);

Tsmi = cell(nT,nnt);
AllCent = cell(nnt,1);

for i = 1:nnt
    M = L == nlu(i);
    Mm = max(M,[],[1 2]);
    ls = find(Mm,1,'Last');
    for j = 1:nT
        Mi = M(:,:,j);
        if any(Mi(:))
            SMI = zeros(nT,2);
            for k = 1:nT
                Ri = Rc(:,:,k);
                SMI(k,1) = mean(double(Ri(Mi)));
                Gi = Gc(:,:,k);
                SMI(k,2) = mean(double(Gi(Mi)));
            end
            Tsmi{j,i} = SMI;
            if j == ls
                SL = regionprops(double(Mi),'Centroid');
                AllCent{i} = SL.Centroid;
            end
        end
    end
    
    progressText(i / nnt);
end

T = cell(nnt,3);
for i = 1:nnt
    Tt = Tsmi(:,i);
    tfe = cellfun(@isempty,Tt);
    k = find(~tfe,1);
    T{i,1} = L == nlu(i);
    T{i,2} = k;
    T{i,3} = reshape(vertcat(Tt{:}),nT,sum(~tfe),2);
end


%% Computing time delay (difference between time to reach 50 percent
% of normalized intensity and slope

maxintrg = zeros(nnt,2);
AllHGC = cell(nnt,1);
AllWIDTH = cell(nnt,1);
ALLARATE = cell(nnt,1);
ALLDRATE = cell(nnt,1);
ALLSHAPE = cell(nnt,1);
for i = 1:nnt
    Ttt = T{i,3};
    maxintrg(i,:) = max(squeeze(Ttt(:,1,:)),[],1);
    tf=Ttt(1,:,1) >= kr | Ttt(1,:,2) >= kg;
    Ttt(:,tf,:) = [];
    if ~isempty(Ttt)
        Tm = squeeze(mean(Ttt,2));
        T{i,3} = Tm;
        allahg = NaN(1,1);
        araterg = NaN(1,2);
        draterg = NaN(1,2);
        widthrg = NaN(1,2);
        
        %%% SMOOTHING DATA %%%
        
        tr = rescale(smoothdata(Tm(:,1),'movmean'));
        tg = rescale(smoothdata(Tm(:,2),'movmean'));

        %%%%%%%%%%%%%%%%%%%%%%%
        
        [hr, sr, wr] = findThalf2(tr,0.5);
        [hg, sg, wg] = findThalf2(tg,0.5);

        if sr == "bell" && sg == "bell"
            allahg(1) = hg(1) - hr(1);
        elseif (sr =="bell" && sg == "up") ||...
                (sr == "up" && sg == "bell") ||...
                (sr == "up" && sg == "up")
            allahg(1) = hg(1) - hr(1);
        end
        
        widthrg(1,1) = wr;
        widthrg(1,2) = wg;
        
        if sr == "bell" || sr == "up"
            [rrlow, rvalow] = findThalf(tr,0.2,'ascend');
            [rrhigh, rvahigh] = findThalf(tr,0.8,'ascend');
            araterg(1,1) = (rvahigh - rvalow) / (rrhigh - rrlow);
        else
            araterg(1,1) = NaN;
        end
        if sg== "bell" || sg =="up"
            [rglow, gvalow] = findThalf(tg,0.2,'ascend');
            [rghigh, gvahigh] = findThalf(tg,0.8,'ascend');
            araterg(1,2) = (gvahigh - gvalow) / (rghigh - rglow);
        else
            araterg(1,2) = NaN;
        end
        
        if sr == "bell" || sr == "down"
            [drrlow, drvalow] = findThalf(tr,0.2,'descend');
            [drrhigh, drvahigh] = findThalf(tr,0.8,'descend');
            draterg(1,1) = (drvahigh - drvalow) / (drrhigh - drrlow);
        else
            draterg(1,1) = NaN;
        end
        if sg == "bell" || sg == "down"
            [drglow, dgvalow] = findThalf(tg,0.2,'descend');
            [drghigh, dgvahigh] = findThalf(tg,0.8,'descend');
            draterg(1,2) = (dgvahigh - dgvalow) / (drghigh - drglow);
        else
            draterg(1,2) = NaN;
        end
        shaperg = [sr sg];
    else
        Tm = [];
        T{i,3} = Tm;
        allahg = NaN(1,1);
        araterg = NaN(1,2);
        draterg = NaN(1,2);
        widthrg = NaN(1,2);
        shaperg = NaN(1,2);
    end
    
    AllHGC{i,1} = allahg * framerate;
    AllWIDTH{i,1} = widthrg * framerate;
    ALLARATE{i,1} = araterg / framerate;
    ALLDRATE{i,1} = draterg / framerate;
    ALLSHAPE{i,1} = shaperg;
end

maxintrg = maxintrg ./ double([zr zg]);

%% Convert to table for export

tf = cellfun(@isempty,T(:,3));
T(tf,:) = [];
AllHGC(tf,:) = [];
AllWIDTH(tf,:) = [];
ALLARATE(tf,:) = [];
ALLDRATE(tf,:) = [];
ALLSHAPE(tf,:) = [];
maxintrg(tf,:) = [];
AllCent(tf) = [];

T = cell2table([T AllHGC num2cell(vertcat(AllWIDTH{:}))...
    num2cell(vertcat(ALLARATE{:})) num2cell(vertcat(ALLDRATE{:}))...
    num2cell(maxintrg) num2cell(vertcat(ALLSHAPE{:}))],...
    'VariableNames',{'FA_Mask','StartingFrame','Intensities',...
    'Shift','Rwidth','Gwidth','aRateR','aRateG','dRateR','dRateG',...
    'strengthR','strengthG','ShapeR','ShapeG'});

ft = size(T,1);
FF = T.FA_Mask;
MASK = false(size(Rc));
for i = 1:ft
    BW = FF{i};
    MASK(BW) = true;
end

RWS = MASK;
Data.FA = T;
Data.Binary = RWS;
Data.Thresholds = tvalues;

%% Rendering segmentation and tracking

B = cell(nT,1);
for i = 1:nT
    B{i} = bwboundaries(RWS(:,:,i));
end

col = colormap(cool(nT));

imshow(RWS(:,:,end),'InitialMagnification',300);
hold on;
for j = 1:nT
    if ~isempty(B{j})
        BB = B{j};
        for k = 1:length(BB)
            bb = BB{k};
            plot(gca,bb(:,2),bb(:,1),'Color',col(j,:),'Linewidth',1);
        end
    end
end

for i = 1:ft
    text(AllCent{i,1}(1,1),AllCent{i,1}(1,2),int2str(i),'Color','y','FontSize',16);
end
BR = bwboundaries(MASKCELL);

for j = 1:length(BR)
    bbr = BR{j};
    plot(gca,bbr(:,2),bbr(:,1),'-r','Linewidth',2);
end

RGB = print('-RGBImage');
[a, b, ~] = size(RGB);
Rf = imresize(Rc(:,:,end),[a b]);
Gf = imresize(Gc(:,:,end),[a b]);
RGBf = [RGB cat(3,uint8(imadjust(Rf) ./ 2^8),uint8(imadjust(Gf) ./ 2^8),zeros(a,b,'like',RGB))];
figure, imshow(RGBf);

%% Plotting
Shift = [Data.FA.Shift];
SF = [Data.FA.StartingFrame];
ShiftSF = [Shift SF];


figure; swarmchart(ones(numel(Shift),1),Shift,72,SF,...
    'filled','MarkerEdgeColor','k');
set(gcf,'Position',[1212 198 213 482])
ylabel('time shift ch2vsch1 (seconds)')
xticklabels([]);
colormap cool; colorbar

end

%%%%%%%%%%%%%%%%%%%%%%%%% Sub Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ht, realhv]=findThalf(trace,halfvalue,strupdown)

if strcmp(strupdown,'ascend')

    k=find(trace>halfvalue);
    if ~isempty(k) && k(1)>1
        FirstCross=k(1);
        AnteCross=k(1)-1;
        ht=interp1([trace(AnteCross) trace(FirstCross)],[AnteCross FirstCross],halfvalue);
        realhv=halfvalue;
    else
        ht=NaN;realhv=NaN;
    end

elseif strcmp(strupdown,'descend')

    k=find(trace>halfvalue);
    if ~isempty(k) && k(end)<numel(trace)
        AnteCross=k(end);
        LastCross=k(end)+1;
        if trace(AnteCross)>trace(LastCross)
            ht=interp1([trace(AnteCross) trace(LastCross)],[AnteCross LastCross],halfvalue);
            realhv=halfvalue;
        end
    else
        ht=NaN;realhv=NaN;
    end
end

end


function [halftime, shape, width]=findThalf2(trace,halfvalue)


halftime=NaN;
shape='undefined';
width=NaN;
ntrace=trace-halfvalue;
sntrace=sign(ntrace);
dsntrace=diff(sntrace);

[~,mi]=max(ntrace);


[id, ~, direc]=find(dsntrace);
if numel(id)>2 && numel(id)<=4
    dmi=diff(sign(mi-id));
    imm=find(dmi==-2);
    if ~isempty(imm)
        id=[id(imm);id(imm+1)];
        direc=[direc(imm);direc(imm+1)];
    elseif mi>id(end)
        id=id(end);
        direc=direc(end);
    elseif mi<id(1)
        id=id(1);
        direc=direc(1);
    end

end
rx=[];
for i=1:numel(id)
    rx=[rx;id(i);id(i)+1]; %#ok<AGROW> 
end
ry=ntrace(rx);
rx=reshape(rx,2,[]);
ry=reshape(ry,2,[]);
if ~isempty(id)
    halftime=zeros(numel(id),1);
end
for i=1:numel(id)
    halftime(i)=interp1(ry(:,i),rx(:,i),0);
end

ninter=numel(direc);

switch ninter
    case 1
        if direc>0
            shape='up';
            width=NaN;
        else
            shape='down';
            width=NaN;
        end
    case 2
        if direc(1)>0
            shape='bell';
            width=halftime(2)-halftime(1);
        else
            shape='ushape';
            width=NaN;
        end
end

shape=categorical(cellstr(shape));

end



