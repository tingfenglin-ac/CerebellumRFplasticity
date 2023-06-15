% 1. Add the path of the current script to the search path
% 2. Change the current folder to "PFstimulation_spatial pattern"
% 3. ctrl+enter to run each section

%% Section 1: Open files
% Select all file in the subfolders WT, L7SK2, and then CaMKII
% within the folder of "PFstimulation_spatial pattern"

% ctrl+enter to run this section


load('002_004_001_ball.mat')%10 sec
time=time(1:310)-time(154);

clearvars -except time
FN=[];
FP=[];
FNlist={};
FPlist={};
condlist=[];

f=1;
while f
    [FileName,FolderPath] = uigetfile('*.mat','Select files', 'Multiselect', 'on');

    if FolderPath==0;f=0;end

    if iscell(FileName)
        NewAddFile=size(FileName,2);
    elseif FileName~=0
        NewAddFile=1;
        if strfind(FileName,'SpikeAnalist')
            load([FolderPath,FileName])
            NewAddFile=size(FN,1);
        end
    else
        NewAddFile=0;
    end
    %%
    if NewAddFile~=0;
        for fnumber=1:NewAddFile
            if iscell(FileName)
                FNlist=cat(1,FNlist,FileName{fnumber});
                FPlist=cat(1,FPlist,FolderPath);
            elseif strfind(FileName,'SpikeAnalist')
                FNlist=cat(1,FNlist,FN{fnumber});
                FPlist=cat(1,FPlist,FP{fnumber});
            else
                FNlist=cat(1,FNlist,FileName);
                FPlist=cat(1,FPlist,FolderPath);
            end
            condlist=cat(1,condlist,f);
        end
        f=f+1;
    end
end

% import data
for i=1:size(FNlist,1)
    data(i)=load([FPlist{i},FNlist{i}]);
end

% analysize the calcium response across dendrites
co=lines;
co(6,:)=[.5 .5 .5];
cnum=[6 2 4 5];
for i=1:size(FNlist,1)
    % for i=66

    linethink=3;
    n=8;%filter bin
    data(i).yvalue=data(i).yscale.*([1:size(data(i).corrF_pre,1)]-1)';

    timewindow=[0.16 0.32];
    data(i).fPKS_pre=medfilt1(flipud(mean(data(i).corrF_pre(:,time>timewindow(1) & time<=timewindow(2)),2)),n);
    data(i).fPKS_early=medfilt1(flipud(mean(data(i).corrF_early(:,time>timewindow(1) & time<=timewindow(2)),2)),n);
    data(i).fPKS_late=medfilt1(flipud(mean(data(i).corrF_late(:,time>timewindow(1) & time<=timewindow(2)),2)),n);

    % figure
    % subplot(1,3,[2 3]);
    % hold on
    % plot(data(i).fPKS_pre,data(i).yvalue,'linewidth',linethink,'color',c(condlist(i)+1,:,1))
    % plot(data(i).fPKS_early,data(i).yvalue,'linewidth',linethink,'color',c(condlist(i)+1,:,2))
    % plot(data(i).fPKS_late,data(i).yvalue,'linewidth',linethink,'color',c(condlist(i)+1,:,3))
    % set(gca,'FontSize',20)
    % set(gcf,'color',[1 1 1])
    % yticks(0:50:200)
    % ylim([min(data(i).yvalue) max(data(i).yvalue)])
    % xlabel('\DeltaF/F')
    % ylabel('Dendritic length (\mum)')

    % imagewd = getframe(gcf);
    % imwrite(imagewd.cdata, [FNlist{i} '_dF.tiff']);
    % close all

end

% fill area of normalized response
for i=1:size(FNlist,1)
    pre=data(i).fPKS_pre/max(data(i).fPKS_pre);
    early=data(i).fPKS_early/max(data(i).fPKS_pre);
    late=data(i).fPKS_late/max(data(i).fPKS_pre);
    Fscale=100;
    imsize=[Fscale*3 size(data(i).yvalue,1)];

    mask_pre = poly2mask(Fscale*[0;pre;0], flipud([1 1:length(data(i).yvalue) length(data(i).yvalue)]'), imsize(2), imsize(1));
    mask_early = poly2mask(Fscale*[0;early;0], flipud([1 1:length(data(i).yvalue) length(data(i).yvalue)]'), imsize(2), imsize(1));
    mask_late = poly2mask(Fscale*[0;late;0], flipud([1 1:length(data(i).yvalue) length(data(i).yvalue)]'), imsize(2), imsize(1));
    % imagesc(data(i).mask_pre)
    % imagesc(data(i).mask_early)
    % imagesc(data(i).mask_late)

    data(i).NDlength_pre=sum(mask_pre,1);
    data(i).NDlength_early=sum(mask_early,1);
    data(i).NDlength_late=sum(mask_late,1);

    mask_pre = poly2mask([0:imsize(1)-1 0], [sum(mask_pre,1) 0], imsize(2), imsize(1));
    mask_early = poly2mask([0:imsize(1)-1 0], [sum(mask_early,1) 0], imsize(2), imsize(1));
    mask_late = poly2mask([0:imsize(1)-1 0], [sum(mask_late,1) 0], imsize(2), imsize(1));

    data(i).cumulativeF_pre=zeros(300,1);
    data(i).cumulativeF_early=zeros(300,1);
    data(i).cumulativeF_late=zeros(300,1);
    data(i).cumulativeF_pre(1:length(sum(mask_pre,2)))=sum(mask_pre,2);
    data(i).cumulativeF_early(1:length(sum(mask_pre,2)))=sum(mask_early,2);
    data(i).cumulativeF_late(1:length(sum(mask_pre,2)))=sum(mask_late,2);
end

%% Section 2: plot indiviual calcium response across dendrites
% ctrl+enter to run this section

[c]=grad_co;
for i=5
    linethink=2.5;
    figure
    hold on
    plot(data(i).fPKS_pre/max(data(i).fPKS_pre),data(i).yvalue,'linewidth',linethink,'color',c(condlist(i)+1,:,1))
    plot(data(i).fPKS_early/max(data(i).fPKS_pre),data(i).yvalue,'linewidth',linethink,'color',c(condlist(i)+1,:,2))
    % plot(data(i).fPKS_late/max(data(i).fPKS_pre),data(i).yvalue,'linewidth',linethink,'color',c(condlist(i)+1,:,2))
    set(gca,'FontSize',24)
    set(gcf,'color',[1 1 1])
    xticks(0:3)
    yticks(0:50:200)
    ylim([min(data(i).yvalue) max(data(i).yvalue)])
    xlim([0 2.5])
    ylabel('Dendritic location (\mum)')
end

%% Section 3: plot individual dendritic length
% ctrl+enter to run this section

% i=5 for WT; i=6 for L7SK2; i=139 for CaMKII TT305/6VA
for i=5
    figure;
    hold on
    x=((1:imsize(1))-1)/Fscale;
    y1=data(i).yscale.*data(i).NDlength_pre;y1(y1==0)=-1;
    y2=data(i).yscale.*data(i).NDlength_early;y2(y2==0)=-1;
    y3=data(i).yscale.*data(i).NDlength_late;y3(y3==0)=-1;
    idx=find(y1<y3);
    filly_pos=[y1 fliplr(y3)];filly_pos(idx)=y3(idx);
    fill([x fliplr(x)],[y1 fliplr(y2)],[0.97 0.91 0.81],'EdgeColor','none')
    % fill([x fliplr(x)],[y1 fliplr(y3)],[0.97 0.91 0.81],'EdgeColor','none')
    fill([x fliplr(x)],filly_pos,[182 208 226]/255,'EdgeColor','none')
    plot(x,y1,'linewidth',linethink,'color',c(condlist(i)+1,:,1))
    plot(x,y2,'linewidth',linethink,'color',c(condlist(i)+1,:,2))
    % plot(x,y3,'linewidth',linethink,'color',c(condlist(i)+1,:,2))
    ylim([min(data(i).yvalue) max(data(i).yvalue)])
    xlim([0 2.5])
    xticks(0:1:2)
    yticks(0:50:200)
    set(gca,'FontSize',24)
    set(gcf,'color',[1 1 1])
    ylabel('Cummulative length (\mum)')
end

%% Section 4: plot individual DEALTA dendritic length
% ctrl+enter to run this section

% i=5 for WT; i=6 for L7SK2; i=139 for CaMKII TT305/6VA
for i=5
    figure;
    % subplot(1,3,[2 3]);
    hold on
    x=((1:imsize)-1)/Fscale;
    y2=data(i).yscale.*(data(i).NDlength_early-data(i).NDlength_pre);
    y3=data(i).yscale.*(data(i).NDlength_late-data(i).NDlength_pre);
    idx=find(y2>0);
    filly_pos=[y2 0 0];filly_pos(idx)=0;
    idx=find(y3>0);
    % filly_pos=[y3 0 0];filly_pos(idx)=0;
    fill([x max(x) min(x)],[y3 0 0],[0.97 0.91 0.81],'EdgeColor','none')
    fill([x max(x) min(x)],filly_pos,[182 208 226]/255,'EdgeColor','none')

    % plot(x,y2,'linewidth',linethink,'color',c(condlist(i)+1,:,2))
    plot(x,y3,'linewidth',linethink,'color',c(condlist(i)+1,:,3))
    plot([max(x) min(x)],[0 0],'linewidth',linethink,'color',c(condlist(i)+1,:,1))

    xlim([0 2.5])
    set(gca,'FontSize',24)
    set(gcf,'color',[1 1 1])
    ylabel('\Deltadendritic length (\mum)')
    % ylim([-10 50])
    yticks(-200:40:100)
    xlim([0 2.5])
    xticks(0:3)
    ylim([-20 85])
    % imagewd = getframe(gcf);
    % imwrite(imagewd.cdata, [FNlist{i} '_DeltaNormDlength.tiff']);
    % close all

end
%% Section 5: Populatino average of Dendritic length across cells
% ctrl+enter to run this section

wt_NDlength_mat=[];
L7SK2_NDlength_mat=[];
CaMKII_NDlength_mat=[];

wt_NDlength_mat(:,:,1)=vertcat(data(condlist==1).NDlength_pre);
wt_NDlength_mat(:,:,2)=vertcat(data(condlist==1).NDlength_early);
wt_NDlength_mat(:,:,3)=vertcat(data(condlist==1).NDlength_late);

L7SK2_NDlength_mat(:,:,1)=vertcat(data(condlist==2).NDlength_pre);
L7SK2_NDlength_mat(:,:,2)=vertcat(data(condlist==2).NDlength_early);
L7SK2_NDlength_mat(:,:,3)=vertcat(data(condlist==2).NDlength_late);

CaMKII_NDlength_mat(:,:,1)=vertcat(data(condlist==3).NDlength_pre);
CaMKII_NDlength_mat(:,:,2)=vertcat(data(condlist==3).NDlength_early);
CaMKII_NDlength_mat(:,:,3)=vertcat(data(condlist==3).NDlength_late);

%Plot

% wt mean
linethink=2.5;
[co]=grad_co;

% WT
wt_mean=mean(wt_NDlength_mat,1);
wt_ser=nanstd(wt_NDlength_mat,0,1)./sqrt(size(wt_NDlength_mat,1));
x=((1:size(wt_NDlength_mat,2))-1)/Fscale;
% plot Normalized dendritic response area
figure
% subplot(1,3,[2 3]);
hold on
idx1_wt=find(wt_mean(:,:,1)~=0);
idx2_wt=find(wt_mean(:,:,2)~=0);
idx3_wt=find(wt_mean(:,:,3)~=0);
if length(idx1_wt)<length(x);idx1_wt=[idx1_wt idx1_wt(end)+1];end;% make the end of array 0
if length(idx2_wt)<length(x);idx2_wt=[idx2_wt idx2_wt(end)+1];end;
if length(idx3_wt)<length(x);idx3_wt=[idx3_wt idx3_wt(end)+1];end;


idx=find(wt_mean(:,:,2)>wt_mean(:,:,1));
filly_pos=[wt_mean(:,:,1) fliplr(wt_mean(:,:,2))];filly_pos(idx)=wt_mean(:,idx,2);
fill([x fliplr(x)],[wt_mean(:,:,1) fliplr(wt_mean(:,:,2))],[0.97 0.91 0.81],'EdgeColor','none')
% idx=find(wt_mean(:,:,3)>wt_mean(:,:,1));
% filly_pos=[wt_mean(:,:,1) fliplr(wt_mean(:,:,3))];filly_pos(idx)=wt_mean(:,idx,3);
% fill([x fliplr(x)],[wt_mean(:,:,1) fliplr(wt_mean(:,:,2))],[0.97 0.91 0.81],'EdgeColor','none')

fill([x fliplr(x)],filly_pos,[182 208 226]/255,'EdgeColor','none')
ErrArea_Smooth(x(idx1_wt)',wt_mean(:,idx1_wt,1)',wt_ser(:,idx1_wt,1)',...
    [co(2,:,1) 0.5]);
ErrArea_Smooth(x(idx2_wt)',wt_mean(:,idx2_wt,2)',wt_ser(:,idx2_wt,2)',...
    [co(2,:,2) .5]);
% ErrArea_Smooth(x(idx3_wt)',wt_mean(:,idx3_wt,3)',wt_ser(:,idx3_wt,3)',...
%     [co(2,:,2) .5]);
plot(x(idx1_wt),wt_mean(:,idx1_wt,1),'linewidth',linethink,'color',co(2,:,1))
plot(x(idx2_wt),wt_mean(:,idx2_wt,2),'linewidth',linethink,'color',co(2,:,2))
% plot(x(idx3_wt),wt_mean(:,idx3_wt,3),'linewidth',linethink,'color',co(2,:,2))
set(gca,'FontSize',24)
set(gcf,'color',[1 1 1])
ylabel('Sorted pixels (\mum)')
ylim([0 140])
xticks(0:3)
xlim([0 2.5])
% xmag=[0.7 0.8];ymag=[20 50];
% xlim(xmag);xticks(xmag);
% ylim(ymag);yticks(ymag);

% L7SK2 mean
L7SK2_mean=mean(L7SK2_NDlength_mat,1);
L7SK2_ser=nanstd(L7SK2_NDlength_mat,0,1)./sqrt(size(L7SK2_NDlength_mat,1));
% plot Normalized dendritic response area
figure
% subplot(1,3,[2 3]);
hold on
idx1_L7SK2=find(L7SK2_mean(:,:,1)~=0);
idx2_L7SK2=find(L7SK2_mean(:,:,2)~=0);
idx3_L7SK2=find(L7SK2_mean(:,:,3)~=0);
if length(idx1_L7SK2)<length(x);idx1_L7SK2=[idx1_L7SK2 idx1_L7SK2(end)+1];end;% make the end of array 0
if length(idx2_L7SK2)<length(x);idx2_L7SK2=[idx2_L7SK2 idx2_L7SK2(end)+1];end;
if length(idx3_L7SK2)<length(x);idx3_L7SK2=[idx3_L7SK2 idx3_L7SK2(end)+1];end;

idx=find(L7SK2_mean(:,:,2)>L7SK2_mean(:,:,1));
filly_pos=[L7SK2_mean(:,:,1) fliplr(L7SK2_mean(:,:,2))];filly_pos(idx)=L7SK2_mean(:,idx,2);
fill([x fliplr(x)],[L7SK2_mean(:,:,1) fliplr(L7SK2_mean(:,:,2))],[0.97 0.91 0.81],'EdgeColor','none')
% idx=find(L7SK2_mean(:,:,3)>L7SK2_mean(:,:,1));
% filly_pos=[L7SK2_mean(:,:,1) fliplr(L7SK2_mean(:,:,3))];filly_pos(idx)=L7SK2_mean(:,idx,3);
% fill([x fliplr(x)],[L7SK2_mean(:,:,1) fliplr(L7SK2_mean(:,:,2))],[0.97 0.91 0.81],'EdgeColor','none')

fill([x fliplr(x)],filly_pos,[182 208 226]/255,'EdgeColor','none')
ErrArea_Smooth(x(idx1_L7SK2)',L7SK2_mean(:,idx1_L7SK2,1)',L7SK2_ser(:,idx1_L7SK2,1)',...
    [co(3,:,1) 0.5]);
ErrArea_Smooth(x(idx2_L7SK2)',L7SK2_mean(:,idx2_L7SK2,2)',L7SK2_ser(:,idx2_L7SK2,2)',...
    [co(3,:,3) .5]);
% ErrArea_Smooth(x(idx3_L7SK2)',L7SK2_mean(:,idx3_L7SK2,3)',L7SK2_ser(:,idx3_L7SK2,3)',...
%     [co(3,:,3) .5]);
plot(x(idx1_L7SK2),L7SK2_mean(:,idx1_L7SK2,1),'linewidth',linethink,'color',co(3,:,1))
plot(x(idx2_L7SK2),L7SK2_mean(:,idx2_L7SK2,2),'linewidth',linethink,'color',co(3,:,3))
% plot(x(idx3_L7SK2),L7SK2_mean(:,idx3_L7SK2,3),'linewidth',linethink,'color',co(3,:,3))
set(gca,'FontSize',24)
set(gcf,'color',[1 1 1])
ylabel('Dendritic length (\mum)')
ylim([0 140])
xticks(0:3)
xlim([0 2.5])


% CaMKII mean
CaMKII_mean=mean(CaMKII_NDlength_mat,1);
CaMKII_ser=nanstd(CaMKII_NDlength_mat,0,1)./sqrt(size(CaMKII_NDlength_mat,1));
% plot Normalized dendritic response area
figure
% subplot(1,3,[2 3]);
hold on
idx1_CaMKII=find(CaMKII_mean(:,:,1)~=0);
idx2_CaMKII=find(CaMKII_mean(:,:,2)~=0);
idx3_CaMKII=find(CaMKII_mean(:,:,3)~=0);
if length(idx1_CaMKII)<length(x);idx1_CaMKII=[idx1_CaMKII idx1_CaMKII(end)+1];end;% make the end of array 0
if length(idx2_CaMKII)<length(x);idx2_CaMKII=[idx2_CaMKII idx2_CaMKII(end)+1];end;
if length(idx3_CaMKII)<length(x);idx3_CaMKII=[idx3_CaMKII idx3_CaMKII(end)+1];end;


idx=find(CaMKII_mean(:,:,2)>CaMKII_mean(:,:,1));
filly_pos=[CaMKII_mean(:,:,1) fliplr(CaMKII_mean(:,:,2))];filly_pos(idx)=CaMKII_mean(:,idx,2);
fill([x fliplr(x)],[CaMKII_mean(:,:,1) fliplr(CaMKII_mean(:,:,2))],[0.97 0.91 0.81],'EdgeColor','none')
% idx=find(CaMKII_mean(:,:,3)>CaMKII_mean(:,:,1));
% filly_pos=[CaMKII_mean(:,:,1) fliplr(CaMKII_mean(:,:,3))];filly_pos(idx)=CaMKII_mean(:,idx,3);
% fill([x fliplr(x)],[CaMKII_mean(:,:,1) fliplr(CaMKII_mean(:,:,2))],[0.97 0.91 0.81],'EdgeColor','none')

fill([x fliplr(x)],filly_pos,[182 208 226]/255,'EdgeColor','none')
ErrArea_Smooth(x(idx1_CaMKII)',CaMKII_mean(:,idx1_CaMKII,1)',CaMKII_ser(:,idx1_CaMKII,1)',...
    [co(4,:,1) 0.5]);
ErrArea_Smooth(x(idx2_CaMKII)',CaMKII_mean(:,idx2_CaMKII,2)',CaMKII_ser(:,idx2_CaMKII,2)',...
    [co(4,:,3) .5]);
% ErrArea_Smooth(x(idx3_CaMKII)',CaMKII_mean(:,idx3_CaMKII,3)',CaMKII_ser(:,idx3_CaMKII,3)',...
%     [co(4,:,3) .5]);

plot(x(idx1_CaMKII),CaMKII_mean(:,idx1_CaMKII,1),'linewidth',linethink,'color',co(4,:,1))
plot(x(idx2_CaMKII),CaMKII_mean(:,idx2_CaMKII,2),'linewidth',linethink,'color',co(4,:,3))
% plot(x(idx3_CaMKII),CaMKII_mean(:,idx3_CaMKII,3),'linewidth',linethink,'color',co(4,:,3))
set(gca,'FontSize',24)
set(gcf,'color',[1 1 1])
ylabel('Dendritic length (\mum)')
ylim([0 140])
xticks(0:3)
xlim([0 2.5])
xlim([0 2.5])

% Magnification

% wt mean
linethink=2.5;
[co]=grad_co;

% WT
wt_mean=mean(wt_NDlength_mat,1);
wt_ser=nanstd(wt_NDlength_mat,0,1)./sqrt(size(wt_NDlength_mat,1));
x=((1:size(wt_NDlength_mat,2))-1)/Fscale;
% plot Normalized dendritic response area
figure
% subplot(1,3,[2 3]);
hold on
idx1_wt=find(wt_mean(:,:,1)~=0);
idx2_wt=find(wt_mean(:,:,2)~=0);
idx3_wt=find(wt_mean(:,:,3)~=0);
if length(idx1_wt)<length(x);idx1_wt=[idx1_wt idx1_wt(end)+1];end;% make the end of array 0
if length(idx2_wt)<length(x);idx2_wt=[idx2_wt idx2_wt(end)+1];end;
if length(idx3_wt)<length(x);idx3_wt=[idx3_wt idx3_wt(end)+1];end;


% idx=find(wt_mean(:,:,2)>wt_mean(:,:,1));
% filly_pos=[wt_mean(:,:,1) fliplr(wt_mean(:,:,2))];filly_pos(idx)=wt_mean(:,idx,2);
% fill([x fliplr(x)],[wt_mean(:,:,1) fliplr(wt_mean(:,:,2))],[0.97 0.91 0.81],'EdgeColor','none')
idx=find(wt_mean(:,:,3)>wt_mean(:,:,1));
filly_pos=[wt_mean(:,:,1) fliplr(wt_mean(:,:,3))];filly_pos(idx)=wt_mean(:,idx,3);
fill([x fliplr(x)],[wt_mean(:,:,1) fliplr(wt_mean(:,:,2))],[0.97 0.91 0.81],'EdgeColor','none')

fill([x fliplr(x)],filly_pos,[182 208 226]/255,'EdgeColor','none')
ErrArea_Smooth(x(idx1_wt)',wt_mean(:,idx1_wt,1)',wt_ser(:,idx1_wt,1)',...
    [co(2,:,1) 0.5]);
ErrArea_Smooth(x(idx2_wt)',wt_mean(:,idx2_wt,2)',wt_ser(:,idx2_wt,2)',...
    [co(2,:,2) .5]);
% ErrArea_Smooth(x(idx3_wt)',wt_mean(:,idx3_wt,3)',wt_ser(:,idx3_wt,3)',...
%     [co(2,:,2) .5]);
plot(x(idx1_wt),wt_mean(:,idx1_wt,1),'linewidth',linethink,'color',co(2,:,1))
plot(x(idx2_wt),wt_mean(:,idx2_wt,2),'linewidth',linethink,'color',co(2,:,2))
% plot(x(idx3_wt),wt_mean(:,idx3_wt,3),'linewidth',linethink,'color',co(2,:,2))
set(gca,'FontSize',24)
set(gcf,'color',[1 1 1])
ylabel('Sorted pixels (\mum)')
ylim([0 140])
xticks(0:3)
xlim([0 2.5])
ylim([55 65])
xlim([0 2.5]/14+0.42)


% L7SK2 mean
L7SK2_mean=mean(L7SK2_NDlength_mat,1);
L7SK2_ser=nanstd(L7SK2_NDlength_mat,0,1)./sqrt(size(L7SK2_NDlength_mat,1));
% plot Normalized dendritic response area
figure
% subplot(1,3,[2 3]);
hold on
idx1_L7SK2=find(L7SK2_mean(:,:,1)~=0);
idx2_L7SK2=find(L7SK2_mean(:,:,2)~=0);
idx3_L7SK2=find(L7SK2_mean(:,:,3)~=0);
if length(idx1_L7SK2)<length(x);idx1_L7SK2=[idx1_L7SK2 idx1_L7SK2(end)+1];end;% make the end of array 0
if length(idx2_L7SK2)<length(x);idx2_L7SK2=[idx2_L7SK2 idx2_L7SK2(end)+1];end;
if length(idx3_L7SK2)<length(x);idx3_L7SK2=[idx3_L7SK2 idx3_L7SK2(end)+1];end;

% idx=find(L7SK2_mean(:,:,2)>L7SK2_mean(:,:,1));
% filly_pos=[L7SK2_mean(:,:,1) fliplr(L7SK2_mean(:,:,2))];filly_pos(idx)=L7SK2_mean(:,idx,2);
% fill([x fliplr(x)],[L7SK2_mean(:,:,1) fliplr(L7SK2_mean(:,:,2))],[0.97 0.91 0.81],'EdgeColor','none')
idx=find(L7SK2_mean(:,:,3)>L7SK2_mean(:,:,1));
filly_pos=[L7SK2_mean(:,:,1) fliplr(L7SK2_mean(:,:,3))];filly_pos(idx)=L7SK2_mean(:,idx,3);
fill([x fliplr(x)],[L7SK2_mean(:,:,1) fliplr(L7SK2_mean(:,:,2))],[0.97 0.91 0.81],'EdgeColor','none')

fill([x fliplr(x)],filly_pos,[182 208 226]/255,'EdgeColor','none')
ErrArea_Smooth(x(idx1_L7SK2)',L7SK2_mean(:,idx1_L7SK2,1)',L7SK2_ser(:,idx1_L7SK2,1)',...
    [co(3,:,1) 0.5]);
ErrArea_Smooth(x(idx2_L7SK2)',L7SK2_mean(:,idx2_L7SK2,2)',L7SK2_ser(:,idx2_L7SK2,2)',...
    [co(3,:,3) .5]);
% ErrArea_Smooth(x(idx3_L7SK2)',L7SK2_mean(:,idx3_L7SK2,3)',L7SK2_ser(:,idx3_L7SK2,3)',...
%     [co(3,:,3) .5]);
plot(x(idx1_L7SK2),L7SK2_mean(:,idx1_L7SK2,1),'linewidth',linethink,'color',co(3,:,1))
plot(x(idx2_L7SK2),L7SK2_mean(:,idx2_L7SK2,2),'linewidth',linethink,'color',co(3,:,3))
% plot(x(idx3_L7SK2),L7SK2_mean(:,idx3_L7SK2,3),'linewidth',linethink,'color',co(3,:,3))
set(gca,'FontSize',24)
set(gcf,'color',[1 1 1])
ylabel('Dendritic length (\mum)')
ylim([0 140])
xticks(0:3)
xlim([0 2.5])
xlim([0 2.5])
ylim([55 65])
xlim([0 2.5]/14+0.32)

% CaMKII mean
CaMKII_mean=mean(CaMKII_NDlength_mat,1);
CaMKII_ser=nanstd(CaMKII_NDlength_mat,0,1)./sqrt(size(CaMKII_NDlength_mat,1));
% plot Normalized dendritic response area
figure
% subplot(1,3,[2 3]);
hold on
idx1_CaMKII=find(CaMKII_mean(:,:,1)~=0);
idx2_CaMKII=find(CaMKII_mean(:,:,2)~=0);
idx3_CaMKII=find(CaMKII_mean(:,:,3)~=0);
if length(idx1_CaMKII)<length(x);idx1_CaMKII=[idx1_CaMKII idx1_CaMKII(end)+1];end;% make the end of array 0
if length(idx2_CaMKII)<length(x);idx2_CaMKII=[idx2_CaMKII idx2_CaMKII(end)+1];end;
if length(idx3_CaMKII)<length(x);idx3_CaMKII=[idx3_CaMKII idx3_CaMKII(end)+1];end;


idx=find(CaMKII_mean(:,:,2)>CaMKII_mean(:,:,1));
filly_pos=[CaMKII_mean(:,:,1) fliplr(CaMKII_mean(:,:,2))];filly_pos(idx)=CaMKII_mean(:,idx,2);
fill([x fliplr(x)],[CaMKII_mean(:,:,1) fliplr(CaMKII_mean(:,:,2))],[0.97 0.91 0.81],'EdgeColor','none')
% idx=find(CaMKII_mean(:,:,3)>CaMKII_mean(:,:,1));
% filly_pos=[CaMKII_mean(:,:,1) fliplr(CaMKII_mean(:,:,3))];filly_pos(idx)=CaMKII_mean(:,idx,3);
% fill([x fliplr(x)],[CaMKII_mean(:,:,1) fliplr(CaMKII_mean(:,:,2))],[0.97 0.91 0.81],'EdgeColor','none')

fill([x fliplr(x)],filly_pos,[182 208 226]/255,'EdgeColor','none')
ErrArea_Smooth(x(idx1_CaMKII)',CaMKII_mean(:,idx1_CaMKII,1)',CaMKII_ser(:,idx1_CaMKII,1)',...
    [co(4,:,1) 0.5]);
ErrArea_Smooth(x(idx2_CaMKII)',CaMKII_mean(:,idx2_CaMKII,2)',CaMKII_ser(:,idx2_CaMKII,2)',...
    [co(4,:,3) .5]);
% ErrArea_Smooth(x(idx3_CaMKII)',CaMKII_mean(:,idx3_CaMKII,2)',CaMKII_ser(:,idx3_CaMKII,2)',...
%     [co(4,:,3) .5]);

plot(x(idx1_CaMKII),CaMKII_mean(:,idx1_CaMKII,1),'linewidth',linethink,'color',co(4,:,1))
plot(x(idx2_CaMKII),CaMKII_mean(:,idx2_CaMKII,2),'linewidth',linethink,'color',co(4,:,3))
% plot(x(idx3_CaMKII),CaMKII_mean(:,idx3_CaMKII,3),'linewidth',linethink,'color',co(4,:,3))
set(gca,'FontSize',24)
set(gcf,'color',[1 1 1])
ylabel('Dendritic length (\mum)')
ylim([0 140])
xticks(0:3)
xlim([0 2.5])
xlim([0 2.5])
ylim([55 65])
xlim([0 2.5]/14+0.32)

%% Section 6: plot Delta dendritic length (off-hotspot vs hotspot)
% ctrl+enter to run this section

% Statistics
HS_edg=[0.30 0.67 0.9];% criteria for hot spot
non_HS=x>HS_edg(1) & x<=HS_edg(2);% index for non-hotspot
HS=x>HS_edg(3);% index for hotspot

%SumArea_idx=[time hotspot genotype]
%Early
%Non-hotspot
%WT
SumArea=mean(wt_NDlength_mat(:,non_HS,2)-wt_NDlength_mat(:,non_HS,1),2);% wt early F<=1
SumArea_idx=repmat([1 1 1],size(wt_NDlength_mat,1),1);
%L7SK2
SumArea=[SumArea;mean(L7SK2_NDlength_mat(:,non_HS,2)-L7SK2_NDlength_mat(:,non_HS,1),2)];% L7SK2 early F<=1
SumArea_idx=[SumArea_idx;repmat([1 1 2],size(L7SK2_NDlength_mat,1),1)];
%CaMKII
SumArea=[SumArea;mean(CaMKII_NDlength_mat(:,non_HS,2)-CaMKII_NDlength_mat(:,non_HS,1),2)];% CaMKII early F<=1
SumArea_idx=[SumArea_idx;repmat([1 1 3],size(CaMKII_NDlength_mat,1),1)];

%Early
%hotspot
%WT
SumArea=[SumArea;mean(wt_NDlength_mat(:,HS,2)-wt_NDlength_mat(:,HS,1),2)];% wt early F>1
SumArea_idx=[SumArea_idx;repmat([1 2 1],size(wt_NDlength_mat,1),1)];
%L7SK2
SumArea=[SumArea;mean(L7SK2_NDlength_mat(:,HS,2)-L7SK2_NDlength_mat(:,HS,1),2)];% L7SK2 early F>1
SumArea_idx=[SumArea_idx;repmat([1 2 2],size(L7SK2_NDlength_mat,1),1)];
%CaMKII
SumArea=[SumArea;mean(CaMKII_NDlength_mat(:,HS,2)-CaMKII_NDlength_mat(:,HS,1),2)];% CaMKII early F>1
SumArea_idx=[SumArea_idx;repmat([1 2 3],size(CaMKII_NDlength_mat,1),1)];

%Late
%Non-hotspot
%WT
SumArea=[SumArea;mean(wt_NDlength_mat(:,non_HS,3)-wt_NDlength_mat(:,non_HS,1),2)];% wt late F<=1
SumArea_idx=[SumArea_idx;repmat([2 1 1],size(wt_NDlength_mat,1),1)];
%L7SK2
SumArea=[SumArea;mean(L7SK2_NDlength_mat(:,non_HS,3)-L7SK2_NDlength_mat(:,non_HS,1),2)];% L7SK2 late F<=1
SumArea_idx=[SumArea_idx;repmat([2 1 2],size(L7SK2_NDlength_mat,1),1)];
%CaMKII
SumArea=[SumArea;mean(CaMKII_NDlength_mat(:,non_HS,3)-CaMKII_NDlength_mat(:,non_HS,1),2)];% CaMKII late F<=1
SumArea_idx=[SumArea_idx;repmat([2 1 3],size(CaMKII_NDlength_mat,1),1)];

%Late
%hotspot
%WT
SumArea=[SumArea;mean(wt_NDlength_mat(:,HS,3)-wt_NDlength_mat(:,HS,1),2)];% wt early F>1
SumArea_idx=[SumArea_idx;repmat([2 2 1],size(wt_NDlength_mat,1),1)];
%L7SK2
SumArea=[SumArea;mean(L7SK2_NDlength_mat(:,HS,3)-L7SK2_NDlength_mat(:,HS,1),2)];% L7SK2 early F>1
SumArea_idx=[SumArea_idx;repmat([2 2 2],size(L7SK2_NDlength_mat,1),1)];
%CaMKII
SumArea=[SumArea;mean(CaMKII_NDlength_mat(:,HS,3)-CaMKII_NDlength_mat(:,HS,1),2)];% CaMKII early F>1
SumArea_idx=[SumArea_idx;repmat([2 2 3],size(CaMKII_NDlength_mat,1),1)];

[co]=grad_co;
% plot
space=0.22;
MarkerSize=50;
meanSize=5;
jrange=0.4;

figure;hold on
SumArea_x=(SumArea_idx(SumArea_idx(:,3)==1,1)-1)*8+(SumArea_idx(SumArea_idx(:,3)==1,2)-1)*4+1;
scatter(SumArea_x-space-0.5*jrange+jrange*rand(length(SumArea_x),1),SumArea(SumArea_idx(:,3)==1),MarkerSize,...
    'MarkerFaceColor','flat',...
    'MarkerEdgeColor','flat',...
    'MarkerFaceAlpha',0.2,...
    'MarkerEdgeAlpha',0.2,...
    'CData',co(2,:,2));
SumArea_x=(SumArea_idx(SumArea_idx(:,3)==2,1)-1)*8+(SumArea_idx(SumArea_idx(:,3)==2,2)-1)*4+2;
scatter(SumArea_x-space-0.5*jrange+jrange*rand(length(SumArea_x),1),SumArea(SumArea_idx(:,3)==2),MarkerSize,...
    'MarkerFaceColor','flat',...
    'MarkerEdgeColor','flat',...
    'MarkerFaceAlpha',0.2,...
    'MarkerEdgeAlpha',0.2,...
    'CData',co(3,:,2));
SumArea_x=(SumArea_idx(SumArea_idx(:,3)==3,1)-1)*8+(SumArea_idx(SumArea_idx(:,3)==3,2)-1)*4+3;
scatter(SumArea_x-space-0.5*jrange+jrange*rand(length(SumArea_x),1),SumArea(SumArea_idx(:,3)==3),MarkerSize,...
    'MarkerFaceColor','flat',...
    'MarkerEdgeColor','flat',...
    'MarkerFaceAlpha',0.2,...
    'MarkerEdgeAlpha',0.2,...
    'CData',co(4,:,2));

% 1-way ANOVA stat between on and off hotspot
gcom=[1 1;1 2;2 1;2 2];
clearvars p stats c m
for i=1:4
    stat_idx=SumArea_idx(:,1)==gcom(i,1) & SumArea_idx(:,2)==gcom(i,2);
    [p(i).anova1,~,stats(i).anova1] = anova1(SumArea(stat_idx),SumArea_idx(stat_idx,3),"off");
    [c(i).anova1 m(i).anova1]=multcompare(stats(i).anova1,'Display','off');
    arrayfun(@(a) errorbar((i-1)*4+a+space,m(i).anova1(a,1),m(i).anova1(a,2),'o-',...
        'color',co(a+1,:,2),...
        'CapSize',20,...
        'MarkerSize',meanSize,...
        'LineWidth',3),1:3,'uni',0);

    % plot significance of anova1
    Ysig=50;Yincre=5;sigincre=1;

    sigidx=find(c(i).anova1(:,end)<0.05 & c(i).anova1(:,end)>=0.01)';
    if sigidx
        for sn=sigidx
            plot((i-1)*4+c(i).anova1(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',2)
            text((i-1)*4+mean(c(i).anova1(sn,1:2)),Ysig+sigincre,'*','FontSize',24,...
                'HorizontalAlignment','center',...
                'color',[0 0 0]);
            Ysig=Ysig+Yincre;
        end
    end
    sigidx=find(c(i).anova1(:,end)<0.01 & c(i).anova1(:,end)>=0.001)';
    if sigidx

        for sn=sigidx
            plot((i-1)*4+c(i).anova1(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',2)
            text((i-1)*4+mean(c(i).anova1(sn,1:2)),Ysig+sigincre,'**','FontSize',24,...
                'HorizontalAlignment','center',...
                'color',[0 0 0]);
            Ysig=Ysig+Yincre;
        end
    end
    sigidx=find(c(i).anova1(:,end)<0.001)';
    if sigidx
        for sn=sigidx
            plot((i-1)*4+c(i).anova1(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',2)
            text((i-1)*4+mean(c(i).anova1(sn,1:2)),Ysig+sigincre,'***','FontSize',24,...
                'HorizontalAlignment','center',...
                'color',[0 0 0]);
            Ysig=Ysig+Yincre;
        end
    end
end
ylabel('\Deltadendritic length (\mum)')
yticks([-40:40:200])
xticks([2 6 10 14])
xticklabels({'Off (early-pre)' 'Hotspot(early-pre)' 'Off(late-pre)' 'Hotspot(late-pre)'})
set(gca,'FontSize',24)
set(gcf,'color',[1 1 1])

% plot population Delta D length

% WT
Delta_wt_mean=mean(cat(3,wt_NDlength_mat(:,:,2)-wt_NDlength_mat(:,:,1),wt_NDlength_mat(:,:,3)-wt_NDlength_mat(:,:,1)),1);
Delta_wt_ser=nanstd(cat(3,wt_NDlength_mat(:,:,2)-wt_NDlength_mat(:,:,1),wt_NDlength_mat(:,:,3)-wt_NDlength_mat(:,:,1)),0,1)./sqrt(size(wt_NDlength_mat,1));

% plot Normalized dendritic response area
figure
hold on

% early-pre
idx=find(Delta_wt_mean(:,:,1)>0);
filly_pos=[Delta_wt_mean(:,:,1) 0 0];filly_pos(idx)=0;
fill([x max(x) min(x)],[Delta_wt_mean(:,:,1) 0 0],[0.97 0.91 0.81],'EdgeColor','none')
% late-pre
% idx=find(Delta_wt_mean(:,:,2)>0);
% filly_pos=[Delta_wt_mean(:,:,2) 0 0];filly_pos(idx)=0;
% fill([x max(x) min(x)],[Delta_wt_mean(:,:,2) 0 0],[0.97 0.91 0.81],'EdgeColor','none')

fill([x max(x) min(x)],filly_pos,[182 208 226]/255,'EdgeColor','none')
% early-pre
ErrArea_Smooth(x',Delta_wt_mean(:,:,1)',Delta_wt_ser(:,:,1)',...
    [co(2,:,2) .5]);
% late-pre
% ErrArea_Smooth(x',Delta_wt_mean(:,:,2)',Delta_wt_ser(:,:,2)',...
%     [co(2,:,2) .5]);

% early-pre
plot(x,Delta_wt_mean(:,:,1),'linewidth',linethink,'color',co(2,:,2))
% late-pre
% plot(x,Delta_wt_mean(:,:,2),'linewidth',linethink,'color',co(2,:,2))
plot([max(x) min(x)],[0 0],'linewidth',linethink,'color',co(2,:,1))
set(gca,'FontSize',24)
set(gcf,'color',[1 1 1])
ylabel('\Deltadendritic length (\mum)')
ylim([-5 40])
yticks(0:20:60)
xticks(0:3)
xlim([0 2.5])

% L7SK2
Delta_L7SK2_mean=mean(cat(3,L7SK2_NDlength_mat(:,:,2)-L7SK2_NDlength_mat(:,:,1),L7SK2_NDlength_mat(:,:,3)-L7SK2_NDlength_mat(:,:,1)),1);
Delta_L7SK2_ser=nanstd(cat(3,L7SK2_NDlength_mat(:,:,2)-L7SK2_NDlength_mat(:,:,1),L7SK2_NDlength_mat(:,:,3)-L7SK2_NDlength_mat(:,:,1)),0,1)./sqrt(size(L7SK2_NDlength_mat,1));

% plot Normalized dendritic response area
figure
hold on
% early-pre
idx=find(Delta_L7SK2_mean(:,:,1)>0);
filly_pos=[Delta_L7SK2_mean(:,:,1) 0 0];filly_pos(idx)=0;
fill([x max(x) min(x)],[Delta_L7SK2_mean(:,:,1) 0 0],[0.97 0.91 0.81],'EdgeColor','none')
% late-pre
% idx=find(Delta_L7SK2_mean(:,:,2)>0);
% filly_pos=[Delta_L7SK2_mean(:,:,2) 0 0];filly_pos(idx)=0;
% fill([x max(x) min(x)],[Delta_L7SK2_mean(:,:,2) 0 0],[0.97 0.91 0.81],'EdgeColor','none')

fill([x max(x) min(x)],filly_pos,[182 208 226]/255,'EdgeColor','none')
% early-pre
ErrArea_Smooth(x',Delta_L7SK2_mean(:,:,1)',Delta_L7SK2_ser(:,:,1)',...
    [co(3,:,3) .6]);
% late-pre
% ErrArea_Smooth(x',Delta_L7SK2_mean(:,:,2)',Delta_L7SK2_ser(:,:,2)',...
%     [co(3,:,3) .6]);

% early-pre
plot(x,Delta_L7SK2_mean(:,:,1),'linewidth',linethink,'color',co(3,:,3))
% late-pre
% plot(x,Delta_L7SK2_mean(:,:,2),'linewidth',linethink,'color',co(3,:,3))
plot([max(x) min(x)],[0 0],'linewidth',linethink,'color',co(3,:,1))
set(gca,'FontSize',24)
set(gcf,'color',[1 1 1])
ylabel('\Deltadendritic length (\mum)')
ylim([-5 40])
yticks(0:20:60)
xticks(0:3)
xlim([0 2.5])

% CaMKII
Delta_CaMKII_mean=mean(cat(3,CaMKII_NDlength_mat(:,:,2)-CaMKII_NDlength_mat(:,:,1),CaMKII_NDlength_mat(:,:,3)-CaMKII_NDlength_mat(:,:,1)),1);
Delta_CaMKII_ser=nanstd(cat(3,CaMKII_NDlength_mat(:,:,2)-CaMKII_NDlength_mat(:,:,1),CaMKII_NDlength_mat(:,:,3)-CaMKII_NDlength_mat(:,:,1)),0,1)./sqrt(size(CaMKII_NDlength_mat,1));

% plot Normalized dendritic response area
figure
hold on
% early-pre
idx=find(Delta_CaMKII_mean(:,:,1)>0);
filly_pos=[Delta_CaMKII_mean(:,:,1) 0 0];filly_pos(idx)=0;
fill([x max(x) min(x)],[Delta_CaMKII_mean(:,:,1) 0 0],[0.97 0.91 0.81],'EdgeColor','none')
% late-pre
% idx=find(Delta_CaMKII_mean(:,:,2)>0);
% filly_pos=[Delta_CaMKII_mean(:,:,2) 0 0];filly_pos(idx)=0;
% fill([x max(x) min(x)],[Delta_CaMKII_mean(:,:,2) 0 0],[0.97 0.91 0.81],'EdgeColor','none')

fill([x max(x) min(x)],filly_pos,[182 208 226]/255,'EdgeColor','none')
% early-pre
ErrArea_Smooth(x',Delta_CaMKII_mean(:,:,1)',Delta_CaMKII_ser(:,:,1)',...
    [co(4,:,3) .6]);
% late-pre
% ErrArea_Smooth(x',Delta_CaMKII_mean(:,:,2)',Delta_CaMKII_ser(:,:,2)',...
%     [co(4,:,3) .6]);
% early-pre
plot(x,Delta_CaMKII_mean(:,:,1),'linewidth',linethink,'color',co(4,:,3))
% late-pre
% plot(x,Delta_CaMKII_mean(:,:,2),'linewidth',linethink,'color',co(4,:,3))
plot([max(x) min(x)],[0 0],'linewidth',linethink,'color',co(4,:,1))
set(gca,'FontSize',24)
set(gcf,'color',[1 1 1])
ylabel('\Deltadendritic length (\mum)')
ylim([-5 40])
yticks(0:20:60)
xticks(0:3)
xlim([0 2.5])

% Early-pre vs late-pre

% WT
Delta_wt_mean=mean(cat(3,wt_NDlength_mat(:,:,2)-wt_NDlength_mat(:,:,1),wt_NDlength_mat(:,:,3)-wt_NDlength_mat(:,:,1)),1);
Delta_wt_ser=nanstd(cat(3,wt_NDlength_mat(:,:,2)-wt_NDlength_mat(:,:,1),wt_NDlength_mat(:,:,3)-wt_NDlength_mat(:,:,1)),0,1)./sqrt(size(wt_NDlength_mat,1));
% plot Normalized dendritic response area
figure
% subplot(1,3,[2 3]);
hold on

% idx=find(Delta_wt_mean(:,:,1)>0);
% filly_pos=[Delta_wt_mean(:,:,1) 0 0];filly_pos(idx)=0;
% fill([x max(x) min(x)],[Delta_wt_mean(:,:,1) 0 0],[0.97 0.91 0.81],'EdgeColor','none')
idx=find(Delta_wt_mean(:,:,2)>0);
filly_pos=[Delta_wt_mean(:,:,2) 0 0];filly_pos(idx)=0;
% fill([x max(x) min(x)],[Delta_wt_mean(:,:,2) 0 0],[0.97 0.91 0.81],'EdgeColor','none')
%
% fill([x max(x) min(x)],filly_pos,[182 208 226]/255,'EdgeColor','none')
ErrArea_Smooth(x',Delta_wt_mean(:,:,1)',Delta_wt_ser(:,:,1)',...
    [co(2,:,1) .5]);
ErrArea_Smooth(x',Delta_wt_mean(:,:,2)',Delta_wt_ser(:,:,2)',...
    [co(2,:,2) .5]);
% plot(x,wt_NDlength_mat(:,:,2)-wt_NDlength_mat(:,:,1),'linewidth',1,'color',[[236 169 140]/255 0.3])
% plot(x,wt_NDlength_mat(:,:,3)-wt_NDlength_mat(:,:,1),'linewidth',1,'color',[[247 221 209]/255 0.3])

plot(x,Delta_wt_mean(:,:,1),'linewidth',linethink,'color',co(2,:,1))
plot(x,Delta_wt_mean(:,:,2),'linewidth',linethink,'color',co(2,:,2))
% plot([max(x) min(x)],[0 0],'linewidth',linethink,'color',co(2,:,1))
set(gca,'FontSize',24)
set(gcf,'color',[1 1 1])
ylabel('\Deltadendritic length (\mum)')
ylim([-5 40])
yticks(0:20:60)
xticks(0:3)
xlim([0 2.5])

% L7SK2

Delta_L7SK2_mean=mean(cat(3,L7SK2_NDlength_mat(:,:,2)-L7SK2_NDlength_mat(:,:,1),L7SK2_NDlength_mat(:,:,3)-L7SK2_NDlength_mat(:,:,1)),1);
Delta_L7SK2_ser=nanstd(cat(3,L7SK2_NDlength_mat(:,:,2)-L7SK2_NDlength_mat(:,:,1),L7SK2_NDlength_mat(:,:,3)-L7SK2_NDlength_mat(:,:,1)),0,1)./sqrt(size(L7SK2_NDlength_mat,1));
% plot Normalized dendritic response area
figure
% subplot(1,3,[2 3]);
hold on
% idx=find(Delta_L7SK2_mean(:,:,1)>0);
% filly_pos=[Delta_L7SK2_mean(:,:,1) 0 0];filly_pos(idx)=0;
% fill([x max(x) min(x)],[Delta_L7SK2_mean(:,:,1) 0 0],[0.97 0.91 0.81],'EdgeColor','none')
idx=find(Delta_L7SK2_mean(:,:,2)>0);
% filly_pos=[Delta_L7SK2_mean(:,:,2) 0 0];filly_pos(idx)=0;
% fill([x max(x) min(x)],[Delta_L7SK2_mean(:,:,2) 0 0],[0.97 0.91 0.81],'EdgeColor','none')

fill([x max(x) min(x)],filly_pos,[182 208 226]/255,'EdgeColor','none')
% plot(x,L7SK2_NDlength_mat(:,:,2)-L7SK2_NDlength_mat(:,:,1),'linewidth',1,'color',[[190 151 195]/255 0.3])
% plot(x,L7SK2_NDlength_mat(:,:,3)-L7SK2_NDlength_mat(:,:,1),'linewidth',1,'color',[[229 213 232]/255 0.3])

ErrArea_Smooth(x',Delta_L7SK2_mean(:,:,1)',Delta_L7SK2_ser(:,:,1)',...
    [co(3,:,1) .6]);
ErrArea_Smooth(x',Delta_L7SK2_mean(:,:,2)',Delta_L7SK2_ser(:,:,2)',...
    [co(3,:,3) .6]);

plot(x,Delta_L7SK2_mean(:,:,1),'linewidth',linethink,'color',co(3,:,1))
plot(x,Delta_L7SK2_mean(:,:,2),'linewidth',linethink,'color',co(3,:,3))
% plot([max(x) min(x)],[0 0],'linewidth',linethink,'color',co(3,:,1))
set(gca,'FontSize',24)
set(gcf,'color',[1 1 1])
ylabel('\Deltadendritic length (\mum)')
ylim([-5 20])
yticks(0:10:60)
xticks(0:3)
xlim([0 2.5])

% CaMKII

Delta_CaMKII_mean=mean(cat(3,CaMKII_NDlength_mat(:,:,2)-CaMKII_NDlength_mat(:,:,1),CaMKII_NDlength_mat(:,:,3)-CaMKII_NDlength_mat(:,:,1)),1);
Delta_CaMKII_ser=nanstd(cat(3,CaMKII_NDlength_mat(:,:,2)-CaMKII_NDlength_mat(:,:,1),CaMKII_NDlength_mat(:,:,3)-CaMKII_NDlength_mat(:,:,1)),0,1)./sqrt(size(CaMKII_NDlength_mat,1));
% plot Normalized dendritic response area
figure
% subplot(1,3,[2 3]);
hold on
% idx=find(Delta_CaMKII_mean(:,:,1)>0);
% filly_pos=[Delta_CaMKII_mean(:,:,1) 0 0];filly_pos(idx)=0;
% fill([x max(x) min(x)],[Delta_CaMKII_mean(:,:,1) 0 0],[0.97 0.91 0.81],'EdgeColor','none')
idx=find(Delta_CaMKII_mean(:,:,2)>0);
filly_pos=[Delta_CaMKII_mean(:,:,2) 0 0];filly_pos(idx)=0;
% fill([x max(x) min(x)],[Delta_CaMKII_mean(:,:,2) 0 0],[0.97 0.91 0.81],'EdgeColor','none')
%
% fill([x max(x) min(x)],filly_pos,[182 208 226]/255,'EdgeColor','none')
% plot(x,CaMKII_NDlength_mat(:,:,2)-CaMKII_NDlength_mat(:,:,1),'linewidth',1,'color',[[160 197 110]/255 0.3])
% plot(x,CaMKII_NDlength_mat(:,:,3)-CaMKII_NDlength_mat(:,:,1),'linewidth',1,'color',[[187 217 151]/255 0.3])

ErrArea_Smooth(x',Delta_CaMKII_mean(:,:,1)',Delta_CaMKII_ser(:,:,1)',...
    [co(4,:,1) .6]);
ErrArea_Smooth(x',Delta_CaMKII_mean(:,:,2)',Delta_CaMKII_ser(:,:,2)',...
    [co(4,:,3) .6]);
plot(x,Delta_CaMKII_mean(:,:,1),'linewidth',linethink,'color',co(4,:,1))
plot(x,Delta_CaMKII_mean(:,:,2),'linewidth',linethink,'color',co(4,:,3))
% plot([max(x) min(x)],[0 0],'linewidth',linethink,'color',[0.5 .5 .5])
set(gca,'FontSize',24)
set(gcf,'color',[1 1 1])
ylabel('\Deltadendritic length (\mum)')
ylim([-5 20])
yticks(0:10:60)
xticks(0:3)
xlim([0 2.5])

%% Section7: stat for distribution shift of DELTA dendritic length
% ctrl+enter to run this section

% statistics
samplesize=1;
offset=min([Delta_wt_mean],[],"all");

%Early-pre
%WT
tempHSumArea=arrayfun(@(a) repmat(x(a),round(samplesize*(Delta_wt_mean(:,a,1)-offset)),1),1:length(x),'uni',0);% wt early
HSumArea=vertcat(tempHSumArea{:});
HSumArea_idx=repmat([1 1],length(vertcat(tempHSumArea{:})),1);
%L7SK2
tempHSumArea=arrayfun(@(a) repmat(x(a),round(samplesize*(Delta_L7SK2_mean(:,a,1)-offset)),1),1:length(x),'uni',0);% wt early
HSumArea=[HSumArea;vertcat(tempHSumArea{:})];
HSumArea_idx=[HSumArea_idx;repmat([1 2],length(vertcat(tempHSumArea{:})),1)];
%CaMKII
tempHSumArea=arrayfun(@(a) repmat(x(a),round(samplesize*(Delta_CaMKII_mean(:,a,1)-offset)),1),1:length(x),'uni',0);% wt early
HSumArea=[HSumArea;vertcat(tempHSumArea{:})];
HSumArea_idx=[HSumArea_idx;repmat([1 3],length(vertcat(tempHSumArea{:})),1)];

%Late-pre
%WT
tempHSumArea=arrayfun(@(a) repmat(x(a),round(samplesize*(Delta_wt_mean(:,a,2)-offset)),1),1:length(x),'uni',0);% wt early
HSumArea=[HSumArea;vertcat(tempHSumArea{:})];
HSumArea_idx=[HSumArea_idx;repmat([2 1],length(vertcat(tempHSumArea{:})),1)];
%L7SK2
tempHSumArea=arrayfun(@(a) repmat(x(a),round(samplesize*(Delta_L7SK2_mean(:,a,2)-offset)),1),1:length(x),'uni',0);% wt early
HSumArea=[HSumArea;vertcat(tempHSumArea{:})];
HSumArea_idx=[HSumArea_idx;repmat([2 2],length(vertcat(tempHSumArea{:})),1)];
%CaMKII
tempHSumArea=arrayfun(@(a) repmat(x(a),round(samplesize*(Delta_CaMKII_mean(:,a,2)-offset)),1),1:length(x),'uni',0);% wt early
HSumArea=[HSumArea;vertcat(tempHSumArea{:})];
HSumArea_idx=[HSumArea_idx;repmat([2 3],length(vertcat(tempHSumArea{:})),1)];

% 1-way ANOVA stat between on and off hotspot
clearvars p stats c m
for i=1:2
    stat_idx=HSumArea_idx(:,1)==i;
    [p(i).anova1,~,stats(i).anova1] = kruskalwallis(HSumArea(stat_idx),HSumArea_idx(stat_idx,2));
    [c(i).anova1 m(i).anova1]=multcompare(stats(i).anova1,'Display','off',"CType","dunn-sidak");
end

% Plot
figure
hold on

yhorzbar=32;
yincre=3;
arrayfun(@(a) errorbar(median(HSumArea(HSumArea_idx(:,1)==1 & HSumArea_idx(:,2)==a)),...
    yhorzbar+(a-1)*yincre,...
    mad(HSumArea(HSumArea_idx(:,1)==1 & HSumArea_idx(:,2)==a)),...
    'horizontal','o-',...
    'color',co(a+1,:,2),...
    'CapSize',10,...
    'MarkerSize',meanSize,...
    'LineWidth',3),1:3,'uni',0);

for i=1
    % plot significance of anova1
    Xsig=max(arrayfun(@(a) median(HSumArea(HSumArea_idx(:,1)==i & HSumArea_idx(:,2)==a))+...
        mad(HSumArea(HSumArea_idx(:,1)==i & HSumArea_idx(:,2)==a)),1:3))+0.1;
    xincre=0.15;sigincre=0.08;yadjust=[-1.4 1.4];
    sigidx=find(c(i).anova1(:,end)<0.05 & c(i).anova1(:,end)>=0.01)';
    if sigidx
        for sn=sigidx
            plot([Xsig Xsig],yhorzbar+(c(i).anova1(sn,1:2)-1)*yincre+yadjust,'color',[0 0 0],'linewidth',2)
            text(Xsig+sigincre,yhorzbar+mean(c(i).anova1(sn,1:2)-1)*yincre,'*','FontSize',24,...
                'VerticalAlignment','middle',...
                'HorizontalAlignment','center',...
                'color',[0 0 0],...
                'Rotation',90);
            Xsig=Xsig+xincre;
        end
    end
    sigidx=find(c(i).anova1(:,end)<0.01 & c(i).anova1(:,end)>=0.001)';
    if sigidx

        for sn=sigidx
            plot([Xsig Xsig],yhorzbar+(c(i).anova1(sn,1:2)-1)*yincre+yadjust,'color',[0 0 0],'linewidth',2)
            text(Xsig+sigincre,yhorzbar+mean(c(i).anova1(sn,1:2)-1)*yincre,'**','FontSize',24,...
                'VerticalAlignment','middle',...
                'HorizontalAlignment','center',...
                'color',[0 0 0],...
                'Rotation',90);
            Xsig=Xsig+xincre;
        end
    end
    sigidx=find(c(i).anova1(:,end)<0.001)';
    if sigidx
        for sn=sigidx
            plot([Xsig Xsig],yhorzbar+(c(i).anova1(sn,1:2)-1)*yincre+yadjust,'color',[0 0 0],'linewidth',2)
            text(Xsig+sigincre,yhorzbar+mean(c(i).anova1(sn,1:2)-1)*yincre,'***','FontSize',24,...
                'VerticalAlignment','middle',...
                'HorizontalAlignment','center',...
                'color',[0 0 0],...
                'Rotation',90);
            Xsig=Xsig+xincre;
        end
    end
end

plot([max(x) min(x)],[0 0],'--','linewidth',1,'color',[0 0 0])
ErrArea_Smooth(x',Delta_wt_mean(:,:,1)',Delta_wt_ser(:,:,1)',...
    [co(2,:,2) .5]);
ErrArea_Smooth(x',Delta_L7SK2_mean(:,:,1)',Delta_L7SK2_ser(:,:,1)',...
    [co(3,:,2) .5]);
ErrArea_Smooth(x',Delta_CaMKII_mean(:,:,1)',Delta_CaMKII_ser(:,:,1)',...
    [co(4,:,2) .5]);

plot(x,Delta_wt_mean(:,:,1),'linewidth',linethink,'color',co(2,:,2))
plot(x,Delta_L7SK2_mean(:,:,1),'linewidth',linethink,'color',co(3,:,2))
plot(x,Delta_CaMKII_mean(:,:,1),'linewidth',linethink,'color',co(4,:,2))

set(gca,'FontSize',24)
set(gcf,'color',[1 1 1])
xlabel('Normalized \DeltaF/F')
ylabel('\Deltadendritic length (\mum)')
ylim([-5 40])
yticks(0:20:60)
xticks(0:3)
xlim([0 2.5])

% plot All Delta Normalized dendritic response area
figure
hold on

yhorzbar=41;
yincre=3;
arrayfun(@(a) errorbar(median(HSumArea(HSumArea_idx(:,1)==2 & HSumArea_idx(:,2)==a)),...
    yhorzbar+(a-1)*yincre,...
    mad(HSumArea(HSumArea_idx(:,1)==2 & HSumArea_idx(:,2)==a)),...
    'horizontal','o-',...
    'color',co(a+1,:,2),...
    'CapSize',10,...
    'MarkerSize',meanSize,...
    'LineWidth',3),1:3,'uni',0);

for i=2
    % plot significance of anova1
    Xsig=max(arrayfun(@(a) median(HSumArea(HSumArea_idx(:,1)==i & HSumArea_idx(:,2)==a))+...
        mad(HSumArea(HSumArea_idx(:,1)==i & HSumArea_idx(:,2)==a)),1:3))+0.1;
    xincre=0.15;sigincre=0.08;yadjust=[-1.4 1.4];
    sigidx=find(c(i).anova1(:,end)<0.05 & c(i).anova1(:,end)>=0.01)';
    if sigidx
        for sn=sigidx
            plot([Xsig Xsig],yhorzbar+(c(i).anova1(sn,1:2)-1)*yincre+yadjust,'color',[0 0 0],'linewidth',2)
            text(Xsig+sigincre,yhorzbar+mean(c(i).anova1(sn,1:2)-1)*yincre,'*','FontSize',24,...
                'VerticalAlignment','middle',...
                'HorizontalAlignment','center',...
                'color',[0 0 0],...
                'Rotation',90);
            Xsig=Xsig+xincre;
        end
    end
    sigidx=find(c(i).anova1(:,end)<0.01 & c(i).anova1(:,end)>=0.001)';
    if sigidx

        for sn=sigidx
            plot([Xsig Xsig],yhorzbar+(c(i).anova1(sn,1:2)-1)*yincre+yadjust,'color',[0 0 0],'linewidth',2)
            text(Xsig+sigincre,yhorzbar+mean(c(i).anova1(sn,1:2)-1)*yincre,'**','FontSize',24,...
                'VerticalAlignment','middle',...
                'HorizontalAlignment','center',...
                'color',[0 0 0],...
                'Rotation',90);
            Xsig=Xsig+xincre;
        end
    end
    sigidx=find(c(i).anova1(:,end)<0.001)';
    if sigidx
        for sn=sigidx
            plot([Xsig Xsig],yhorzbar+(c(i).anova1(sn,1:2)-1)*yincre+yadjust,'color',[0 0 0],'linewidth',2)
            text(Xsig+sigincre,yhorzbar+mean(c(i).anova1(sn,1:2)-1)*yincre,'***','FontSize',24,...
                'VerticalAlignment','middle',...
                'HorizontalAlignment','center',...
                'color',[0 0 0],...
                'Rotation',90);
            Xsig=Xsig+xincre;
        end
    end
end

plot([max(x) min(x)],[0 0],'--','linewidth',1,'color',[0 0 0])

ErrArea_Smooth(x',Delta_wt_mean(:,:,2)',Delta_wt_ser(:,:,2)',...
    [co(2,:,2) .5]);
ErrArea_Smooth(x',Delta_L7SK2_mean(:,:,2)',Delta_L7SK2_ser(:,:,2)',...
    [co(3,:,2) .5]);
ErrArea_Smooth(x',Delta_CaMKII_mean(:,:,2)',Delta_CaMKII_ser(:,:,2)',...
    [co(4,:,2) .5]);

plot(x,Delta_wt_mean(:,:,2),'linewidth',linethink,'color',co(2,:,2))
plot(x,Delta_L7SK2_mean(:,:,2),'linewidth',linethink,'color',co(3,:,2))
plot(x,Delta_CaMKII_mean(:,:,2),'linewidth',linethink,'color',co(4,:,2))

set(gca,'FontSize',24)
set(gcf,'color',[1 1 1])
xlabel('Normalized \DeltaF/F')
ylabel('\Deltaresponse area (\mum)')
ylim([-5 50])
yticks(0:20:60)
xticks(0:3)
xlim([0 2.5])