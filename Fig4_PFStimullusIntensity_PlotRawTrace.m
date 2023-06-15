% 1. Add the path of the current script to the search path
% 2. Change the current folder to "PFstimulation_stimulus dependence"
% 3. analyze the data within the folder of "PFstimulation_stimulus dependence"
% 4. Using ctrl+enter within each function to run each section

%% Section 1: Open files 
% Select SpikeAnalist file from 10uA_SpikeAnalist.mat, 20uA_SpikeAnalist.mat
% 30uA_SpikeAnalist.mat to 10psi_SpikeAnalist.mat
% ctrl+enter to run this section

clear
FN=[];
FP=[];
FNlist={};

FPlist={};
condlist=[];

f=1;
while f
    [FileName,FolderPath] = uigetfile({'*SpikeAna.mat;*SpikeAnalist.mat'},'Select SpikeAnalist', 'Multiselect', 'on');
    if FolderPath==0;f=0;end
    if iscell(FileName)
        NewAddFile=size(FileName,2);
    elseif FileName~=0
        NewAddFile=1;
        if strfind(FileName,'SpikeAnalist')
            load([FolderPath,FileName])
            NewAddFile=size(FileName,2);
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
clear data
for i=1:size(FNlist,1)
    data(i)=load([FPlist{i},FNlist{i}]);
end

% pull the data out
pksLOCS=arrayfun(@(x) data(x).pksLOCS,1:length(data),'uni',0);
time=arrayfun(@(x) data(x).time,1:length(data),'uni',0);
smoothBC_signal=arrayfun(@(x) data(x).smoothBC_signal,1:length(data),'uni',0);
ROIlabel=arrayfun(@(a) ['Cell ' num2str(a)],1:size(smoothBC_signal{1},2),'uni',0);

%% Section 2: plot average calcium traces of each stimulus conditions
% Randomly shuffling a time window within the range of SPana_window
% (default=[-4.5 -0.5]) to generate calcium trace without PF stimulation.
% A shuffle_idx.mat is saved to reproduce the randomly
% generated calcium traces shown in the paper. To reproduce exactly the
% same figure, please comment line 100, load shuffle_idx in the folder 
% of Airpuff_stimulus dependence, and then run this section.

function ave_overlap_trace(Marker)
%% ctrl+enter to run this section
increvalue=0;
co=lines;
co(6,:)=[.5 .5 .5];
cnum=[6 6 6 6];% when ploting WT tetanus group
delay=1.6; %delay of each spike
display_window=[-0.2 1.3]; % display of the time window for each spike
T=time{1}-time{1}(data(1).stim(1));
window_idx=T>display_window(1) & T<display_window(2);
linestyle={'-' '-' '-' '-'};
SPana_window=[-4.5 -0.5];
SPana_idx=find(T>SPana_window(1) & T<SPana_window(2));

y=[];
counts=cell(size(condlist));
for cond=1:max(condlist)
    condidx=find(condlist==cond)';
    ind_activ_cell=arrayfun(@(a) nanmean(smoothBC_signal{a},3),condidx,'UniformOutput',0);
    ind_activ{cond}=horzcat(ind_activ_cell{:});
end
T=time{1}-time{1}(data(1).stim(1));

figure; hold on
for cond=1:max(condlist)
    ave_activ{cond}=nanmean(ind_activ{cond},2);
    std_activ{cond}=nanstd(ind_activ{cond},0,2)./sqrt(size(ind_activ{cond},2));

    suffle_idx{cond}=arrayfun(@(a) SPana_idx(randperm(length(SPana_idx)-sum(window_idx)+1,1)),1:size(ind_activ{cond},2));
    suffle_activ{cond}=arrayfun(@(a,b) ind_activ{cond}(a:a+sum(window_idx)-1,b),suffle_idx{cond},1:length(suffle_idx{cond}),'uni',0)
    suffle_activ{cond}=horzcat(suffle_activ{cond}{:});
    suffle_ave_activ{cond}=nanmean(suffle_activ{cond},2);
    suffle_std_activ{cond}=nanstd(suffle_activ{cond},0,2)./sqrt(size(suffle_activ{cond},2));
end

% draw stim based on data
stimT=[T(data(1).stim)';T(data(1).stim)'];
stimT=stimT(:);
Stim_signal=[0 5 5 0]';
anaT=[0.05 0.4;0.05 0.4];
anaT=anaT(:);
ana_signal=[-1 2 2 -1]';

%plotting stimulus
b=area(stimT,Stim_signal,...
    'FaceColor',[255 96 0]/255,...
    'EdgeColor','none',...
    'FaceAlpha',0.4);
%plotting analysis window
b=area(anaT,Stim_signal,...
    'FaceColor',[249, 191, 69]./255,...
    'EdgeColor','none',...
    'FaceAlpha',0.4);
%plotting error
incre=increvalue.*(1:size(ind_activ{cond},2));
ErrArea_Smooth(T(window_idx),suffle_ave_activ{cond},...
    suffle_std_activ{cond},...
    [co(cnum(cond),:) .5]);
%plotting data
plot(T(window_idx),suffle_ave_activ{cond},linestyle{cond},...
    'color',[co(cnum(cond),:)],...
    'LineWidth',1.5);

for cond=1:max(condlist)
    %plotting stimulus
    b=area(stimT+(cond-1).*delay+delay,Stim_signal,-10,...
        'FaceColor',[255 96 0]/255,...
        'EdgeColor','none',...
        'FaceAlpha',0.4);

    %plotting analysis window
    b=area(anaT+(cond-1).*delay+delay,Stim_signal,...
        'FaceColor',[249, 191, 69]./255,...
        'EdgeColor','none',...
        'FaceAlpha',0.4);

    %plotting error
    incre=increvalue.*(1:size(ind_activ{cond},2));
    ErrArea_Smooth(T(window_idx)+(cond-1).*delay+delay,ave_activ{cond}(window_idx),...
        std_activ{cond}(window_idx),...
        [co(cnum(cond),:) .5]);
end

for cond=1:max(condlist)
    %     subplot(1,5,cond)
    incre=increvalue.*(1:size(ind_activ{cond},2));
    plot(T(window_idx)+(cond-1).*delay+delay,ave_activ{cond}(window_idx),linestyle{cond},...
        'color',[co(cnum(cond),:)],...
        'LineWidth',1.5);

    box off
    b.BaseLine.Visible = 'off'
    xlim([-0.3 7.9])
    ylim([0.0 0.15])
    ylim([0.0 2.5]-0.1)
    xlabel('Time (sec)')
    % ylabel('Amplitude (dF/F0)')
    xticks(-10:10)
    set(gca,'FontSize',18,'TickDir','out')
    set(gcf,'color',[1 1 1])
end
end

%% Section 3: plot average cells of responsive trace and spontaneous spike

function ave_responsive(Marker)
%%
LOCS=pksLOCS;
co=lines;
co(6,:)=[.5 .5 .5];
cnum=[6 6 6 6];% when ploting WT tetanus group
delay=1.6; %delay of each spike
display_window=[-0.2 1.3]; % display of the time window for each spike
T=time{1}-time{1}(data(1).stim(1));
window_idx=T>display_window(1) & T<display_window(2);
linestyle={'-' '-' '-' '-'};
ana_window=[0.05 0.4];
[~,shiftreference]=min(abs(T-mean(ana_window)));
SPana_window=[-4.5 -0.5];
SPana_idx=find(T>SPana_window(1) & T<SPana_window(2));

% spontaneous spike
SP_displaywindow=display_window-0.1;
spike_indx=[min(find(T>SP_displaywindow(1))):max(find(T<SP_displaywindow(2)))]-find(T==0); % index range relative to spike onset or spike peak

ind_rise_time=[];
spon_spike_idx=[];
spon_spike=[];
for i=1:length(LOCS)
    ind_rise_time{i}=cellfun(@(a) T(a),LOCS{i},'uni',0);
    spon_spike_idx{i}=cellfun(@(a,b) b(a+SP_displaywindow(1)>SPana_window(1) & a+SP_displaywindow(2)<SPana_window(2)),...
        ind_rise_time{i},LOCS{i},'uni',0);%index of spontaneous spike onset or peak

    for cellnum=1:size(spon_spike_idx{i},2)
        for sweep=1:size(spon_spike_idx{i},1)
            each_sweep_signal=smoothBC_signal{i}(:,cellnum,sweep)';% signal for each sweep
            if length(spon_spike_idx{i}{sweep,cellnum})~=0
                spon_spike{i}{sweep,cellnum}=nanmean(each_sweep_signal(spon_spike_idx{i}{sweep,cellnum}+spike_indx),1)';%average spike signal for each sweep
            end
        end
    end
end
%
ave_spon_spike=[];
for i=1:length(spon_spike_idx)
    for cellnum=1:size(spon_spike_idx{i},2)
        ave_spon_spike{i}{cellnum}=nanmean(horzcat(spon_spike{i}{:,cellnum}),2);
    end
end

ave_spon_spike=horzcat(ave_spon_spike{:});
ave_spon_spike=horzcat(ave_spon_spike{:});
SP_ave_activ=nanmean(ave_spon_spike,2);
SP_std_activ=nanstd(ave_spon_spike,0,2)./sqrt(size(ave_spon_spike,2));

%plotting spontaneous SP
figure; hold on
ErrArea_Smooth(T(window_idx),SP_ave_activ,...
    SP_std_activ,...
    [co(cnum(cond),:) .5]);
%plotting data
plot(T(window_idx),SP_ave_activ,linestyle{cond},...
    'color',[co(cnum(cond),:)],...
    'LineWidth',1.5);
% tactile response
y=[];
counts=cell(size(condlist));
for cond=1:max(condlist)
    condidx=find(condlist==cond)';
    ind_activ_cell=arrayfun(@(a) nanmean(smoothBC_signal{a},3),condidx,'UniformOutput',0);
    ind_std_cell=arrayfun(@(a) nanstd(smoothBC_signal{a},0,3),condidx,'UniformOutput',0);

    ind_activ{cond}=horzcat(ind_activ_cell{:});
    ind_std{cond}=horzcat(ind_std_cell{:});
end
T=time{1}-time{1}(data(1).stim(1));

for cond=1:max(condlist)
    condidx=find(condlist==cond)';
    mean_activ=[];
    for i=1:length(condidx)
        idx=condidx(i);
        for roi=1:size(LOCS{idx},2)
            % because the histcounts include left left edge, so
            % convert the sign, and then left-right flip.
            % because the histcounts include the right edge of the
            % last bin, so add one more bin at right
            counts{idx}{roi,1} = cellfun(@(a) logical(histcounts(fliplr(-T(a)),fliplr(-[ana_window(1) ana_window]))),LOCS{idx}(:,roi),'uni',0);
            counts{idx}{roi,1} = cellfun(@(a) a(1),counts{idx}{roi,1});

            % the shift to Ca signal to align the peak
            xshift{idx}{roi,1} = cellfun(@(a) min(a(T(a)>ana_window(1) & T(a)<=ana_window(2)))-shiftreference,LOCS{idx}(:,roi),'uni',0);
            %             xshift{idx}{roi,1} = cellfun(@(a) a(1),xshift{idx}{roi,1});
            activ=smoothBC_signal{idx}(:,roi,:);
            for sh=1:length(xshift{idx}{roi,1})
                if ~isempty(xshift{idx}{roi,1}{sh})
                    if xshift{idx}{roi,1}{sh}>=0
                        alignIDX=[xshift{idx}{roi,1}{sh}+1:length(activ),1:xshift{idx}{roi,1}{sh}]';
                    else
                        alignIDX=[length(activ)+1+xshift{idx}{roi,1}{sh}:length(activ),1:length(activ)+xshift{idx}{roi,1}{sh}]';
                    end
                    activ(:,:,sh)=activ(alignIDX,:,sh);
                end
            end


            activ=smoothBC_signal{idx}(:,roi,:);
            mean_activ{idx}(:,roi,:)=nanmean(activ(:,:,counts{idx}{roi,1}),3);
        end
    end
    ind_activ{cond}=horzcat(mean_activ{:});
end

for cond=1:max(condlist)
    ave_activ{cond}=nanmean(ind_activ{cond},2);
    std_activ{cond}=nanstd(ind_activ{cond},0,2)./sqrt(size(ind_activ{cond},2));
end

% draw stim based on data
stimT=[T(data(1).stim)';T(data(1).stim)'];
stimT=stimT(:);
Stim_signal=[0 5 5 0]';
anaT=[0.05 0.4;0.05 0.4];
anaT=anaT(:);
ana_signal=[0 5 5 0]';

for cond=1:max(condlist)
    %plotting stimulus
    b=area(stimT+(cond-1).*delay+delay,Stim_signal,-10,...
        'FaceColor',[255 96 0]/255,...
        'EdgeColor','none',...
        'FaceAlpha',0.4);

    %plotting analysis window
    b=area(anaT+(cond-1).*delay+delay,Stim_signal,...
        'FaceColor',[249, 191, 69]./255,...
        'EdgeColor','none',...
        'FaceAlpha',0.4);

    %plotting error
    ErrArea_Smooth(T(window_idx)+(cond-1).*delay+delay,ave_activ{cond}(window_idx),...
        std_activ{cond}(window_idx),...
        [co(cnum(cond),:) .5]);
end

for cond=1:max(condlist)
    hold on
    plot(T(window_idx)+(cond-1).*delay+delay,ave_activ{cond}(window_idx),linestyle{cond},...
        'color',[co(cnum(cond),:)],...
        'LineWidth',1.5);
    box off
    b.BaseLine.Visible = 'off'
    xlim([-0.3 7.9])
    ylim([0.0 2.5]-0.1)
    xlabel('Time (sec)')
    % ylabel('Amplitude (dF/F0)')
    xticks(-10:10)
    set(gca,'FontSize',18,'TickDir','out')
    set(gcf,'color',[1 1 1])
end
end

%% Section 4: plot average cells of NON-responsive trace and spontaneous spike

function ave_NONresponsive(Marker)
%%
LOCS=pksLOCS;
co=lines;
co(6,:)=[.5 .5 .5];
cnum=[6 6 6 6];% when ploting WT tetanus group
delay=1.6; %delay of each spike
display_window=[-0.2 1.3]; % display of the time window for each spike
T=time{1}-time{1}(data(1).stim(1));
window_idx=T>display_window(1) & T<display_window(2);
linestyle={'-' '-' '-' '-'};
ana_window=[0.05 0.4];
SPana_window=[-4.5 -0.5];
SPana_idx=find(T>SPana_window(1) & T<SPana_window(2));

% spontaneous spike
spike_indx=[min(find(T>display_window(1))):max(find(T<display_window(2)))]-find(T==0); % index range relative to spike onset or spike peak

ind_rise_time=[];
spon_spike_idx=[];
spon_spike=[];
for i=1:length(LOCS)
    ind_rise_time{i}=cellfun(@(a) T(a),LOCS{i},'uni',0);
    spon_spike_idx{i}=cellfun(@(a,b) b(a+display_window(1)>SPana_window(1) & a+display_window(2)<SPana_window(2)),...
        ind_rise_time{i},LOCS{i},'uni',0);%index of spontaneous spike onset or peak

    for cellnum=1:size(spon_spike_idx{i},2)
        for sweep=1:size(spon_spike_idx{i},1)
            each_sweep_signal=smoothBC_signal{i}(:,cellnum,sweep)';% signal for each sweep
            if length(spon_spike_idx{i}{sweep,cellnum})~=0
                spon_spike{i}{sweep,cellnum}=nanmean(each_sweep_signal(spon_spike_idx{i}{sweep,cellnum}+spike_indx),1)';%average spike signal for each sweep
            end
        end
    end
end
%
ave_spon_spike=[];
for i=1:length(spon_spike_idx)
    for cellnum=1:size(spon_spike_idx{i},2)
        ave_spon_spike{i}{cellnum}=nanmean(horzcat(spon_spike{i}{:,cellnum}),2);
    end
end

ave_spon_spike=horzcat(ave_spon_spike{:});
ave_spon_spike=horzcat(ave_spon_spike{:});
SP_ave_activ=nanmean(ave_spon_spike,2);
SP_std_activ=nanstd(ave_spon_spike,0,2)./sqrt(size(ave_spon_spike,2));

%plotting spontaneous SP
figure; hold on

% tactile response
y=[];
counts=cell(size(condlist));
for cond=1:max(condlist)
    condidx=find(condlist==cond)';
    ind_activ_cell=arrayfun(@(a) nanmean(smoothBC_signal{a},3),condidx,'UniformOutput',0);
    ind_std_cell=arrayfun(@(a) nanstd(smoothBC_signal{a},0,3),condidx,'UniformOutput',0);

    ind_activ{cond}=horzcat(ind_activ_cell{:});
    ind_std{cond}=horzcat(ind_std_cell{:});
end
T=time{1}-time{1}(data(1).stim(1));

for cond=1:max(condlist)
    condidx=find(condlist==cond)';
    mean_activ=[];
    for i=1:length(condidx)
        idx=condidx(i);
        for roi=1:size(LOCS{idx},2)
            % because the histcounts include left left edge, so
            % convert the sign, and then left-right flip.
            % because the histcounts include the right edge of the
            % last bin, so add one more bin at right
            counts{idx}{roi,1} = cellfun(@(a) logical(~histcounts(fliplr(-T(a)),fliplr(-[ana_window(1) ana_window]))),LOCS{idx}(:,roi),'uni',0);
            counts{idx}{roi,1} = cellfun(@(a) a(1),counts{idx}{roi,1});

            activ=smoothBC_signal{idx}(:,roi,:);
            mean_activ{idx}(:,roi,:)=nanmean(activ(:,:,counts{idx}{roi,1}),3);
        end
    end
    ind_activ{cond}=horzcat(mean_activ{:});
end

for cond=1:max(condlist)
    ave_activ{cond}=nanmean(ind_activ{cond},2);
    std_activ{cond}=nanstd(ind_activ{cond},0,2)./sqrt(size(ind_activ{cond},2));
end

% draw stim based on data
stimT=[T(data(1).stim)';T(data(1).stim)'];
stimT=stimT(:);
Stim_signal=[0 5 5 0]';
anaT=[0.05 0.4;0.05 0.4];
anaT=anaT(:);
ana_signal=[-5 5 5 -5]';

for cond=1:max(condlist)
    %plotting stimulus
    b=area(stimT+(cond-1).*delay,Stim_signal,-10,...
        'FaceColor',[255 96 0]/255,...
        'EdgeColor','none',...
        'FaceAlpha',0.4);

    %plotting analysis window
    b=area(anaT+(cond-1).*delay,Stim_signal,...
        'FaceColor',[249, 191, 69]./255,...
        'EdgeColor','none',...
        'FaceAlpha',0.4);

    %plotting error
    ErrArea_Smooth(T(window_idx)+(cond-1).*delay,ave_activ{cond}(window_idx),...
        std_activ{cond}(window_idx),...
        [co(cnum(cond),:) .5]);
end

for cond=1:max(condlist)
    hold on
    plot(T(window_idx)+(cond-1).*delay,ave_activ{cond}(window_idx),linestyle{cond},...
        'color',[co(cnum(cond),:)],...
        'LineWidth',1.5);
    box off
    b.BaseLine.Visible = 'off'
    xlim([-0.3 7.9])
    ylim([0.0 2.5]-0.1)
    xlabel('Time (sec)')
    % ylabel('Amplitude (dF/F0)')
    xticks(-10:10)
    set(gca,'FontSize',18,'TickDir','out')
    set(gcf,'color',[1 1 1])
end
end
