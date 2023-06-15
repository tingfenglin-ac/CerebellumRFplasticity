% 1. Place this code in the data folder of
% "Dendrite-Soma-AIS simultaneous recording"
% 2. ctrl+enter to run each section

%% Section 1: Open files 
% Select all file in the subfolders Dendrite, Soma, and then AIS
% within the folder of "Dendrite-Soma-AIS simultaneous recording"

% ctrl+enter to run this section

clear
FN=[];
FP=[];
FNlist={};
FPlist={};
condlist=[];

f=1;
while f
    [FileName,FolderPath] = uigetfile({'*.mat'},'Select SpikeData files', 'Multiselect', 'on');
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
time=arrayfun(@(x) data(x).time,1:length(data),'uni',0);
smoothBC_signal=arrayfun(@(x) data(x).smoothBC_signal,1:length(data),'uni',0);
riseLOCS=arrayfun(@(x) data(x).riseLOCS,1:length(data),'uni',0);
upLOCS=arrayfun(@(x) data(x).upLOCS,1:length(data),'uni',0);
ind_rise_LOCS=cell(size(riseLOCS));
for i=1:length(data)
    ind_rise_LOCS{i}=cellfun(@(x,y) x(~ismember(x,y)),riseLOCS{i},upLOCS{i},'uni',0);
end
ROIlabel=arrayfun(@(a) ['Cell ' num2str(a)],1:size(smoothBC_signal{1},2),'uni',0);

%% Section 2: Average calcium signals for PF response of all cells (mean and sem)

function ave_overlap_trace(Marker)
%% ctrl+enter to run this section

increvalue=0;
co=lines;
co=[.5 .5 .5];
co=repmat(co,6,1);
cnum=[4 4 4 4];% when ploting WT tetanus group
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
% anaT=anaT(:);
% ana_signal=[-1 2 2 -1]';



for cond=1:max(condlist)
    figure
    subplot(1,2,1)
    hold on
    %plotting stimulus
    b=area(stimT,Stim_signal,-10,...
        'FaceColor',[255 96 0]/255,...
        'EdgeColor','none',...
        'FaceAlpha',0.4);
    %plotting analysis window
    anaT=[0.05 0.05 0.4 0.4]';
    b=area(anaT,Stim_signal,...
        'FaceColor',[249, 191, 69]./255,...
        'EdgeColor','none',...
        'FaceAlpha',0.4);
    %     %plotting analysis window
    %     anaT=[0.2 0.2 1 1]';
    %     b=area(anaT,Stim_signal,...
    %         'FaceColor',[0 0.4470 0.7410],...
    %         'EdgeColor','none',...
    %         'FaceAlpha',0.3);


    %     subplot(1,5,cond)
    incre=increvalue.*(1:size(ind_activ{cond},2));
    ErrArea_Smooth(T,ave_activ{cond},...
        std_activ{cond},...
        [co(cnum(cond),:) .5]);
    plot(T,ave_activ{cond},linestyle{cond},...
        'color',[co(cnum(cond),:)],...
        'LineWidth',1.5);

    box off
    b.BaseLine.Visible = 'off'
    xlim([-0.5 1.5])
    ylim([0.0 0.15])
    ylim([0.0 1.5]-0.1)
    %     ylim([0.0 1.5]-0.1)
    xlabel('Time (sec)')
    % ylabel('Amplitude (dF/F0)')
    yticks(-10:0.1:10)
    xticks(-10:10)
    set(gca,'FontSize',18,'TickDir','out')
    set(gcf,'color',[1 1 1])
end
end

%% Section 3: Average calcium signals for PF response of indivdual cells

function ind_reorder_trace(Marker)
%% ctrl+enter to run this section

increvalue=0.1;
yscale=[0 4.5];
co=[.5 .5 .5];
co=repmat(co,6,1);
cnum=[1 2 4 5];
amps=cell(size(condlist));

%
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

% reorder
cond_idx=1; % arrange the order based on condition of cond_indx
mean_range=[0.05 0.4];
range_idx=T>mean_range(1) & T<mean_range(2);
Mean_activ=mean(ind_activ{cond_idx}(range_idx,:),1);
[B,reorder_idx] = sort(Mean_activ);


for cond=1:max(condlist)
    incre=increvalue.*(1:size(ind_activ{cond},2));

    figure
    hold on
    %plotting stimulus
    b=area(stimT,Stim_signal,...
        'FaceColor',[255 96 0]/255,...
        'EdgeColor','none',...
        'FaceAlpha',0.4);
    %plotting analysis window
    anaT=[0.05 0.05 0.4 0.4]';
    b=area(anaT,Stim_signal,...
        'FaceColor',[249, 191, 69]./255,...
        'EdgeColor','none',...
        'FaceAlpha',0.4);

    for roi=1:size(ind_activ{cond},2)
        plot(T,ind_activ{cond}(:,reorder_idx(roi))+incre(roi),...
            'color',co(cond,:),...
            'LineWidth',1.5);
        ErrArea_Smooth(T,ind_activ{cond}(:,reorder_idx(roi))+incre(roi),...
            ind_std{cond}(:,reorder_idx(roi)),...
            co(cond,:));
    end

    % draw stim based on data
    stimT=[T(data(1).stim)';T(data(1).stim)'];
    stimT=stimT(:);
    Stim_signal=[0 50 50 0]';
    anaT=[0.05 0.2;0.05 0.2];
    anaT=anaT(:);


    %%
    xlim([-1 3])
    ylim(yscale)
    xlabel('Time (sec)')
    % ylabel('Amplitude (dF/F0)')
    %             xticks(-10:10)
    set(gca,'FontSize',18)
    set(gcf,'color',[1 1 1])
end
end

%% Section 4: Average calcium signals for spontaneous spikes of indivdual cells

function ave_spon_trace_for_ind_cell_reorder(Marker)
%% ctrl+enter to run this section

co=[.5 .5 .5];
co=repmat(co,6,1);
T=time{1}-time{1}(data(1).stim(1));
increvalue=.05;
yshift=1*increvalue;
yscale=[0 2.5];

sponwindow=[-5 0]; % time window of spontaneous activity
display_window=[-0.5 2]; % display of the time window for each spike
spike_indx=[min(find(T>display_window(1)))-find(T==0)+1:max(find(T<display_window(2)))-find(T==0)-1]; % index range relative to spike onset or spike peak

recordnum=size(smoothBC_signal,2)./3;%number of recording per condition (condition means PF, AIS and soma)
smoothBC_signal_reshape=[];
for cond=1:max(condlist)
    smoothBC_signal_reshape=[smoothBC_signal_reshape;smoothBC_signal((cond-1).*recordnum+[1:recordnum;])];
end

ind_rise_time=[];
spon_spike_idx=[];
spon_spike=[];
for i=1:recordnum
    ind_rise_time{i}=cellfun(@(a) T(a),ind_rise_LOCS{i},'uni',0);
    spon_spike_idx{i}=cellfun(@(a,b) b(a+display_window(1)>sponwindow(1) & a+display_window(2)<sponwindow(2)),...
        ind_rise_time{i},ind_rise_LOCS{i},'uni',0);%index of spontaneous spike onset or peak

    for cond=1:max(condlist)
        for cellnum=1:size(spon_spike_idx{i},2)
            for sweep=1:size(spon_spike_idx{i},1)
                each_sweep_signal=smoothBC_signal_reshape{cond,i}(:,cellnum,sweep)';% signal for each sweep
                if length(spon_spike_idx{i}{sweep,cellnum})~=0
                    spon_spike{cond,i}{cellnum,sweep}=nanmean(each_sweep_signal(spon_spike_idx{i}{sweep,cellnum}+spike_indx),1)';%average spike signal for each sweep
                end
            end
        end
    end
end

ave_spon_spike=[];
ste_spon_spike=[];
for i=1:recordnum
    for cond=1:max(condlist)
        for cellnum=1:size(spon_spike_idx{i},2)
            ave_spon_spike{i,cond}{cellnum}=nanmean(horzcat(spon_spike{cond,i}{cellnum,:}),2);
            ste_spon_spike{i,cond}{cellnum}=nanstd(horzcat(spon_spike{cond,i}{cellnum,:}),0,2)./sqrt(size(horzcat(spon_spike{cond,i}{cellnum,:}),2));
        end
    end
end
ave_spon_spike=horzcat(ave_spon_spike{:});
ave_spon_spike=horzcat(ave_spon_spike{:});
ste_spon_spike=horzcat(ste_spon_spike{:});
ste_spon_spike=horzcat(ste_spon_spike{:});

spon_spike_perCell=[];
spon_spike_err_perCell=[];
for cond=1:max(condlist)
    spon_spike_perCell{cond}=ave_spon_spike(:,(cond-1)*size(ind_activ{cond},2)+[1:size(ind_activ{cond},2)]);
    spon_spike_err_perCell{cond}=ste_spon_spike(:,(cond-1)*size(ind_activ{cond},2)+[1:size(ind_activ{cond},2)]);
end

% reorder
cond_idx=1; % arrange the order based on condition of cond_indx
mean_range=[0 0.2];
sponSpike_T=T(spike_indx+data(1).stim(1));
range_idx=sponSpike_T>mean_range(1) & sponSpike_T<mean_range(2);
Mean_activ=mean(spon_spike_perCell{cond_idx}(range_idx,:),1);
[~,closestIndex] = min(abs(sponSpike_T-0));%microadjust the y value

[B,reorder_idx] = sort(Mean_activ);

for cond=1:max(condlist)
    incre=increvalue.*(1:size(ind_activ{cond},2));

    figure
    hold on

    %plotting onset
    %     b=area(stimT,Stim_signal,...
    %         'FaceColor',[255 96 0]/255,...
    %         'EdgeColor','none',...
    %         'FaceAlpha',0.4);
    anaT=[0.05 0.05 0.4 0.4]';
    b=area(anaT,Stim_signal,...
        'FaceColor',[249, 191, 69]./255,...
        'EdgeColor','none',...
        'FaceAlpha',0.4);
    %     %plotting analysis window
    %     anaT=[0.2 0.2 1 1]';
    %     b=area(anaT,Stim_signal,...
    %         'FaceColor',[0 0.4470 0.7410],...
    %         'EdgeColor','none',...
    %         'FaceAlpha',0.3);

    for roi=1:size(ind_activ{cond},2)
        plot(sponSpike_T,spon_spike_perCell{cond}(:,reorder_idx(roi))+incre(roi)-spon_spike_perCell{cond}(closestIndex,reorder_idx(roi))+yshift,...
            'color',co(cond,:),...
            'LineWidth',1.5);
        ErrArea_Smooth(sponSpike_T,spon_spike_perCell{cond}(:,reorder_idx(roi))+incre(roi)-spon_spike_perCell{cond}(closestIndex,reorder_idx(roi))+yshift,...
            spon_spike_err_perCell{cond}(:,reorder_idx(roi)),...
            co(cond,:));
    end


    %%
    xlim([-1 3])
    ylim(yscale)
    xlabel('Time (sec)')
    % ylabel('Amplitude (dF/F0)')
    %             xticks(-10:10)
    set(gca,'FontSize',18)
    set(gcf,'color',[1 1 1])
end
end

%% Section 5: Average calcium signals for PF response of all cells (mean and sem)

function ave_stake_trace_of_spontaneousSpike(Marker)
%% ctrl+enter to run this section

for cond=1:max(condlist)
    ave_activ{cond}=nanmean(spon_spike_perCell{cond},2);
    std_activ{cond}=nanstd(spon_spike_perCell{cond},0,2)./sqrt(size(spon_spike_perCell{cond},2));
end

for cond=1:max(condlist)
    incre=increvalue.*(1:size(ind_activ{cond},2));
    figure
    subplot(1,2,1)
    hold on
    anaT=[0.05 0.05 0.4 0.4]';
    b=area(anaT,Stim_signal,...
        'FaceColor',[249, 191, 69]./255,...
        'EdgeColor','none',...
        'FaceAlpha',0.4);
    plot(sponSpike_T,ave_activ{cond},...
        'color',co(cond,:),...
        'LineWidth',1.5);
    ErrArea_Smooth(sponSpike_T,ave_activ{cond},...
        std_activ{cond},...
        co(cond,:));

    %%
    xlim([-0.5 1.5])
    ylim([-0 .3])
    xlabel('Time (sec)')
    xticks(-10:10)
    set(gca,'FontSize',18)
    set(gcf,'color',[1 1 1])
end
end

%% Section 6: max correlation of PF response

function correlation_max(Marker)
%% ctrl+enter to run this section

figure
xp=0.3;
yp=0.39;
co=lines;
cnum=[1 2 4 5];
xscale=[0 3];
yscale=[0 0.4001];
ytickscale=-1:0.2:2;

% assign the x variable
cond_idx=1; % 1=dendrite; 2=soma; 3=AIS
mean_range=[0.05 0.4];%time window for analysis
range_idx=T>mean_range(1) & T<=mean_range(2);
cond_activ_1=maxk(ind_activ{cond_idx}(range_idx,:),1);%pick up the maximum within time window

% assign the y variable
cond_idx=2;% 1=dendrite; 2=soma; 3=AIS
mean_range=[0.05 0.4];;%time window for analysis
range_idx=T>mean_range(1) & T<=mean_range(2);
cond_activ_2=maxk(ind_activ{cond_idx}(range_idx,:),1);%pick up the maximum within time window

% assign the y variable
cond_idx=3;% 1=dendrite; 2=soma; 3=AIS
mean_range=[0.05 0.4];;%time window for analysis
range_idx=T>mean_range(1) & T<=mean_range(2);
cond_activ_3=maxk(ind_activ{cond_idx}(range_idx,:),1);%pick up the maximum within time window

%linear fit pre
var1=cond_activ_1;var2=cond_activ_3;
mdl = fitlm(var1,var2,'Intercept',true);
beta = mdl.Coefficients.Estimate;
Xnew = linspace(min(cond_activ_1)-0.1, max(cond_activ_1)+0.1, 1000)';
[ypred,yci] = predict(mdl, Xnew,'Alpha',0.05);

hold on
plot(Xnew, ypred, '-','LineWidth',1.6,'Color',[.3 .3 .3 .3]);

marksize=40;
[R,P] = corrcoef(var1,var2);
scatter(var1,var2,marksize,...
    'MarkerEdgeColor','flat',...
    'MarkerEdgeAlpha',0.5,...
    'MarkerFaceAlpha',0.5,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor',[.3 .3 .3],...
    'linewidth',1.5);
xlabel('Dendrite (\DeltaF/F)');
ylabel('AIS (\DeltaF/F)');
text(xp,yp,['R = ' num2str(R(2),'%.3f') '; {\it p} = ' num2str(P(2),'%.3f')],...
    'FontSize',20);

xlim(xscale)
ylim(yscale)
yticks(ytickscale)
xticks(-2:4)
axis square
set(gca,'Fontsize',20,'linewidth',1.5)
set(gcf,'color',[1 1 1])
    end

%% Section 7: max correlation of spontaneous spike

function Spon_correlation_max(Marker)
%% ctrl+enter to run this section

figure
xp=0.1;
yp=0.39;
co=lines;
cnum=[1 2 4 5];
xscale=[0 0.8000001];
yscale=[0 0.4001];
ytickscale=-1:0.2:2;

% assign the x variable
cond_idx=1; % 1=dendrite; 2=soma; 3=AIS
mean_range=[0.05 0.4];%time window for analysis
range_idx=sponSpike_T>mean_range(1) & sponSpike_T<=mean_range(2);
cond_activ_1=maxk(spon_spike_perCell{cond_idx}(range_idx,:),1);%pick up the maximum within time window

% assign the y variable
cond_idx=2;% 1=dendrite; 2=soma; 3=AIS
mean_range=[0.05 0.4];;%time window for analysis
range_idx=sponSpike_T>mean_range(1) & sponSpike_T<=mean_range(2);
cond_activ_2=maxk(spon_spike_perCell{cond_idx}(range_idx,:),1);%pick up the maximum within time window

% assign the y variable
cond_idx=3;% 1=dendrite; 2=soma; 3=AIS
mean_range=[0.05 0.4];;%time window for analysis
range_idx=sponSpike_T>mean_range(1) & sponSpike_T<=mean_range(2);
cond_activ_3=maxk(spon_spike_perCell{cond_idx}(range_idx,:),1);%pick up the maximum within time window

%linear fit pre
var1=cond_activ_1;var2=cond_activ_3;
mdl = fitlm(var1,var2,'Intercept',true);
beta = mdl.Coefficients.Estimate;
Xnew = linspace(min(cond_activ_1)-0.1, max(cond_activ_1)+0.1, 1000)';
[ypred,yci] = predict(mdl, Xnew,'Alpha',0.05);

hold on
plot(Xnew, ypred, '-','LineWidth',1.6,'Color',[.3 .3 .3 .3]);

marksize=40;
[R,P] = corrcoef(var1,var2);
scatter(var1,var2,marksize,...
    'MarkerEdgeColor','flat',...
    'MarkerEdgeAlpha',0.5,...
    'MarkerFaceAlpha',0.5,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor',[.3 .3 .3],...
    'linewidth',1.5);
xlabel('Dendrite (\DeltaF/F)');
ylabel('AIS (\DeltaF/F)');
text(xp,yp,['R = ' num2str(R(2),'%.3f') '; {\it p} = ' num2str(P(2),'%.3f')],...
    'FontSize',20);

xlim(xscale)
ylim(yscale)
yticks(ytickscale)
xticks(-2:0.4:4)
axis square
set(gca,'Fontsize',20,'linewidth',1.5)
set(gcf,'color',[1 1 1])
    end
