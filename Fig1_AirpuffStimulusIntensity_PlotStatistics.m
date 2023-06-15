% 1. Add the path of the current script to the search path
% 2. Change the current folder to "Airpuff_stimulus dependence"
% 3. analyze the data within the folder of "Airpuff_stimulus dependence"
% 4. Using ctrl+enter within each function to run each section

function PlotSpikeAmp_StimIntensity

%% Section 1: Open files
% Select SpikeAnalist file from the weakest intensity to the strongest
% intensity (i.e., 04psi_SpikeAnalist.mat, 06psi_SpikeAnalist.mat,
% 08psi_SpikeAnalist.mat, to 10psi_SpikeAnalist.mat)
% ctrl+enter to run this section

clear
mn=1;% mouse#

FN=[];
FP=[];
FNlist={};
FPlist={};
condlist=[];
f=1;
filefolder=pwd;
while f
    [FileName,FolderPath] = uigetfile({'*SpikeAna.mat;*SpikeAnalist.mat'},'Select SpikeData files', 'Multiselect', 'on',filefolder);
    filefolder=FolderPath;
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
clear data;
for i=1:size(FNlist,1)
    data(i)=load([FPlist{i},FNlist{i}]);
end

% pull the data out
pksLOCS=arrayfun(@(x) data(x).pksLOCS,1:length(data),'uni',0);
time=arrayfun(@(x) data(x).time,1:length(data),'uni',0);
smoothBC_signal=arrayfun(@(x) data(x).smoothBC_signal,1:length(data),'uni',0);

% categorize the accumulated signal (during bust) from individual spike
% signal
StimInt=[];
for i=1:length(data)
StimInt(i)=str2num(FNlist{i}(end-17:end-16));
end
ROIlabel=arrayfun(@(a) ['Cell ' num2str(a)],1:size(smoothBC_signal{1},2),'uni',0);

mouse(mn).smoothBC_signal=smoothBC_signal;
mouse(mn).pksLOCS=pksLOCS;
mouse(mn).FNlist=FNlist;
mouse(mn).StimInt=StimInt;
mouse(mn).condlist=condlist;
mouse(mn).time=time;
mn = mn+1;

% Timewindow
ana_window=[0 0.35];% timewindow for response spike detection
SPana_window=[-4.5 -0.5];%timewindow for spontaneous activity

%% Section 2: plot statistics for probability
    function prob_state(Marker)
        %% ctrl+enter to run this section
        % plotting the probability within "window" sec after stim

        prob_by_cond=cell(size(mouse));
        StimInten=cell(size(mouse));
        for mn=1:length(mouse)
            LOCS=mouse(mn).pksLOCS;
            condlist=mouse(mn).condlist;
            counts=cell(size(condlist));
            SPcounts=cell(size(condlist));
            prob=cell(size(condlist));
            SPprob=cell(size(condlist));
            prob_by_cond{mn}=cell(max(condlist),1);
            SPprob_by_cond{mn}=cell(max(condlist),1);
            time=mouse(mn).time;
            y=[];
            for cond=1:max(condlist)
                %% histogram based on time
                %             cellnum=sum(cellfun(@(a) size(a,2),riseLOCS(condlist==cond)));
                condidx=find(condlist==cond)';
                StimInten{mn}=horzcat(StimInten{mn},mouse(mn).StimInt(condidx(1)));
                %%
                for i=1:length(condidx)
                    % average clls
                    idx=condidx(i);
                    ISinterval=mean(diff(time{idx}));
                    st_idx=data(1).stim(1);
                    T=time{idx}-time{idx}(st_idx+1);
                    for roi=1:size(LOCS{idx},2)

                        % because the histcounts include left left edge, so
                        % convert the sign, and then left-right flip.
                        % because the histcounts include the right edge of the
                        % last bin, so add one more bin at right
                        counts{idx}{roi,1} = cellfun(@(a) logical(histcounts(fliplr(-T(a)),fliplr(-[ana_window(1) ana_window]))),LOCS{idx}(:,roi),'uni',0);
                        counts{idx}{roi,1} = cellfun(@(a) a(1),counts{idx}{roi,1},'uni',0);
                        prob{idx}(roi,1) = sum(vertcat(counts{idx}{roi}{:}),1)./length(LOCS{idx}(:,roi));

                        SPcounts{idx}{roi,1} = cellfun(@(a) (histcounts(fliplr(-T(a)),fliplr(-[SPana_window(1) SPana_window]))./diff(SPana_window)).*diff(ana_window),LOCS{idx}(:,roi),'uni',0);
                        SPcounts{idx}{roi,1} = cellfun(@(a) a(1),SPcounts{idx}{roi,1});
                    end
                    prob_by_cond{mn}{cond}=vertcat(prob_by_cond{mn}{cond},prob{idx});
                    SPprob_by_cond{mn}{cond}=vertcat(SPprob_by_cond{mn}{cond},SPcounts{idx});
                end
            end
            SPprob_by_cond{mn}=horzcat(SPprob_by_cond{mn}{:});
            SPprob_by_cond{mn}=arrayfun(@(a) vertcat(SPprob_by_cond{mn}{a,:}),(1:size(SPprob_by_cond{mn},1))','UniformOutput',0);
            SPprob_by_cond{mn}=arrayfun(@(a) nanmean(SPprob_by_cond{mn}{a}),(1:size(SPprob_by_cond{mn},1))');
        end

        % statistics
        % 1-way ANOVA stat
        ydata=[];
        g2=[];
        for mn=1:length(mouse)
            ydata=[ydata;SPprob_by_cond{mn}];
            g2=[g2;repmat(0,length(SPprob_by_cond{mn}),1)];%condition
            for cond=1:length(prob_by_cond{mn})
                ydata=[ydata;prob_by_cond{mn}{cond}];
                g2=[g2;repmat(StimInten{mn}(cond),length(prob_by_cond{mn}{cond}),1)];%condition
            end
        end
        [p.anova,tb,stats.anova] = anova1(ydata,g2,'display','off');
        c.anova=multcompare(stats.anova,'Dimension',[1,2],'Display','off');
        g2mat=repmat(1:length(prob_by_cond{mn})+1,2,1);
        groupindex=[repmat(1:2,1,length(prob_by_cond{mn})+1);g2mat(1:end)]';%index of combinations and the corresponding groups
        [R,P]=arrayfun(@(a) corrcoef(g2,ydata),1:2,'uni',0);

        % plot response-stimulus curve
        co=lines;
        co(6,:)=[.5 .5 .5];
        cnum=[6 4 5];
        sigY=1.0;
        L_L_dist=0.017; %distance between significant lines
        L_D_dist=0.005; %distance between significant line and star symbol

        stim=[0 StimInten{1}];
        figure;hold on;
        x=1:5;
        y=[];
        ysem=[];
        for mn=1:length(mouse)
            y=[y mean([SPprob_by_cond{mn} horzcat(prob_by_cond{mn}{:})])'];
            ysem=[ysem nanstd([SPprob_by_cond{mn} horzcat(prob_by_cond{mn}{:})],1)'./sqrt(length(SPprob_by_cond{mn}))];
            errorbar(x,y(:,mn),ysem(:,mn),'o-',...
                'color',co(cnum(mn),:),...
                'CapSize',20,...
                'MarkerSize',5,...
                'LineWidth',3);
        end

        % significance
        sigX=0.5;
        sigY=0.55;
        dsigY=0.035;
        formatSpec = '%.3f';
        FS=16;

        if p.anova(1)<0.001
            text(sigX,sigY+dsigY,'{\it p} < 0.001','FontSize',FS)
        else
            text(sigX,sigY+dsigY,['{\it p} = ' num2str(p.anova(1),formatSpec)],'FontSize',FS)
        end

        pos=get(gcf,'position');
        pos(3)=350;
        set(gcf,'Position',pos);
        xticks(x)
        xticklabels({'SP' '4' '6' '8' '10'})
        xlim([0 6])
        yticks(0:0.1:10)
        ylim([0 0.25]+0.36)
        ylabel('Probability')
        xlabel('Airpuff stimulation (psi)')
        set(gca,'FontSize',18)
        set(gcf,'color',[1 1 1])

        % plot significance of anova over stimulus conditions
        Ysig=1;Yincre=0.05;
        gcomb=[1 1;2 2];% combinations of groups to compare
        sig=[];
        for mn=1
            sig{mn}=cell(5);
            idx=(1:5)';%find index of a genotype
            for col=1:5%plot significance of each stimulus intencity
                for row=find(stim~=(stim(col)))
                    if sum(ismember(c.anova(:,[1 2]),[idx(col) idx(row)],'rows'))
                        pvalue=c.anova(ismember(c.anova(:,[1 2]),[idx(col) idx(row)],'rows'),end);
                        if pvalue>=0.05
                            sig{mn}{col,row}=num2str(pvalue,'%1.3f');
                        elseif pvalue<0.05 & pvalue>=0.01
                            sig{mn}{col,row}=[num2str(pvalue,'%1.3f') '*'];
                        elseif pvalue<0.01 & pvalue>=0.001
                            sig{mn}{col,row}=[num2str(pvalue,'%1.3f') '**'];
                        else
                            sig{mn}{col,row}='<0.001***';
                        end
                    end
                end
            end
        end
    end

%% Section 3: plot statistics for normalized probability
    function Normalized_prob_withSP_state(Marker)
        %% ctrl+enter to run this section
        % plotting the probability within "window" sec after stim
        prob_by_cond=cell(size(mouse));
        StimInten=cell(size(mouse));
        for mn=1:length(mouse)
            LOCS=mouse(mn).pksLOCS;
            condlist=mouse(mn).condlist;
            counts=cell(size(condlist));
            SPcounts=cell(size(condlist));
            prob=cell(size(condlist));
            SPprob=cell(size(condlist));
            prob_by_cond{mn}=cell(max(condlist),1);
            SPprob_by_cond{mn}=cell(max(condlist),1);
            Norm_prob_by_cond{mn}=cell(max(condlist),1);
            SPNorm_prob_by_cond{mn}=cell(max(condlist),1);

            time=mouse(mn).time;
            y=[];
            for cond=1:max(condlist)
                %% histogram based on time
                %             cellnum=sum(cellfun(@(a) size(a,2),riseLOCS(condlist==cond)));
                condidx=find(condlist==cond)';
                StimInten{mn}=horzcat(StimInten{mn},mouse(mn).StimInt(condidx(1)));
                %%
                for i=1:length(condidx)
                    % average clls
                    idx=condidx(i);
                    ISinterval=mean(diff(time{idx}));
                    st_idx=data(1).stim(1);
                    T=time{idx}-time{idx}(st_idx+1);
                    for roi=1:size(LOCS{idx},2)

                        % because the histcounts include left left edge, so
                        % convert the sign, and then left-right flip.
                        % because the histcounts include the right edge of the
                        % last bin, so add one more bin at right
                        counts{idx}{roi,1} = cellfun(@(a) logical(histcounts(fliplr(-T(a)),fliplr(-[ana_window(1) ana_window]))),LOCS{idx}(:,roi),'uni',0);
                        counts{idx}{roi,1} = cellfun(@(a) a(1),counts{idx}{roi,1},'uni',0);
                        prob{idx}(roi,1) = sum(vertcat(counts{idx}{roi}{:}),1)./length(LOCS{idx}(:,roi));

                        SPcounts{idx}{roi,1} = cellfun(@(a) (histcounts(fliplr(-T(a)),fliplr(-[SPana_window(1) SPana_window]))./diff(SPana_window)).*diff(ana_window),LOCS{idx}(:,roi),'uni',0);
                        SPcounts{idx}{roi,1} = cellfun(@(a) a(1),SPcounts{idx}{roi,1});
                    end
                    prob_by_cond{mn}{cond}=vertcat(prob_by_cond{mn}{cond},prob{idx});
                    SPprob_by_cond{mn}{cond}=vertcat(SPprob_by_cond{mn}{cond},SPcounts{idx});
                end
            end
            SPprob_by_cond{mn}=horzcat(SPprob_by_cond{mn}{:});
            SPprob_by_cond{mn}=arrayfun(@(a) vertcat(SPprob_by_cond{mn}{a,:}),(1:size(SPprob_by_cond{mn},1))','UniformOutput',0);
            SPprob_by_cond{mn}=arrayfun(@(a) nanmean(SPprob_by_cond{mn}{a}),(1:size(SPprob_by_cond{mn},1))');
            SPNorm_prob_by_cond{mn}=SPprob_by_cond{mn}./SPprob_by_cond{mn};
            Norm_prob_by_cond{mn}=arrayfun(@(a) prob_by_cond{mn}{a}./SPprob_by_cond{mn},1:length(prob_by_cond{mn}),'UniformOutput',0);
        end

        % statistics
        % 1-way ANOVA stat
        ydata=[];
        g2=[];
        for mn=1:length(mouse)
            ydata=[ydata;SPNorm_prob_by_cond{mn}];
            g2=[g2;repmat(0,length(SPNorm_prob_by_cond{mn}),1)];%condition
            for cond=1:length(prob_by_cond{mn})
                ydata=[ydata;Norm_prob_by_cond{mn}{cond}];
                g2=[g2;repmat(StimInten{mn}(cond),length(Norm_prob_by_cond{mn}{cond}),1)];%condition
            end
        end
        [p.anova,tb,stats.anova] = anova1(ydata,g2);
        c.anova=multcompare(stats.anova,'Dimension',[1,2],'Display','off');
        g2mat=repmat(1:length(prob_by_cond{mn})+1,2,1);
        groupindex=[repmat(1:2,1,length(prob_by_cond{mn})+1);g2mat(1:end)]';%index of combinations and the corresponding groups
        [R,P]=corrcoef(g2,ydata);

        co=lines;
        co(6,:)=[.5 .5 .5];
        cnum=[6 4 5];
        % plot response-stimulus curve
        sigY=1.0;
        L_L_dist=0.017; %distance between significant lines
        L_D_dist=0.005; %distance between significant line and star symbol

        stim=[0 StimInten{1}];
        figure;hold on;
        x=1:5;
        y=[];
        ysem=[];
        for mn=1:length(mouse)
            y=[y mean([SPNorm_prob_by_cond{mn} horzcat(Norm_prob_by_cond{mn}{:})])'];
            ysem=[ysem nanstd([SPNorm_prob_by_cond{mn} horzcat(Norm_prob_by_cond{mn}{:})],1)'./sqrt(length(SPNorm_prob_by_cond{mn}))];
            errorbar(x,y(:,mn),ysem(:,mn),'o-',...
                'color',co(cnum(mn),:),...
                'CapSize',20,...
                'MarkerSize',5,...
                'LineWidth',3);
        end

        %significance
        sigX=0.5;
        sigY=1.65;
        formatSpec = '%.3f';
        FS=16;
        if p.anova<0.001
            text(sigX,sigY,'{\it p} < 0.001','FontSize',FS)
        else
            text(sigX,sigY,['{\it p} = ' num2str(p.anova,formatSpec)],'FontSize',FS)
        end

        pos=get(gcf,'position');
        pos(3)=350;
        set(gcf,'Position',pos);
        xticks(x)
        xticklabels({'SP' '4' '6' '8' '10'})
        xlim([0 6])
        yticks(1:0.3:10)
        ylim([0 0.8]+0.9)
        ylabel('Normalized probability')
        xlabel('Airpuff stimulation (psi)')
        set(gca,'FontSize',18)
        set(gcf,'color',[1 1 1])
    end

%% Section 4: categorize stimulus intensity-dependent cells according to normalized probability

    function Categorized_Normalized_prob_withSP_state(Marker)
        %% ctrl+enter to run this section
        % plotting the probability within "window" sec after stim
        prob_by_cond=cell(size(mouse));
        StimInten=cell(size(mouse));
        for mn=1:length(mouse)
            LOCS=mouse(mn).pksLOCS;
            condlist=mouse(mn).condlist;
            counts=cell(size(condlist));
            SPcounts=cell(size(condlist));
            prob=cell(size(condlist));
            SPprob=cell(size(condlist));
            prob_by_cond{mn}=cell(max(condlist),1);
            SPprob_by_cond{mn}=cell(max(condlist),1);
            Norm_prob_by_cond{mn}=cell(max(condlist),1);
            SPNorm_prob_by_cond{mn}=cell(max(condlist),1);

            time=mouse(mn).time;
            y=[];
            for cond=1:max(condlist)
                %% histogram based on time
                %             cellnum=sum(cellfun(@(a) size(a,2),riseLOCS(condlist==cond)));
                condidx=find(condlist==cond)';
                StimInten{mn}=horzcat(StimInten{mn},mouse(mn).StimInt(condidx(1)));
                %%
                for i=1:length(condidx)
                    % average clls
                    idx=condidx(i);
                    ISinterval=mean(diff(time{idx}));
                    st_idx=data(1).stim(1);
                    T=time{idx}-time{idx}(st_idx+1);
                    for roi=1:size(LOCS{idx},2)

                        % because the histcounts include left left edge, so
                        % convert the sign, and then left-right flip.
                        % because the histcounts include the right edge of the
                        % last bin, so add one more bin at right
                        counts{idx}{roi,1} = cellfun(@(a) logical(histcounts(fliplr(-T(a)),fliplr(-[ana_window(1) ana_window]))),LOCS{idx}(:,roi),'uni',0);
                        counts{idx}{roi,1} = cellfun(@(a) a(1),counts{idx}{roi,1},'uni',0);
                        prob{idx}(roi,1) = sum(vertcat(counts{idx}{roi}{:}),1)./length(LOCS{idx}(:,roi));

                        SPcounts{idx}{roi,1} = cellfun(@(a) (histcounts(fliplr(-T(a)),fliplr(-[SPana_window(1) SPana_window]))./diff(SPana_window)).*diff(ana_window),LOCS{idx}(:,roi),'uni',0);
                        SPcounts{idx}{roi,1} = cellfun(@(a) a(1),SPcounts{idx}{roi,1});
                    end
                    prob_by_cond{mn}{cond}=vertcat(prob_by_cond{mn}{cond},prob{idx});
                    SPprob_by_cond{mn}{cond}=vertcat(SPprob_by_cond{mn}{cond},SPcounts{idx});
                end
            end
            SPprob_by_cond{mn}=horzcat(SPprob_by_cond{mn}{:});
            SPprob_by_cond{mn}=arrayfun(@(a) vertcat(SPprob_by_cond{mn}{a,:}),(1:size(SPprob_by_cond{mn},1))','UniformOutput',0);
            SPprob_by_cond{mn}=arrayfun(@(a) nanmean(SPprob_by_cond{mn}{a}),(1:size(SPprob_by_cond{mn},1))');
            SPNorm_prob_by_cond{mn}=SPprob_by_cond{mn}./SPprob_by_cond{mn};
            Norm_prob_by_cond{mn}=arrayfun(@(a) prob_by_cond{mn}{a}./SPprob_by_cond{mn},1:length(prob_by_cond{mn}),'UniformOutput',0);
            Norm_prob_by_cond{mn}=horzcat(Norm_prob_by_cond{mn}{:});
        end

        % plot response-stimulus curve
        sigY=1.0;
        L_L_dist=0.017; %distance between significant lines
        L_D_dist=0.005; %distance between significant line and star symbol

        x=1:5;
        for mn=1:length(mouse)
            [R,P] = arrayfun(@(a) corrcoef(1:5,[SPNorm_prob_by_cond{mn}(a) Norm_prob_by_cond{mn}(a,:)]),(1:length(SPNorm_prob_by_cond{mn}))','uni',0);
            RP{mn} = arrayfun(@(a) [R{a}(2) P{a}(2)],(1:length(SPNorm_prob_by_cond{mn}))','uni',0);
            RP{mn}=vertcat(RP{mn}{:});
            corrindex=RP{mn}(:,1)>0.5 & RP{mn}(:,2)<0.5;

            %correlated
            figure;hold on;
            plot(x',[SPNorm_prob_by_cond{mn}(corrindex) Norm_prob_by_cond{mn}(corrindex,:)]','-',...
                'color',[co(cnum(mn),:) 0.3],...
                'MarkerSize',5,...
                'LineWidth',1);
            text(0.5,4.7,[num2str(sum(corrindex)) '/' num2str(length(corrindex))],'FontSize',18)
            pos=get(gcf,'position');
            pos(3)=350;
            set(gcf,'Position',pos);
            xticks(1:5)
            xticklabels({'SP' '4' '6' '8' '10'})
            yticks(0:2:10)
            xlim([0 6])
            ylim([0 5])
            ylabel('Normalized probability')
            xlabel('Airpuff stimulation (psi)')
            set(gca,'FontSize',18)
            set(gcf,'color',[1 1 1])

            %non-correlated
            figure;hold on;
            plot(x',[SPNorm_prob_by_cond{mn}(~corrindex) Norm_prob_by_cond{mn}(~corrindex,:)]','-',...
                'color',[co(cnum(mn),:) 0.3],...
                'MarkerSize',5,...
                'LineWidth',1);
            text(0.5,4.7,[num2str(sum(~corrindex)) '/' num2str(length(corrindex))],'FontSize',18)

            pos=get(gcf,'position');
            pos(3)=350;
            set(gcf,'Position',pos);
            xticks(1:5)
            xticklabels({'SP' '4' '6' '8' '10'})
            yticks(0:2:10)
            xlim([0 6])
            ylim([0 5])
            ylabel('Normalized probability')
            xlabel('Airpuff stimulation (psi)')
            set(gca,'FontSize',18)
            set(gcf,'color',[1 1 1])
        end
    end

%% Section 5: plot statistics for amplitude

    function First_withSP_amp_state(Marker)
        %% ctrl+enter to run this section
        StimInten=cell(size(mouse));
        amp_by_cond=cell(size(mouse));
        for mn=1:length(mouse)
            smoothBC_signal=mouse(mn).smoothBC_signal;
            LOCS=mouse(mn).pksLOCS;
            condlist=mouse(mn).condlist;

            amps=cell(size(condlist));
            SPamps=cell(size(condlist));
            amps_idx=cell(size(condlist));
            SPamps_idx=cell(size(condlist));
            cellmean=cell(size(condlist));
            SPcellmean=cell(size(condlist));
            amp_by_cond{mn}=cell(max(condlist),1);
            SPamp_by_cond{mn}=cell(max(condlist),1);
            Norm_amp_by_cond{mn}=cell(max(condlist),1);
            SPNorm_amp_by_cond{mn}=cell(max(condlist),1);

            time=mouse(mn).time;

            for cond=1:max(condlist)
                %% histogram based on time
                condidx=find(condlist==cond)';
                StimInten{mn}=horzcat(StimInten{mn},mouse(mn).StimInt(condidx(1)));
                %%
                for i=1:length(condidx)
                    % average clls
                    idx=condidx(i);
                    ISinterval=mean(diff(time{idx}));
                    st_idx=data(1).stim(1);
                    T=time{idx}-time{idx}(st_idx+1);

                    for roi=1:size(LOCS{idx},2)
                        activ=smoothBC_signal{idx}(:,roi,:);
                        ana_idx=LOCS{idx}(:,roi); %index of onset or peak depend on analysis

                        amps{idx}{roi,1} = arrayfun(@(a) activ(ana_idx{a}(T(ana_idx{a})>ana_window(1) & T(ana_idx{a})<=ana_window(2)),:,a),...
                            1:length(ana_idx),'uni',0);%filter out the data out of the ana_window
                        amps_idx{idx}{roi,1} = arrayfun(@(a) find(amps{idx}{roi,1}{a}),...
                            1:length(ana_idx),'uni',0);%find all the index of amplitude
                        amps_idx{idx}{roi,1} = arrayfun(@(a) min(amps_idx{idx}{roi,1}{a}),...
                            1:length(ana_idx),'uni',0);%find the first index
                        amps{idx}{roi,1} = arrayfun(@(a) amps{idx}{roi,1}{a}(amps_idx{idx}{roi,1}{a}),...
                            1:length(ana_idx),'uni',0);%first amplitude within the ana_window
                        cellmean{idx}(roi,1) = nanmean(vertcat(amps{idx}{roi}{:}),1);

                        SPamps{idx}{roi,1} = arrayfun(@(a) activ(ana_idx{a}(T(ana_idx{a})>SPana_window(1) & T(ana_idx{a})<=SPana_window(2)),:,a),...
                            1:length(ana_idx),'uni',0);%filter out the data out of the ana_window
                        SPamps_idx{idx}{roi,1} = arrayfun(@(a) find(SPamps{idx}{roi,1}{a}),...
                            1:length(ana_idx),'uni',0);%find all the index of amplitude
                        %                         SPamps_idx{idx}{roi,1} = arrayfun(@(a) min(SPamps_idx{idx}{roi,1}{a}),...
                        %                             1:length(ana_idx),'uni',0);%find the first index
                        SPamps{idx}{roi,1} = arrayfun(@(a) SPamps{idx}{roi,1}{a}(SPamps_idx{idx}{roi,1}{a}),...
                            1:length(ana_idx),'uni',0);%first amplitude within the ana_window
                        SPcellmean{idx}{roi} =vertcat(SPamps{idx}{roi}{:});
                    end
                    amp_by_cond{mn}{cond}=vertcat(amp_by_cond{mn}{cond},cellmean{idx});
                    SPamp_by_cond{mn}{cond}=vertcat(SPamp_by_cond{mn}{cond},SPcellmean{idx}');
                end
            end
            SPamp_by_cond{mn}=horzcat(SPamp_by_cond{mn}{:});
            SPamp_by_cond{mn}=arrayfun(@(a) vertcat(SPamp_by_cond{mn}{a,:}),(1:size(SPamp_by_cond{mn},1))','UniformOutput',0);
            SPamp_by_cond{mn}=arrayfun(@(a) nanmean(SPamp_by_cond{mn}{a}),(1:size(SPamp_by_cond{mn},1))');
            amp_by_cond{mn}=horzcat(amp_by_cond{mn}{:});
        end

        % statistics
        % 1-way ANOVA stat
        ydata=[];
        g2=[];
        for mn=1:length(mouse)
            ydata=[ydata;SPamp_by_cond{mn}];
            g2=[g2;repmat(0,length(SPamp_by_cond{mn}),1)];
            for cond=1:size(amp_by_cond{mn},2)
                ydata=[ydata;amp_by_cond{mn}(:,cond)];
                g2=[g2;repmat(StimInten{mn}(cond),length(amp_by_cond{mn}(:,cond)),1)];%condition
            end
            nvalue{mn}=sum(~isnan([SPamp_by_cond{mn} amp_by_cond{mn}]),1);
        end
        [p.anova,tb,stats.anova] = anova1(ydata,g2,'display','off');
        c.anova=multcompare(stats.anova,'Dimension',[1,2],'Display','off');
        g2mat=repmat(1:length(amp_by_cond{mn})+1,2,1);
        groupindex=[repmat(1:2,1,length(amp_by_cond{mn})+1);g2mat(1:end)]';%index of combinations and the corresponding groups
        [R,P]=corrcoef(g2(~isnan(ydata)),ydata(~isnan(ydata)));

        % plot all
        co=lines;
        co(6,:)=[.5 .5 .5];
        cnum=[6 4 5 1];

        figure;hold on;
        x=1:5;
        y=[];
        ysem=[];
        for mn=1:length(mouse)
            y=[y nanmean([SPamp_by_cond{mn} amp_by_cond{mn}],1)'];
            ysem=[ysem nanstd([SPamp_by_cond{mn} amp_by_cond{mn}],1)'./sqrt(length(SPamp_by_cond{mn}))];
            errorbar(x,y(:,mn),ysem(:,mn),'o-',...
                'color',co(cnum(mn),:),...
                'CapSize',20,...
                'MarkerSize',5,...
                'LineWidth',3);
        end

        %significance
        sigX=0.5;
        sigY=0.395;
        formatSpec = '%.3f';
        FS=16;
        if p.anova<0.001
            text(sigX,sigY,'{\it p} < 0.001','FontSize',FS)
        else
            text(sigX,sigY,['{\it p} = ' num2str(p.anova,formatSpec)],'FontSize',FS)
        end

        pos=get(gcf,'position');
        pos(3)=350;
        set(gcf,'Position',pos);
        xticks(1:5)
        xticklabels({'SP' '4' '6' '8' '10'})
        yticks(0:.05:10)
        xlim([0 6])
        ylim([0 .1]+0.3)
        ylabel('Amplitude (\DeltaF/F)')
        xlabel('Airpuff stimulation (psi)')
        set(gca,'FontSize',18)
        set(gcf,'color',[1 1 1])

        % plot significance of anova over stimulus conditions
        Ysig=1;Yincre=0.05;
        gcomb=[1 1;2 2];% combinations of groups to compare
        sig=[];
        for mn=1
            sig{mn}=cell(5);
            idx=(1:5)';%find index of a genotype
            for col=1:5%plot significance of each stimulus intencity
                for row=find(stim~=(stim(col)))
                    if sum(ismember(c.anova(:,[1 2]),[idx(col) idx(row)],'rows'))
                        pvalue=c.anova(ismember(c.anova(:,[1 2]),[idx(col) idx(row)],'rows'),end);
                        if pvalue>=0.05
                            sig{mn}{col,row}=num2str(pvalue,'%1.3f');
                        elseif pvalue<0.05 & pvalue>=0.01
                            sig{mn}{col,row}=[num2str(pvalue,'%1.3f') '*'];
                        elseif pvalue<0.01 & pvalue>=0.001
                            sig{mn}{col,row}=[num2str(pvalue,'%1.3f') '**'];
                        else
                            sig{mn}{col,row}='<0.001***';
                        end
                    end
                end
            end
        end
    end

%% Section 6: plot statistics for normalized amplitude
    function Normalized_First_amp_state(Marker)
        %% ctrl+enter to run this section
        StimInten=cell(size(mouse));
        amp_by_cond=cell(size(mouse));
        for mn=1:length(mouse)
            smoothBC_signal=mouse(mn).smoothBC_signal;
            LOCS=mouse(mn).pksLOCS;
            condlist=mouse(mn).condlist;

            amps=cell(size(condlist));
            SPamps=cell(size(condlist));
            amps_idx=cell(size(condlist));
            SPamps_idx=cell(size(condlist));
            cellmean=cell(size(condlist));
            SPcellmean=cell(size(condlist));
            amp_by_cond{mn}=cell(max(condlist),1);
            SPamp_by_cond{mn}=cell(max(condlist),1);
            Norm_amp_by_cond{mn}=cell(max(condlist),1);
            SPNorm_amp_by_cond{mn}=cell(max(condlist),1);
            time=mouse(mn).time;

            for cond=1:max(condlist)
                %% histogram based on time
                condidx=find(condlist==cond)';
                StimInten{mn}=horzcat(StimInten{mn},mouse(mn).StimInt(condidx(1)));
                %%
                for i=1:length(condidx)
                    % average clls
                    idx=condidx(i);
                    ISinterval=mean(diff(time{idx}));
                    st_idx=data(1).stim(1);
                    T=time{idx}-time{idx}(st_idx+1);

                    for roi=1:size(LOCS{idx},2)
                        activ=smoothBC_signal{idx}(:,roi,:);
                        ana_idx=LOCS{idx}(:,roi); %index of onset or peak depend on analysis

                        amps{idx}{roi,1} = arrayfun(@(a) activ(ana_idx{a}(T(ana_idx{a})>ana_window(1) & T(ana_idx{a})<=ana_window(2)),:,a),...
                            1:length(ana_idx),'uni',0);%filter out the data out of the ana_window
                        amps_idx{idx}{roi,1} = arrayfun(@(a) find(amps{idx}{roi,1}{a}),...
                            1:length(ana_idx),'uni',0);%find all the index of amplitude
                        amps_idx{idx}{roi,1} = arrayfun(@(a) min(amps_idx{idx}{roi,1}{a}),...
                            1:length(ana_idx),'uni',0);%find the first index
                        amps{idx}{roi,1} = arrayfun(@(a) amps{idx}{roi,1}{a}(amps_idx{idx}{roi,1}{a}),...
                            1:length(ana_idx),'uni',0);%first amplitude within the ana_window
                        cellmean{idx}(roi,1) = nanmean(vertcat(amps{idx}{roi}{:}),1);

                        SPamps{idx}{roi,1} = arrayfun(@(a) activ(ana_idx{a}(T(ana_idx{a})>SPana_window(1) & T(ana_idx{a})<=SPana_window(2)),:,a),...
                            1:length(ana_idx),'uni',0);%filter out the data out of the ana_window
                        SPamps_idx{idx}{roi,1} = arrayfun(@(a) find(SPamps{idx}{roi,1}{a}),...
                            1:length(ana_idx),'uni',0);%find all the index of amplitude
                        %                         SPamps_idx{idx}{roi,1} = arrayfun(@(a) min(SPamps_idx{idx}{roi,1}{a}),...
                        %                             1:length(ana_idx),'uni',0);%find the first index
                        SPamps{idx}{roi,1} = arrayfun(@(a) SPamps{idx}{roi,1}{a}(SPamps_idx{idx}{roi,1}{a}),...
                            1:length(ana_idx),'uni',0);%first amplitude within the ana_window
                        SPcellmean{idx}{roi} =vertcat(SPamps{idx}{roi}{:});
                    end
                    amp_by_cond{mn}{cond}=vertcat(amp_by_cond{mn}{cond},cellmean{idx});
                    SPamp_by_cond{mn}{cond}=vertcat(SPamp_by_cond{mn}{cond},SPcellmean{idx}');
                end
            end
            SPamp_by_cond{mn}=horzcat(SPamp_by_cond{mn}{:});
            SPamp_by_cond{mn}=arrayfun(@(a) vertcat(SPamp_by_cond{mn}{a,:}),(1:size(SPamp_by_cond{mn},1))','UniformOutput',0);
            SPamp_by_cond{mn}=arrayfun(@(a) nanmean(SPamp_by_cond{mn}{a}),(1:size(SPamp_by_cond{mn},1))');
            SPNorm_amp_by_cond{mn}=SPamp_by_cond{mn}./SPamp_by_cond{mn};
            Norm_amp_by_cond{mn}=horzcat(amp_by_cond{mn}{:});
            Norm_amp_by_cond{mn}=Norm_amp_by_cond{mn}./SPamp_by_cond{mn};
        end
        
        % statistics
        % 1-way ANOVA stat
        ydata=[];
        g1=[];
        g2=[];
        for mn=1:length(mouse)
            ydata=[ydata;SPNorm_amp_by_cond{mn}];
            g2=[g2;repmat(0,length(SPNorm_amp_by_cond{mn}),1)];
            for cond=1:size(Norm_amp_by_cond{mn},2)
                ydata=[ydata;Norm_amp_by_cond{mn}(:,cond)];
                g2=[g2;repmat(StimInten{mn}(cond),length(Norm_amp_by_cond{mn}(:,cond)),1)];%condition
            end
            nvalue{mn}=sum(~isnan([SPNorm_amp_by_cond{mn} Norm_amp_by_cond{mn}]),1);
        end
        [p.anova,tb,stats.anova] = anova1(ydata,g2,'display','off');
        c.anova=multcompare(stats.anova,'Dimension',[1,2],'Display','off');
        g2mat=repmat(1:length(amp_by_cond{mn})+1,2,1);
        groupindex=[repmat(1:2,1,length(amp_by_cond{mn})+1);g2mat(1:end)]';%index of combinations and the corresponding groups
        [R,P]=corrcoef(g2(~isnan(ydata)),ydata(~isnan(ydata)));

        % plot all
        co=lines;
        co(6,:)=[.5 .5 .5];
        cnum=[6 4 5 1];
        % plot response-stimulus curve
        sigY=1.0;
        L_L_dist=0.017; %distance between significant lines
        L_D_dist=0.005; %distance between significant line and star symbol

        figure;hold on;
        x=1:5;
        y=[];
        ysem=[];
        for mn=1:length(mouse)
            y=[y nanmean([SPNorm_amp_by_cond{mn} Norm_amp_by_cond{mn}],1)'];
            ysem=[ysem nanstd([SPNorm_amp_by_cond{mn} Norm_amp_by_cond{mn}],1)'./sqrt(length(SPNorm_amp_by_cond{mn}))];
            errorbar(x,y(:,mn),ysem(:,mn),'o-',...
                'color',co(cnum(mn),:),...
                'CapSize',20,...
                'MarkerSize',5,...
                'LineWidth',3);
        end

        %significance
        sigX=0.5;
        sigY=1.15;
        dsigY=0.05;
        formatSpec = '%.3f';
        FS=16;
        if p.anova<0.001
            text(sigX,sigY+dsigY*2,'{\it p} < 0.001','FontSize',FS)
        else
            text(sigX,sigY+dsigY*2,['{\it p} = ' num2str(p.anova,formatSpec)],'FontSize',FS)
        end

        pos=get(gcf,'position');
        pos(3)=350;
        set(gcf,'Position',pos);
        xticks(1:5)
        xticklabels({'SP' '4' '6' '8' '10'})
        xlim([0 6])
        yticks(1:.1:10)
        ylim([0.9 1.22]+0.05)
        ylabel('Normalized amplitude')
        xlabel('Airpuff stimulation (psi)')
        set(gca,'FontSize',18)
        set(gcf,'color',[1 1 1])
    end

%% Section 7: categorize stimulus intensity-dependent cells according to normalized amplitude
    function Categorized_Normalized_First_amp_state(Marker)
        %% ctrl+enter to run this section
        StimInten=cell(size(mouse));
        amp_by_cond=cell(size(mouse));
        for mn=1:length(mouse)
            smoothBC_signal=mouse(mn).smoothBC_signal;
            LOCS=mouse(mn).pksLOCS;
            condlist=mouse(mn).condlist;

            amps=cell(size(condlist));
            SPamps=cell(size(condlist));
            amps_idx=cell(size(condlist));
            SPamps_idx=cell(size(condlist));
            cellmean=cell(size(condlist));
            SPcellmean=cell(size(condlist));
            amp_by_cond{mn}=cell(max(condlist),1);
            SPamp_by_cond{mn}=cell(max(condlist),1);
            Norm_amp_by_cond{mn}=cell(max(condlist),1);
            SPNorm_amp_by_cond{mn}=cell(max(condlist),1);
            time=mouse(mn).time;

            for cond=1:max(condlist)
                %% histogram based on time
                condidx=find(condlist==cond)';
                StimInten{mn}=horzcat(StimInten{mn},mouse(mn).StimInt(condidx(1)));
                %%
                for i=1:length(condidx)
                    % average clls
                    idx=condidx(i);
                    ISinterval=mean(diff(time{idx}));
                    st_idx=data(1).stim(1);
                    T=time{idx}-time{idx}(st_idx+1);

                    for roi=1:size(LOCS{idx},2)
                        activ=smoothBC_signal{idx}(:,roi,:);
                        ana_idx=LOCS{idx}(:,roi); %index of onset or peak depend on analysis

                        amps{idx}{roi,1} = arrayfun(@(a) activ(ana_idx{a}(T(ana_idx{a})>ana_window(1) & T(ana_idx{a})<=ana_window(2)),:,a),...
                            1:length(ana_idx),'uni',0);%filter out the data out of the ana_window
                        amps_idx{idx}{roi,1} = arrayfun(@(a) find(amps{idx}{roi,1}{a}),...
                            1:length(ana_idx),'uni',0);%find all the index of amplitude
                        amps_idx{idx}{roi,1} = arrayfun(@(a) min(amps_idx{idx}{roi,1}{a}),...
                            1:length(ana_idx),'uni',0);%find the first index
                        amps{idx}{roi,1} = arrayfun(@(a) amps{idx}{roi,1}{a}(amps_idx{idx}{roi,1}{a}),...
                            1:length(ana_idx),'uni',0);%first amplitude within the ana_window
                        cellmean{idx}(roi,1) = nanmean(vertcat(amps{idx}{roi}{:}),1);

                        SPamps{idx}{roi,1} = arrayfun(@(a) activ(ana_idx{a}(T(ana_idx{a})>SPana_window(1) & T(ana_idx{a})<=SPana_window(2)),:,a),...
                            1:length(ana_idx),'uni',0);%filter out the data out of the ana_window
                        SPamps_idx{idx}{roi,1} = arrayfun(@(a) find(SPamps{idx}{roi,1}{a}),...
                            1:length(ana_idx),'uni',0);%find all the index of amplitude
                        %                         SPamps_idx{idx}{roi,1} = arrayfun(@(a) min(SPamps_idx{idx}{roi,1}{a}),...
                        %                             1:length(ana_idx),'uni',0);%find the first index
                        SPamps{idx}{roi,1} = arrayfun(@(a) SPamps{idx}{roi,1}{a}(SPamps_idx{idx}{roi,1}{a}),...
                            1:length(ana_idx),'uni',0);%first amplitude within the ana_window
                        SPcellmean{idx}{roi} =vertcat(SPamps{idx}{roi}{:});
                    end
                    amp_by_cond{mn}{cond}=vertcat(amp_by_cond{mn}{cond},cellmean{idx});
                    SPamp_by_cond{mn}{cond}=vertcat(SPamp_by_cond{mn}{cond},SPcellmean{idx}');
                end
            end
            SPamp_by_cond{mn}=horzcat(SPamp_by_cond{mn}{:});
            SPamp_by_cond{mn}=arrayfun(@(a) vertcat(SPamp_by_cond{mn}{a,:}),(1:size(SPamp_by_cond{mn},1))','UniformOutput',0);
            SPamp_by_cond{mn}=arrayfun(@(a) nanmean(SPamp_by_cond{mn}{a}),(1:size(SPamp_by_cond{mn},1))');
            SPNorm_amp_by_cond{mn}=SPamp_by_cond{mn}./SPamp_by_cond{mn};
            Norm_amp_by_cond{mn}=horzcat(amp_by_cond{mn}{:});
            Norm_amp_by_cond{mn}=Norm_amp_by_cond{mn}./SPamp_by_cond{mn};
        end

        % plot response-stimulus curve
        sigY=1.0;
        L_L_dist=0.017; %distance between significant lines
        L_D_dist=0.005; %distance between significant line and star symbol
        % plot response-stimulus curve
        sigY=1.0;
        L_L_dist=0.017; %distance between significant lines
        L_D_dist=0.005; %distance between significant line and star symbol

        x=1:5;
        for mn=1:length(mouse)
            [R,P] = arrayfun(@(a) corrcoef(1:5,[SPNorm_amp_by_cond{mn}(a) Norm_amp_by_cond{mn}(a,:)]),(1:length(SPNorm_amp_by_cond{mn}))','uni',0);
            RP{mn} = arrayfun(@(a) [R{a}(2) P{a}(2)],(1:length(SPNorm_amp_by_cond{mn}))','uni',0);
            RP{mn}=vertcat(RP{mn}{:});
            corrindex=RP{mn}(:,1)>0.5 & RP{mn}(:,2)<0.5;

            %correlated
            figure;hold on;
            plot(x',[SPNorm_amp_by_cond{mn}(corrindex) Norm_amp_by_cond{mn}(corrindex,:)]','-',...
                'color',[co(cnum(mn),:) .3],...
                'MarkerSize',5,...
                'LineWidth',1);
            text(0.5,3.5,[num2str(sum(corrindex)) '/' num2str(length(corrindex))],'FontSize',18)
            pos=get(gcf,'position');
            pos(3)=350;
            set(gcf,'Position',pos);
            xticks(1:5)
            xticklabels({'SP' '4' '6' '8' '10'})
            yticks(0:2:10)
            xlim([0 6])
            ylim([0 4])
            ylabel('Normalized amplitude')
            xlabel('Airpuff stimulation (psi)')
            set(gca,'FontSize',18)
            set(gcf,'color',[1 1 1])

            %non-correlated
            figure;hold on;
            plot(x',[SPNorm_amp_by_cond{mn}(~corrindex) Norm_amp_by_cond{mn}(~corrindex,:)]','-',...
                'color',[co(cnum(mn),:) .3],...
                'MarkerSize',5,...
                'LineWidth',1);
            text(0.5,3.5,[num2str(sum(~corrindex)) '/' num2str(length(corrindex))],'FontSize',18)

            pos=get(gcf,'position');
            pos(3)=350;
            set(gcf,'Position',pos);
            xticks(1:5)
            xticklabels({'SP' '4' '6' '8' '10'})
            yticks(0:2:10)
            xlim([0 6])
            ylim([0 4])
            ylabel('Normalized amplitude')
            xlabel('Airpuff stimulation (psi)')
            set(gca,'FontSize',18)
            set(gcf,'color',[1 1 1])
        end
    end

end