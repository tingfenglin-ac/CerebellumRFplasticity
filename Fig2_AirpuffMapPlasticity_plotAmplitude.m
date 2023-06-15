% 1. Add the path of the current script to the search path
% 2. Change the current folder to "Airpuff_map plasticity\amplitude"
% 3. analyze the data within the folder of "Airpuff_map plasticity\amplitude"
% 4. Using ctrl+enter within each function to run each section

%% Section 1: open files
% ctrl+enter to run this section

clear
co=lines;
co(6,:)=[.5 .5 .5];
cnum=[6 2 4 5];
legendlist=[];
colorlist=[];
markersize=30;
criteria=[0 100];
legendlist={'WT control' 'WT tetanus' 'SK2 KO' 'CaMKII TT305/6VA'};
groupnum=4;
clear GOF

% load control data
load('WTcontrolLocation1.mat');
load('WTcontrolLocation2.mat');
conloc=mean(cat(3,location1,location2),3);
con_data=conloc;

% load WT tetanus data
load('WTtetanus.mat');
tet_data=induced;

% load L7SK2 tetanus
load('SK2tetanus.mat');
SK2_data=induced;

% load CaMKII TT3056VA tetanus
load('CaMKIItetanus.mat');
T305V_data=induced;

%% Section 2: plot original calcium amplitude of each genotype before airpuff tetanization
% ctrl+enter to run this section

con_mean=nanmean(con_data,2);
con_sem=nanstd(con_data,0,2)./sqrt(size(con_data,2));
matconncat=[con_data];
concatgroup=ones(size(con_data));

tet_mean=nanmean(tet_data,2);
tet_sem=nanstd(tet_data,0,2)./sqrt(size(tet_data,2));
matconncat=[matconncat tet_data];
concatgroup=[concatgroup 2.*ones(size(tet_data))];

combine_data=cat(2,con_data,tet_data);
comb_mean=nanmean(combine_data,2);
comb_sem=nanstd(combine_data,0,2)./sqrt(size(combine_data,2));
matconncat=[matconncat combine_data];
concatgroup=[concatgroup 0.*ones(size(combine_data))];

SK2_mean=nanmean(SK2_data,2);
SK2_sem=nanstd(SK2_data,0,2)./sqrt(size(SK2_data,2));
matconncat=[matconncat SK2_data];
concatgroup=[concatgroup 3.*ones(size(SK2_data))];

T305V_mean=nanmean(T305V_data,2);
T305V_sem=nanstd(T305V_data,0,2)./sqrt(size(T305V_data,2));
matconncat=[matconncat T305V_data];
concatgroup=[concatgroup 5.*ones(size(T305V_data))];

% 1-way ANOVA stat
data=matconncat(1,concatgroup(1,:)==0|concatgroup(1,:)==3|concatgroup(1,:)==5);
genotype=concatgroup(1,concatgroup(1,:)==0|concatgroup(1,:)==3|concatgroup(1,:)==5);%genotype
[p.anova1,~,stats.anova1] = anova1(data,genotype);
c.anova1=multcompare(stats.anova1);

% plot all groups
cnum=[6 4 5];

y=[comb_mean(1) SK2_mean(1) T305V_mean(1)];
ysem=[comb_sem(1) SK2_sem(1) T305V_sem(1)];
figure;
hold on;

MarkerSize=100;
meanSize=10;
space=0.17;
jrange=0.3;
arrayfun(@(a) errorbar(a+space,y(a),ysem(a),'o-',...
    'color',co(cnum(a),:),...
    'CapSize',30,...
    'MarkerSize',meanSize,...
    'LineWidth',3),1:3,'uni',0);

compgroup=[0 3 5];
arrayfun(@(a) scatter(a-space-0.5*jrange+jrange*rand(1,length(data(genotype==compgroup(a)))),data(genotype==compgroup(a)),MarkerSize,...
    'MarkerFaceColor','flat',...
    'MarkerEdgeColor','flat',...
    'MarkerFaceAlpha',0.2,...
    'MarkerEdgeAlpha',0.2,...
    'CData',co(cnum(a),:)),1:3,'uni',0);

%significance
sigX=1;
sigY=0.04;
dsigY=0.04;
formatSpec = '%.3f';
FS=16;
Ysig=0.85;Yincre=0.07;sigincre=0.008;

sigidx=find(c.anova1(:,end)<0.05 & c.anova1(:,end)>=0.01)';
if sigidx
    for sn=sigidx
        plot(c.anova1(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(c.anova1(sn,1:2)),Ysig+sigincre,'*','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end
sigidx=find(c.anova1(:,end)<0.01 & c.anova1(:,end)>=0.001)';
if sigidx
    for sn=sigidx
        plot(c.anova1(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(c.anova1(sn,1:2)),Ysig+sigincre,'**','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end
sigidx=find(c.anova1(:,end)<0.001)';
if sigidx
    for sn=sigidx
        plot(c.anova1(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(c.anova1(sn,1:2)),Ysig+sigincre,'***','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end

xticks(1:3)
yticks(0:0.5:1)
xticklabels([]);
ylim([0 1])
xlim([0.5 3.5])
ylabel('Amplitude (\DeltaF/F)')
set(gca,'FontSize',23)
set(gcf,'color',[1 1 1])

%% Section 3: plot normalized amplitude (mean+sem) of all genotypes before and after tetanization
% ctrl+enter to run this section

NZcon_data=con_data(:,con_data(1,:)>0);%remove pre is zero
Ncon_data=NZcon_data./NZcon_data(1,:);%normalize
con_mean=nanmean(Ncon_data,2);
con_sem=nanstd(Ncon_data,0,2)./sqrt(size(Ncon_data,2));
matconncat=[Ncon_data];
concatgroup=ones(size(Ncon_data));

NZtet_data=tet_data(:,tet_data(1,:)>0);%remove pre is zero
Ntet_data=NZtet_data./NZtet_data(1,:);%normalize
tet_mean=nanmean(Ntet_data,2);
tet_sem=nanstd(Ntet_data,0,2)./sqrt(size(Ntet_data,2));
matconncat=[matconncat Ntet_data];
concatgroup=[concatgroup 2.*ones(size(Ntet_data))];

NZSK2_data=SK2_data(:,SK2_data(1,:)>0);%remove pre is zero
NSK2_data=NZSK2_data./NZSK2_data(1,:);%normalize
SK2_mean=nanmean(NSK2_data,2);
SK2_sem=nanstd(NSK2_data,0,2)./sqrt(size(NSK2_data,2));
matconncat=[matconncat NSK2_data];
concatgroup=[concatgroup 3.*ones(size(NSK2_data))];

NZT305V_data=T305V_data(:,T305V_data(1,:)>0);%remove pre is zero
NT305V_data=NZT305V_data./NZT305V_data(1,:);%normalize
T305V_mean=nanmean(NT305V_data,2);
T305V_sem=nanstd(NT305V_data,0,2)./sqrt(size(NT305V_data,2));
matconncat=[matconncat NT305V_data];
concatgroup=[concatgroup 4.*ones(size(NT305V_data))];

% repeated ANOVA stat
data=matconncat(1:end);
subject=repmat(1:length(matconncat),3,1);
subject=subject(:)';
g1=repmat(1:3,1,length(matconncat));%time
g2=concatgroup(1:end);%genotype

genotype_cat=categorical(concatgroup(1,:)');
t2 = table(genotype_cat, matconncat(1, :)', matconncat(2, :)', matconncat(3, :)',...
    'VariableNames', {'genotype', 'pre', 'early', 'late'});
Time = [1 2 3]; % 0: pre, 1: during, 2: post
rm2 = fitrm(t2, 'pre-late ~genotype', ...
    'WithinModel', Time, 'WithinModel', 'separatemeans');
ranovatbl =ranova(rm2);
between_sbj=multcompare(rm2,'genotype','By','Time', 'ComparisonType', 'tukey-kramer');
within_sbj=multcompare(rm2,'Time','By','genotype', 'ComparisonType', 'tukey-kramer');



% plot all
cnum=[6 2 4 5];

x=repmat((1:3)',1,groupnum);
y=[con_mean tet_mean SK2_mean T305V_mean];
ysem=[con_sem tet_sem SK2_sem T305V_sem];
figure;subplot(1,3,[2 3]);hold on;
arrayfun(@(a) errorbar((1:3)',y(:,a),ysem(:,a),'o-',...
    'color',co(cnum(a),:),...
    'CapSize',20,...
    'MarkerSize',5,...
    'LineWidth',3),1:groupnum,'uni',0);

%significance
sigX=1;
sigY=1.3;
dsigY=0.03;
formatSpec = '%.3f';
FS=16;
if ranovatbl.pValue(1)<0.001
    text(sigX,sigY+dsigY*2,'{\it p_s_t_i_m} < 0.001','FontSize',FS)
else
    text(sigX,sigY+dsigY*2,['{\it p_s_t_i_m} = ' num2str(ranovatbl.pValue(1),formatSpec)],'FontSize',FS)
end
if ranovatbl.pValue(2)<0.001
    text(sigX,sigY,'{\it p_s_t_i_m _x _g_e_n_e} < 0.001','FontSize',FS)
else
    text(sigX,sigY,['{\it p_s_t_i_m _x _g_e_n_e} = ' num2str(ranovatbl.pValue(2),formatSpec)],'FontSize',FS)
end

% config
xticks(1:4)
xticklabels({'Pre','Early','Late'});
ylim([0.8 1.2])%spon
ylim([0.7 1.4]+0.01)%1stAmp
yticks(0:0.2:2)
xlim([0.8 3.2])
ylabel('Normalized amplitude')
set(gca,'FontSize',23)
set(gcf,'color',[1 1 1])

% plot multiple compbinations
x=repmat((1:3)',1,groupnum);
y=[con_mean tet_mean SK2_mean T305V_mean];
ysem=[con_sem tet_sem SK2_sem T305V_sem];
gcomb=[1 2;2 3;2 4];% combinations of groups to compare
for i=1:size(gcomb,1)
    figure;subplot(1,3,[2 3]);hold on;
    arrayfun(@(a) errorbar((1:3)',y(:,a),ysem(:,a),'o-',...
        'color',co(cnum(a),:),...
        'CapSize',20,...
        'MarkerSize',5,...
        'LineWidth',3),gcomb(i,:),'uni',0);

    % plot significance of anovan over time
    Ysig=1.23;Yincre=0.035;sigincre=0.004;
    for g=1:2%plot significance of each group
        idx=within_sbj.genotype==num2str(gcomb(i,g));%find index of combinations between time
        sig(1,:)=[1 2 within_sbj.pValue(idx & within_sbj.Time_1==1 & within_sbj.Time_2==2)];% x1 x2 and p-value
        sig(2,:)=[1 3 within_sbj.pValue(idx & within_sbj.Time_1==1 & within_sbj.Time_2==3)];% x1 x2 and p-value
        sig(3,:)=[2 3 within_sbj.pValue(idx & within_sbj.Time_1==2 & within_sbj.Time_2==3)];% x1 x2 and p-value

        sigidx=find(sig(:,end)<0.05 & sig(:,end)>=0.01)';
        if sigidx
            for sn=sigidx
                plot(sig(sn,1:2),[Ysig Ysig],'color',co(cnum(gcomb(i,g)),:),'linewidth',3)
                text(mean(sig(sn,1:2)),Ysig+sigincre,'*','FontSize',24,...
                    'HorizontalAlignment','center',...
                    'color',co(cnum(gcomb(i,g)),:));
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(sig(:,end)<0.01 & sig(:,end)>=0.001)';
        if sigidx
            for sn=sigidx
                plot(sig(sn,1:2),[Ysig Ysig],'color',co(cnum(gcomb(i,g)),:),'linewidth',3)
                text(mean(sig(sn,1:2)),Ysig+sigincre,'**','FontSize',24,...
                    'HorizontalAlignment','center',...
                    'color',co(cnum(gcomb(i,g)),:));
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(sig(:,end)<0.001)';
        if sigidx
            for sn=sigidx
                plot(sig(sn,1:2),[Ysig Ysig],'color',co(cnum(gcomb(i,g)),:),'linewidth',3)
                text(mean(sig(sn,1:2)),Ysig+sigincre,'***','FontSize',24,...
                    'HorizontalAlignment','center',...
                    'color',co(cnum(gcomb(i,g)),:));
                Ysig=Ysig+Yincre;
            end
        end
    end

    % plot significance of anova between groups
    for t=1:3%plot significance of each time
        idx=find(between_sbj.genotype_1==num2str(gcomb(i,1)) &...
            between_sbj.genotype_2==num2str(gcomb(i,2)) &...
            between_sbj.Time==t);%find index of combinations between time
        bw_pvalue=between_sbj.pValue(idx);

        if bw_pvalue<0.05 & bw_pvalue>=0.01
            if bw_pvalue~=0
                plot([t t]+0.15,y(t,gcomb(i,:)),'color',[0 0 0],'linewidth',3)
                t=text(t+0.3,mean(y(t,gcomb(i,:))),'*','FontSize',24,...
                    'HorizontalAlignment','center')
                t.Rotation=90;
            end
        end
        if bw_pvalue<0.01 & bw_pvalue>=0.001
            if bw_pvalue~=0
                plot([t t]+0.15,y(t,gcomb(i,:)),'color',[0 0 0],'linewidth',3)
                t=text(t+0.3,mean(y(t,gcomb(i,:))),'**','FontSize',24,...
                    'HorizontalAlignment','center')
                t.Rotation=90;
            end
        end
        if bw_pvalue<0.001
            if bw_pvalue~=0
                plot([t t]+0.15,y(t,gcomb(i,:)),'color',[0 0 0],'linewidth',3)
                t=text(t+0.3,mean(y(t,gcomb(i,:))),'***','FontSize',24,...
                    'HorizontalAlignment','center')
                t.Rotation=90;
            end
        end
    end

    xticks(1:4)
    xticklabels({'Pre','Early','Late'});
    ylim([0.7 1.4]+0.01)%1stAmp
    yticks(0:0.2:2)
    xlim([0.8 3.2])
    ylabel('Normalized amplitude')
    set(gca,'FontSize',23)
    set(gcf,'color',[1 1 1])
end