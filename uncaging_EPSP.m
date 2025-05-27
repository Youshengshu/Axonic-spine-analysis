clear all;clc
fidx=dir('pre AMPA2 20230829 P64_C57_male_S1C4_2 as_good AS_bad leave*.mat');
GroupName={'pre'};%GroupName={'AS'};
 un_Amptemp=[];
            un_Areatemp=[];
            un_DecayTimetemp=[];
for k=1:length(fidx)
    %% Loading data
    filename=fidx(k).name;
  
    load(filename);
    disp([num2str(k),'/',num2str(length(fidx)),'--',filename])
      FilenameRecord{k,1}=filename;
    %% Reading data
    Marker=[];
    waveformV=[];
    varia=who('*Ch*');
    for i=1:length(varia)
        if isfield(eval(varia{i}),'units')
            str=eval([varia{i},'.units']);
            if ismember('pA',str)
                dataI=eval([varia{i} '.values']);
                dt=eval([varia{i} '.interval']);
                Tstart=eval([varia{i} '.start']);
            end
            if ismember('mV',str)
                dataV=eval([varia{i} '.values']);
            end
        end
        str=eval([varia{i},'.title']);
        if strcmp(str,'imaging')
            Marker=eval([varia{i} '.times']);
            Marker=Marker-Tstart;
        end
    end
    dataI=dataI-mean(dataI(1:10));
    clear(varia{:})
    %% align uncage response
    pretime=0.08;
    postime=0.4;
    Marker_AP=[];Marker_uEPSP=[];
    peaktimeRecord=[];
    Vhold=[];Ihold=[];
    for i=1:length(Marker)
        dataidx=round(Marker(i)/dt-pretime/dt:Marker(i)/dt+postime/dt);
        baselinidx=round(Marker(i)/dt-pretime/dt:Marker(i)/dt);
        waveformV(:,i)=dataV(dataidx);
        Vhold=[Vhold;mean(dataV(baselinidx))];
        Ihold=[Ihold;mean(dataI(baselinidx))];
        peaktime=peakfinder(waveformV(:,i)-Vhold(i),30,25, 1, 0).*dt;
        if ~isempty(peaktime)
            Marker_AP=[Marker_AP,i];
            peaktimeRecord=[peaktimeRecord;peaktime];
        else
            Marker_uEPSP=[Marker_uEPSP,i];
        end
    end
    tspan=linspace(-pretime,postime,round(pretime/dt+postime/dt)+1);
    Range=[-60 -40;-75,-60];%膜电位范围
    figure(1),clf
    if ~isempty(Marker_AP)
            subplot(2,1,1)
    plot(tspan,waveformV(:,Marker_AP))
    end   
     if ~isempty(Marker_uEPSP)
          subplot(2,1,2)
    plot(tspan,waveformV(:,Marker_uEPSP))  
     end
   
     
    %% get parameter
    figure(2),clf 
          
    for ii=1:length(Range)% 对不同膜电位下的EPSP进行average fit
        idx=find(Vhold>Range(ii,1)&Vhold<Range(ii,2));
        idx=intersect(idx,Marker_uEPSP);
          
           
        if ~isempty(idx)
         
            for jj=1:length(Marker_uEPSP)
            waveformI_mean=mean(waveformV(:,idx),2);
            FunEPSP=@(tau,A,B,x) A.*exp(-x./tau)+B;
            foEPSP= fitoptions('Method','NonlinearLeastSquares','MaxIter',10000,...
                'Lower',[0,1,-10],...
                'Upper',[2,10000,1000], ...
                'Startpoint',[0.01,50,0]);%一般参数

            un_Ampidx=find(waveformI_mean==max(waveformI_mean),1,'first');
            ftidix=un_Ampidx:un_Ampidx+round(0.35/dt);
            baselinidx=1:round(pretime/dt);
            un_Amp=nanmean(waveformI_mean(un_Ampidx-20:un_Ampidx+20))-mean(waveformI_mean(baselinidx));
            un_Area=trapz(waveformI_mean-mean(waveformI_mean(baselinidx)));
            slope=diff(waveformI_mean)/(dt*1000);
            slope=[slope(1);slope];
           %max_slopeidx=find(slope(ttspan<endtime)==max(slope(ttspan<endtime)),1,'first');
            max_slope=max(slope);
            maxslopeduridx=un_Ampidx-max(baselinidx);
            maxslopedur=maxslopeduridx*dt*1000;

           %min_slopeidx=find(slope(ttspan<endtime)==min(slope(ttspan<endtime)),1,'first');
            min_slope=min(slope);
            % fitting time constant
            x=tspan(ftidix);
            x=x-x(1);
            y=waveformI_mean(ftidix)-mean(waveformI_mean(baselinidx));
            model=fit(x(:),y(:),FunEPSP,foEPSP);
            yy=feval(model,x);
            un_DecayTime=model.tau;
            paraRecord{1}{k}(ii,:)=[un_DecayTime mean(waveformI_mean(baselinidx))];
     
                un_Amptemp(k,:)=mean(un_Amp,2);
                un_Areatemp(k,:)=mean(un_Area,2);
                un_DecayTimetemp(k,:)=mean(un_DecayTime,2);
                un_maxslope(k,:)=mean(max_slope,2);
                un_maxslopedur(k,:)=mean(maxslopedur,2);
                un_minslope(k,:)=mean(min_slope,2);

            subplot(2,1,2)
            plot(tspan,waveformI_mean),hold on
            plot(tspan(ftidix),yy+mean(waveformI_mean(baselinidx)),'k','linewidth',2)
            plot(tspan(un_Ampidx),waveformI_mean(un_Ampidx),'ko','markerfacecolor','w')
            drawnow
            pause
            box off
            end
      
        end
    end
end           
            epsp_paramRecord{1}(1,:)=un_Amptemp;
            epsp_paramRecord{2}(1,:)=un_Areatemp;
            epsp_paramRecord{3}(1,:)=un_DecayTimetemp;
            epsp_paramRecord{4}(1,:)=un_maxslope;
            epsp_paramRecord{5}(1,:)=un_maxslopedur;
            epsp_paramRecord{6}(1,:)=un_minslope;

epsp_ParamName={'cell','un_Amp', 'un_Area', 'un_DecayTime', 'un_Maxslope', 'un_Maxslopedur', 'un_Minslope'};

save('GroupName',...
    'epsp_paramRecord',...
     'epsp_ParamName','FilenameRecord')
%for kk=1:length(fidx)
    xlswrite('Results_2_photo_epsp_preAS with name.xls',FilenameRecord,['epsp_Param_',GroupName{1}],'A2')
    xlswrite('Results_2_photo_epsp_preAS with name.xls',epsp_ParamName,['epsp_Param_',GroupName{1}],'A1')
    xlswrite('Results_2_photo_epsp_preAS with name.xls',[epsp_paramRecord{1}(1,:)',...
                                             epsp_paramRecord{2}(1,:)',...
                                             epsp_paramRecord{3}(1,:)',...
                                             epsp_paramRecord{4}(1,:)',...
                                             epsp_paramRecord{5}(1,:)',...
                                             epsp_paramRecord{6}(1,:)'],['epsp_Param_',GroupName{1}],'B2')
    
%end
    % get AP threshold
% 
%   
%  spkwaveform_temp=[];
%             HWtemp=[];
%             AHPtemp=[];
%             Thresholdtemp=[];
%             Amptemp=[];
%             max_slopetemp=[];
%             min_slopetemp=[];
%             
%             for jj=1:length(Marker_AP)
%    waveformV_AP=waveformV(:,Marker_AP);
%         [spkwaveform,ttspanspk,HW,AHP,Threshold,Amp,max_slope,min_slope,...
%             HWidx,AHPidx, Thresholdidx,max_slopeidx,min_slopeidx]  = AP_Statistic(waveformV_AP(:,jj),peaktimeRecord(jj),tspan,1);
%                 spkwaveform_temp(jj,:,:)=spkwaveform;
%                 HWtemp(:,jj)=HW;
%                 AHPtemp(:,jj)=AHP;
%                 Thresholdtemp(:,jj)=Threshold;
%                 Amptemp(:,jj)=Amp;
%                 max_slopetemp(:,jj)=max_slope;
%                 min_slopetemp(:,jj)=min_slope;
%                 
%      subplot(2,1,1)
%         plot(tspan,waveformV_AP(:,jj),'displayname',num2str(jj)),hold on
%         plot(ttspanspk(HWidx)+peaktimeRecord(jj)-pretime,spkwaveform(HWidx),'ro','DisplayName',num2str(i))
%         plot(ttspanspk(AHPidx)+peaktimeRecord(jj)-pretime,spkwaveform(AHPidx),'ko')
%         plot(ttspanspk(Thresholdidx)+peaktimeRecord(jj)-pretime,spkwaveform(Thresholdidx),'co')
%         axis tight
%         box off
%              pause
%             end
%             spkparamRecord{1}(:,k)=mean(HWtemp,2);
%             spkparamRecord{2}(:,k)=mean(AHPtemp,2);
%             spkparamRecord{3}(:,k)=mean(Thresholdtemp,2);
%             spkparamRecord{4}(:,k)=mean(Amptemp,2);
%             spkparamRecord{5}(:,k)=mean(max_slopetemp,2);
%             spkparamRecord{6}(:,k)=mean(min_slopetemp,2);
%                              
  
% 
% spkParamName={'HW', 'AHP', 'Threshold', 'Amp','max_slope','min_slope'};
% 
% save('GroupName',...
%     'spkparamRecord',...
%      'spkParamName','FilenameRecord')
% for kk=1:length(fidx)
%     xlswrite('Results_one_photo_AP.xls',FilenameRecord{kk},['spkParam',GroupName{kk}],'A2')
%     xlswrite('Results_one_photo_AP.xls',spkParamName,['spkParam_',GroupName{kk}],'B1')
%     xlswrite('Results_one_photo_AP.xls',[spkparamRecord{1}(1,:)',...
%                                              spkparamRecord{2}(1,:)',...
%                                              spkparamRecord{3}(1,:)',...
%                                              spkparamRecord{4}(1,:)',...
%                                              spkparamRecord{5}(1,:)',...
%                                              spkparamRecord{6}(1,:)'],['spkParam_',GroupName{kk}],'B2')
%     
% end