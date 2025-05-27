function [spkwaveform,ttspan,...
    HW, AHP, Threshold, Amp,max_slope,min_slope,...
    HWidx,AHPidx, Thresholdidx,max_slopeidx,min_slopeidx,latencyidx] = AP_Statistic(data,peaktime,tspan,Goodidx)
%INPUT
% data :   a vector
% tspan:  a vector,have the same elememt with data
% NumSpk: successive number of  spike to be analized
%%
dt=tspan(end)-tspan(end-1);

spktime=peaktime(Goodidx);
%
pretime=0.06;
%之前是0.02
postime=0.05;

spkwaveform=[];

for i=1:length(spktime)
    idx=round(spktime(i)/dt-pretime/dt:spktime(i)/dt+postime/dt);
    ttspan=linspace(-pretime,postime,numel(idx));
    spkwaveform=[spkwaveform,data(idx)];
end


%%


for i=1:length(spktime)
    data=spkwaveform(:,i);
    slope=diff(data)/(dt*1000);
    slope=[slope(1);slope];
    
    %  if length(spktime)>1
    if length(spktime)>1&i<length(peaktime) % change
        endtime=(peaktime(i+1)-peaktime(i));
    else
        endtime=postime*0.75;
    end
    
    max_slopeidx(i)=find(slope(ttspan<endtime)==max(slope(ttspan<endtime)),1,'first');
    max_slope(i)=slope(max_slopeidx(i));
    min_slopeidx(i)=find(slope(ttspan<endtime)==min(slope(ttspan<endtime)),1,'first');
    min_slope(i)=slope(min_slopeidx(i));
    %为什么这里不断地利用指数去找，而不是利用迭代？
    idxA=find(slope(:)>50);% 20 mV /ms
    idxB=find(data(:)>-80);
    idxC=intersect(idxA(:),idxB(:));
%        if ~isempty(idxC)
    Thresholdidx(i)=idxC(1);
    Threshold(i)=data(Thresholdidx(i));
      
    AHPidx(i)=find(data(ttspan>0&ttspan<endtime)==min(data(ttspan>0&ttspan<endtime)),1,'first')+sum(ttspan<0);
    AHP(i)=data(Thresholdidx(i))-data(AHPidx(i));
    
    Amp(i)=max(data)-data(Thresholdidx(i));
   
    HWidx(i,:) = eventedgefinder(data,data(Thresholdidx(i))+Amp(i)/2,1/dt,0.0001,0.0001,1,1);
    HW(i)=(HWidx(i,2)-HWidx(i,1))*dt*1000; 
     latencyidx=round((HWidx(i,2)+HWidx(i,1))/2);
%     latency(i)=latencyidx*dt*1000;
%        else
%           Thresholdidx(i)=nan;
%     Threshold(i)=nan;
%       
%     AHPidx(i)=nan;
%     AHP(i)=nan;
%     
%     Amp(i)=nan;
%     HWidx(i,:) =nan;
%     HW(i)=nan;
%       
%        end
end

end
