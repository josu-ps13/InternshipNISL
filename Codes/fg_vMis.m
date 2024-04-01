%% NEURAL INTELLIGENCE SYSTEMS LAB - HANYANG UNIVERSITY
% Josue Perez Sabater - Universidad Carlos III de Madrid

close all;clear;clc

%% Stimulus definition
name={'Single neuron','Von Mises distribution'};
a={[],0,[-135 45],0,[45 -135],0,[]}; %angles of each stimulus type (deg)
stim=3; %stimulus to study
T=4; %duration of each stimulus (s)

%% Perform simulation
load CM N epg
V0=-52;
Vthr=-45;
Vmax=20;
Vmin=-72;
Tap=2;
Fn=5;
Fs=120;
h=.1;
time=0:h/1000:T;
L=length(time);
ANG=linspace(-180,180,9)';
AP=ActPt(Vthr,Vmax,Vmin,Tap,h);

figure
for i=0:1
    vmis=i;
    O=ones(N,L); %noise matrix
    S=StiMA(a,vmis,stim,N,epg,T*1000,h,ANG); %stimulus matrix and angle (deg)
    Vin=InpVI(AP,V0,Fn,Fs,h,O,S); %input voltage (noise+stimulus)

    [Rx,Ry]=find((Vin==Vmax)');
    Rx=time(Rx');Ry=-Ry'+N+1;
    Rx=[Rx;Rx];Ry=[Ry+.5;Ry-.5];

    subplot(2,1,i+1);plot(Rx,Ry,'-k','linewidth',.7)
    title(name{i+1})
    xlabel('Time (s)');ylabel('E-PG','rotation',0,'position',[-.15,23,0])
    xlim([0 T]);ylim(abs([epg(end/2) epg(1)]-N-[.5 1.5]))
    set(gca,'ticklength',[0 0],'ytick',[],'fontsize',20)
    set(findobj('type','fig'),'color','w')
end
%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS USED IN THE SCRIPT %%%%%%%%%%%%%%%%%%%%%%

%% Stimulus - matrix and angle
function[S_]=StiMA(a_,vmis,stim,N,epg,T,h,ANG)
S_=[];STM_=[]; %variables with data of all types

for j=1:length(stim) %"STM" is built for each type of stimulus
    L=length(0:h:T(j));
    S=zeros(N,L);    %stimulus matrix
    STM=nan(2,L);    %stimulus angle - second row needed for competing bars
    a=a_{stim(j)+1}; %angle

    if stim(j)==1 %constant bar
        STM(1,:)=a;

    elseif stim(j)==2 %jumping bar
        t=round(1:(L-1)/length(a):L);
        for i=1:length(a);STM(1,t(i)+1:t(i+1))=a(i);end

    elseif stim(j)==3 %rotating bar
        a=a:a+360;
        t=round(1:(L-1)/length(a):L);
        for i=1:length(a);STM(1,t(i)+1:t(i+1))=a(i);end

    elseif stim(j)==4 %competing bars
        STM=repmat(a',[1,L]);

    elseif stim(j)==5 %periodic
        P=200;D=100; %periodicity and duration (ms)
        t=repmat([true(1,D/h) false(1,(P-D)/h)],[1 ceil(T(j)/P)+1]);t=t(1:L);
        STM(1,t)=a;

    elseif stim(j)==6 %turn
        I=3; %initial neuron
        t=round(1:(L-1)/6:L);
        S([I I+8],t(1):t(2))=1;
        S(1:8,t(3):t(5))=1;
        turn=[t(3) t(5)]*h;end %marks start/end of turn

    if j>1;S(:,1)=[];STM(:,1)=[];end %remove first column to avoid overlapping
    S_=[S_ S];
    STM_=[STM_ wrapTo180(STM)];end

% Distribution of the stimulus according to "vmis"
L=size(S_,2);
k=pi*3/4; %kappa
for t=1:L
    vM=(exp(k*cosd(ANG-STM_(:,t)'))-exp(-k))/(exp(k)-exp(-k)); %von Mises distribution
    vM=sum(vM,2,'omitnan'); %sum two stimuli in competing bars
    S_(epg,t)=[vM;vM];end

if vmis==0 %only neuron with highest input is stimulated
    S_(S_==0)=nan;S_=S_==max(S_);end
end

%% Stimulus/Noise - input voltage and current
function[Vin]=InpVI(AP,V0,Fn,Fs,h,O,S)
[N,L]=size(O);
Vin=V0*ones(N,L);
event=(rand(N,L)<O*Fn*h/1000)+(rand(N,L)<S*Fs*h/1000); %0 = no AP / 1 = AP

for i=1:N
    ind=find(event(i,:)); %index of APs
    if isempty(ind)
    else;ind=ind(logical([1 diff(ind)>AP.t]));end %remove APs that overlap
    
    for t=ind
        m1=min(AP.t,L-t); %remove exceeding AP template
        Vin(i,t:t+m1-1)=AP.V(1:m1);end;end
end

%% Action potential template
function[AP]=ActPt(Vthr,Vmax,Vmin,Tap,h)
t=0:h:Tap; %timespan (ms)
PDF=normpdf(-1+2*t(t<=Tap/2)/Tap);
SIN=sin((t(t>=Tap/2)-Tap/2)*2*pi/Tap+pi/2);
V1=Vthr+(Vmax-Vthr)*(PDF-min(PDF))/(max(PDF)-min(PDF));
V2=Vmin+(Vmax-Vmin)*(SIN-min(SIN))/(max(SIN)-min(SIN));
AP.V=[V1(2:end) V2(2:end)];AP.t=length(AP.V); %action potential (mV)
end
