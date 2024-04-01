%% NEURAL INTELLIGENCE SYSTEMS LAB - HANYANG UNIVERSITY
% Josue Perez Sabater - Universidad Carlos III de Madrid

close all;clear;clc

%% Stimulus definition
name={'Constant bar','Jumping bar','Rotating bar',...
    'Competing bars','Periodic','Turn'}; %type 0 1 2 3 4 5 6
a={[],0,[-135 45],0,[45 -135],0,[]}; %angles of each stimulus type (deg)
stim={1 2 3 4 5 6}'; %stimulus to study
T=2; %duration of each stimulus (s)
vmis=0; %(0) default, (1) stimulus provided with von Mises distribution

%% Perform simulation
load CM N pen epg %connectivity matrix M
V0=-52;     %resting potential (mV)
Vthr=-45;   %threshold potential (mV)
Vmax=20;    %spike peak (mV)
Vmin=-72;   %spike undershoot (mV)
Tap=2;      %AP length (ms)
Fn=5;       %noise frequency (Hz)
Fs=120;     %sensory input frequency (Hz)
h=.1;       %step size (ms)
time=0:h/1000:T; %time axis (s)
L=length(time);
ANG=linspace(-180,180,9)'; %angular positions in half hemisphere (deg)
AP=ActPt(Vthr,Vmax,Vmin,Tap,h); %action potential template

SIM=length(stim); %total number of simulations
figure
for sim=1:SIM
    O=ones(N,L); %noise matrix
    S=StiMA(a,vmis,stim{sim},N,epg,T*1000,h,ANG); %stimulus matrix and angle (deg)
    Vin=InpVI(AP,V0,Fn,Fs,h,O,S); %input voltage (noise+stimulus)

    [Rx,Ry]=find((Vin==Vmax)');
    Rx=time(Rx');Ry=-Ry'+N+1;
    Rx=[Rx;Rx];Ry=[Ry+.5;Ry-.5];

    subplot(3,2,sim);plot(Rx,Ry,'-k','linewidth',.7)
    title(name{sim})
    xlabel('Time (s)');ylabel('E-PG','rotation',0,'position',[-.1,23,0])
    xlim([0 T]);ylim(abs([epg(end/2) epg(1)]-N-[.5 1.5]))
    if stim{sim}==6
        ylim(abs([pen(end) pen(1)]-N-[.5 1.5]))
        ylabel('P-EN','rotation',0,'position',[-.1,51,0]);end
    set(gca,'ticklength',[0 0],'ytick',[],'fontsize',15)
end
set(findobj('type','fig'),'color','w')

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
