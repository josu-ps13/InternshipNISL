%% NEURAL INTELLIGENCE SYSTEMS LAB - HANYANG UNIVERSITY
% Josue Perez Sabater - Universidad Carlos III de Madrid

close all;clear;clc

%% Stimulus definition
name={'dark','constant bar','jumping bar','rotating bar',...
    'competing bars','periodic','turn'}; %type 0 1 2 3 4 5 6
a={[],0,[-135 45],0,[45 -135],0,[]}; %angles of each stimulus type (deg)
stim={0 1 2 3 4 5 6}'; %stimulus to study
T={.5 .5 1 1 1 .5 2}'; %duration of each stimulus (s)
NoS=3;  %number of simulations of each stimulus
comb=0; %(0) default, (1) combine types of stimulus in one simulation
vmis=0; %(0) default, (1) stimulus provided with von Mises distribution
indi=0; %(0) default, (1) plot an individual simulation
rast=0; %(0) default, (1) plot raster
tunc=1; %(0) default, (1) plot tuning curves

if tunc==1;stim={1};T={.5};NoS=30;comb=0;vmis=1;
    trial=10; %number of trials
    Tx=linspace(-180,180,NoS+1);Tx(end)=[]; %stimulus angles for tuning curves (deg)
else;trial=1;end
if comb==1
    stim={[0 1 2 3 4 6]};  %combination of stimuli (several independent
    T={[.5 .5 1 1 1 1]};end %combinations can be defined in different rows)

%% Perform simulation
simu=repelem(stim,NoS,1);
T=repelem(T,NoS,1);
SIM=length(simu); %total number of simulations
for sim=1:SIM %each row is one simulation
    for n=1:trial %each column is one trial
        clc;disp("Running simulation "+sim+"/"+SIM+"...")
        if tunc==1;a{2}=Tx(sim);disp("  > trial "+n);end
        data(sim,n)=goSIM(a,simu{sim},T{sim},vmis);end;end %perform simulation
clc;fprintf("All simulations finished! Plotting error...\n\n")

%% Plot RMSD of all simulations (and tuning curves)
if tunc==1;TunCv(Tx,data);end
RMSD=[data(:,end).RMSD];
figure
fprintf("The error of 'dark' and 'turn' simulations is not computable"+...
    "\nbecause there is no stimulus. The RMSD will appear as 'Nan'.\n")
mu=zeros(size(stim)); %RMSD mean of each stimulus type
leg=cell(size(stim)); %name of each stimulus type (legend)
for i=1:length(stim)
    sim=(i-1)*NoS+1:i*NoS; %simulations of same stimulus type
    mu(i)=mean(RMSD(sim));
    bar(sim,RMSD(sim));hold on
    if comb==1;leg{i}="combination "+i; %combination of stimulus types
    else;leg{i}=name{stim{i}+1};end %single type of stimulus
    disp("The average RMSD of "+leg{i}+" is "+num2str(round(mu(i)))+".");end

title('Error in each simulation')
xlabel('Simulation #');ylabel('RMSD (degrees)')
xlim([0 SIM+1]);legend(leg)

%% Plot results of an individual simulation (and raster plots)
% Run this section with indi=1 to choose and plot the results of individual simulations
if indi==1
    sim=input("\nChoose a simulation to plot: ");
    time=data(sim).time;    rh=data(sim).rh;
    time2=data(sim).time2;  HD3=data(sim).HD3;
    time3=data(sim).time3;  E=data(sim).E;
    STM=data(sim).STM;      RMSD=data(sim).RMSD;
    turn=data(sim).turn;    
    disp("The root-mean-square deviation is "+num2str(round(RMSD))+" degrees.")
    disp("Plotting simulation "+sim+" ("+leg{ceil(sim/NoS)}+")...")
    if rast==1;RasPl(sim,data(sim).rast);end

    figure;subplot(311) %plot network activity
    imagesc(rh)
    title("Network activity (simulation "+sim+")")
    ylabel('E-PG neurons')
    set(gca,'xtick',[],'ytick',[])

    subplot(312) %plot stimulus angle and head direction
    plot(time,STM,'b');hold on %stimulus
    plot(time3,HD3,'m') %head direction
    xline(turn/1000,'--b')
    title('Stimulus and head direction angle');
    legend('Stimulus','','Head direction')
    ylabel('Angle (degrees)');set(gca,'ydir','reverse')
    xlim([0 sum(T{sim})]);ylim([-180 180])

    subplot(313) %plot simulation error
    bar(time2,E)
    title("Simulation error (RMSD = "+num2str(round(RMSD))+" degrees)")
    xlabel('Time (s)');ylabel('Angle (degrees)');end

set(findobj('type','fig'),'color','w')

%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS USED IN THE SCRIPT %%%%%%%%%%%%%%%%%%%%%%

%% Perform simulation
function[data]=goSIM(a,stim,T,vmis)
%% Connectivity matrix, parameters and templates
load CM N M c pen peg epg dt7 %connectivity matrix M
Rm=8;       %membrane resistance (MOhm)
Cm=1;       %membrane capacitance (uF)
V0=-52;     %resting potential (mV)
Vthr=-45;   %threshold potential (mV)
Vmax=20;    %spike peak (mV)
Vmin=-72;   %spike undershoot (mV)
Tap=2;      %AP length (ms)
Ipsc=5;     %PSC amplitude (nA)
Tpsc=5;     %half-life of PSC decay (ms)
Fn=5;       %noise frequency (Hz)
Fs=120;     %sensory input frequency (Hz)
h=.1;       %step size (ms)
time=0:h/1000:sum(T); %time axis (s)
le=24;      %boxcar length (ms)
le2=50;     %boxcar2 length (elements)
L=length(time);

ANG=linspace(-180,180,9)';      %angular positions in half hemisphere (deg)
if length(epg)==16;ANG=ANG(1:end-1);end
AP=ActPt(Vthr,Vmax,Vmin,Tap,h); %action potential template
PC=PstCr(Ipsc,Tpsc,h);          %post-synaptic current template

%% Noise, stimulus and other inputs
A=repmat([ANG+11.25;ANG-11.25]*pi/180,[1,L]); %angle matrix (rad)
O=ones(N,L); %noise matrix (to all neurons)
[S,STM,turn]=StiMA(a,vmis,stim,N,epg,T*1000,h,ANG); %stimulus matrix and angle (deg)
[Vin,Iin]=InpVI(AP,PC,V0,Fn,Fs,h,O,S); %input voltage and current (noise+stimulus)
Iout=zeros(N,L); %post-synaptic (output) current from other neurons
Iect=0; %ectopic current

%% Membrane potential
% t: time index
% i: neuron index
% j: 0 if graded potential, 1 if action potential
% k: index within action potential vector
V(1:N,1)=V0; %initial condition
i=1:N;j=zeros(1,N);k=ones(1,N);
for t=1:L-1
    i0=i(j==0); %neurons not firing AP
    i1=i(j==1); %neurons firing AP

    dVdt=@(V) ((V0-V)/Rm+Iin(i0,t)+(Iout(:,t)'*M(:,i0))'+Iect)/Cm;
    V(i0,t+1)=V(i0,t)+dVdt(V(i0,t))*h+normrnd(0,3e-7,size(i0))'; %graded potential (mV)
    V(i1,t+1)=AP.V(k(i1)); %action potential (mV)
    k(i1)=k(i1)+1;

    GP2AP=i0(V(i0,t+1)>=Vthr); %neurons triggering an AP
    V(GP2AP,t+1)=Vthr;
    m=min(PC.t,L-t); %template exceeding time dimension is removed
    Iout(GP2AP,t+1:t+m)=Iout(GP2AP,t+1:t+m)+PC.I(1:m); %post-synaptic (output) current (nA)
    j(GP2AP)=1;

    AP2GP=i1(k(i1)>AP.t); %neurons returning to GP
    j(AP2GP)=0;k(AP2GP)=1;
end

%% Estimate head direction
Vap=V==Vmax; %APs in all neurons
th=Vap(epg,:).*A;th(Vap(epg,:)==0)=nan; %theta: vector angle (AP position, rad)
rh=zeros(N,L); %rho: vector magnitude (AP frequency, Hz)
boxcar=1000*ones(1,le/h)/le;
for i=1:N;rh(i,:)=conv(Vap(i,:),boxcar,'same');end
[x,y]=pol2cart(th,rh(epg,:));
[HD,~]=cart2pol(sum(x,'omitnan'),sum(y,'omitnan')); %head direction (rad)

% Removal of zeros and smoothing
time2=time(HD~=0);HD2=HD(HD~=0); %time and HD without zeros
boxcar2=ones(1,le2)/le2;
[x,y]=pol2cart(HD2,1);
[HD2,~]=cart2pol(conv(x,boxcar2,'same'),conv(y,boxcar2,'same'));
HD2=HD2*180/pi; %head direction (deg)

% Plot line separation
t=find(abs(diff(HD2))>340);
time3=time2;HD3=HD2; %time and HD with "nan"
for i=1:length(t)
    time3=[time3(1:t(i)) nan time3(t(i)+1:end)];
    HD3=[HD3(1:t(i)) nan HD3(t(i)+1:end)];t=t+1;end

%% Compute simulation error
E=reshape(STM([HD;HD]~=0),[2,length(HD2)])-[HD2;HD2]; %error w.r.t. both bars
[~,m]=min(abs(E)); %determines the bar with smaller error
E=E(sub2ind([2,length(m)],m,1:length(m))); %error w.r.t. nearest bar
E=wrapTo180(E); %simulation error (deg)
RMSD=sqrt(mean(E.^2,'omitnan')); %root-mean-square deviation (deg)

data=struct('time',time,'time2',time2,'time3',time3,'STM',STM,...
    'turn',turn,'rh',rh(epg,:),'HD3',HD3,'E',E,'RMSD',RMSD);
data.rast=struct('N',N,'c',c,'pen',pen,'epg',epg,'stim',stim,'T',T,...
    'Vmax',Vmax,'time',time,'Vin',Vin,'V',V);
data.tunc=struct('pen',pen,'peg',peg,'epg',epg,'dt7',dt7,'rh',rh);
end

%% Stimulus - matrix and angle
function[S_,STM_,turn_]=StiMA(a_,vmis,stim,N,epg,T,h,ANG)
S_=[];STM_=[];turn_=[]; %variables with data of all types

for j=1:length(stim) %"STM" is built for each type of stimulus
    L=length(0:h:T(j));
    S=zeros(N,L);    %stimulus matrix
    STM=nan(2,L);    %stimulus angle - second row needed for competing bars
    turn=nan;        %start/end of stimulus "turn"
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
        P=100;D=50; %periodicity and duration (ms)
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
    STM_=[STM_ wrapTo180(STM)];
    turn_=[turn_ turn+sum(T(1:j-1))];end

% Distribution of the stimulus according to "vmis"
L=size(S_,2);
k=pi*3/4; %kappa
for t=1:L
    vM=(exp(k*cosd(ANG-STM_(:,t)'))-exp(-k))/(exp(k)-exp(-k)); %von Mises distribution
    vM=sum(vM,2,'omitnan'); %sum two stimuli in competing bars
    S_(epg,t)=[vM;vM];end

if vmis==0 %only neuron with highest input is stimulated
    S_(S_==0)=nan;S_=S_==max(S_);
    [~,m]=min(abs(STM_(:)'-ANG));
    nans=isnan(STM_(:)');
    STM_(~nans)=ANG(m(~nans)); %stimulus in "STM_" is rounded to values in "ANG"
    STM_=reshape(STM_,[2,L]);end
end

%% Stimulus/Noise - input voltage and current
function[Vin,Iin]=InpVI(AP,PC,V0,Fn,Fs,h,O,S)
[N,L]=size(O);
Vin=V0*ones(N,L);
Iin=zeros(N,L);
event=(rand(N,L)<O*Fn*h/1000)+(rand(N,L)<S*Fs*h/1000); %0 = no AP / 1 = AP

for i=1:N
    ind=find(event(i,:)); %index of APs
    if isempty(ind)
    else;ind=ind(logical([1 diff(ind)>AP.t]));end %remove APs that overlap
    
    for t=ind
        m1=min(AP.t,L-t); %remove exceeding AP template
        m2=min(PC.t,L-t); %remove exceeding PSC template
        Vin(i,t:t+m1-1)=AP.V(1:m1);
        Iin(i,t:t+m2-1)=Iin(i,t:t+m2-1)+PC.I(1:m2);end;end
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

%% Post-synaptic current template
function[PC]=PstCr(Ipsc,Tpsc,h)
t=0:h:2+7*Tpsc; %timespan (ms)
SIN=sin(t(t<=2)*pi/2-pi/2);
EXP=2.^(-(t(t>=2)-2)/Tpsc);
I1=Ipsc*(SIN-min(SIN))/(max(SIN)-min(SIN));
I2=Ipsc*(EXP-min(EXP))/(max(EXP)-min(EXP));
PC.I=[I1(2:end) I2(2:end)];PC.t=length(PC.I); %post-synaptic current (nA)
end

%% Raster plot
function[]=RasPl(sim,data)
disp('Raster plot is being computed. This may take a few moments...')
N=data.N;       T=sum(data.T);
c=data.c;       time=data.time;
epg=data.epg;   Vmax=data.Vmax;
pen=data.pen;   Vin=data.Vin;
stim=data.stim; V=data.V;

[Rx,Ry]=find((Vin==Vmax)'); %input (stimulus+noise) raster
Rx=time(Rx');Ry=-Ry'+N+1;
Rx=[Rx;Rx];Ry=[Ry+.5;Ry-.5];

figure;subplot(211)
plot(Rx,Ry,'-k','linewidth',.7)
title('Stimulus');xlabel('Time (s)')
xlim([0 T]);ylim(abs([epg(end) epg(1)]-N-[.5 1.5]))
if any(stim==6);ylim(abs([pen(end) pen(1)]-N-[.5 1.5]));end
yticks([5.5 19.5 36.5 52.5]);yticklabels({'Pint','E-PG','P-EG','P-EN'})
set(gca,'ticklength',[0 0])

[Rx,Ry]=find((V==Vmax)');%network raster
color=c(Ry);
Rx=time(Rx');Ry=-Ry'+N+1;
Rx=[Rx;Rx];Ry=[Ry+.5;Ry-.5];

subplot(212)
for i=1:size(Rx,2)
    plot(Rx(:,i),Ry(:,i),color(i),'linewidth',.7);hold on;end
title("Network activity (simulation "+sim+")")
xlabel('Time (s)');xlim([0 T]);ylim([.5 N+.5]);
yticks([5.5 19.5 36.5 52.5]);yticklabels({'Pint','E-PG','P-EG','P-EN'})
set(gca,'ticklength',[0 0])
end

%% Tuning curves
function[]=TunCv(Tx,data)
pen=data(1,1).tunc.pen;
peg=data(1,1).tunc.peg;
epg=data(1,1).tunc.epg;
dt7=data(1,1).tunc.dt7;
I=[pen(4) peg(4) epg(5) dt7(3)]; %neurons for tuning curves (one of each type)
name={'P-EN','P-EG','E-PG','Delta7'};

figure
for i=1:4
    Ty=zeros(size(data));
    for x=1:size(data,1) %for all stimulus angles
        for n=1:size(data,2) %for all trials
            Ty(x,n)=mean(data(x,n).tunc.rh(I(i),:)); %firing rate
        end;end
    err=std(Ty,0,2); %error for each angle
    Ty=mean(Ty,2); %firing rate for each angle

    % The plot is shifted to set the peak of the curve at the center
    center=Tx(find(Ty==max(Ty),1));
    if i==4;center=Tx(find(Ty==min(Ty),1));end
    Tx=wrapTo180(Tx-center);
    [Tx,j]=sort(Tx);Ty=Ty(j);
    subplot(2,2,i);errorbar(Tx,Ty,err)
    title("Tuning curve of "+name{i}+" neurons")
    xlabel('Stimulus angle (degrees)');ylabel('Firing rate (Hz)')
    ylim([0 400]);end
end
