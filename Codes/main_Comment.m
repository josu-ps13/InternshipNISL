%% NEURAL INTELLIGENCE SYSTEMS LAB - HANYANG UNIVERSITY
% Josue Perez Sabater - Universidad Carlos III de Madrid

% Code explanation
% The script "main.m" performs a simulation of the CX of an insect. The
% stimulus is defined at the beginning, with options to deliver complex and
% numerous stimuli. Then the simulation is performed using the connectivity
% matrix of an insect, loaded in the variable "CM.mat". Euler's method is
% used to develop the membrane potentials along time. The resulting network
% activity is processed to obtain the head direction of the insect, which
% is compared to the original stimulus and the error is computed. Once the
% simulation has concluded, a figure shows the root-mean-square deviation
% of each simulation and a prompt in the command window allows to choose
% one specific simulation to further study it. Additional plots may be
% displayed.

close all;clear;clc

%% Stimulus definition
% Variables "name" and "a" contain the names and angles of each stimulus
% type. The type of stimulus provided to the network is defined by "stim".
% 0 -> dark: only random noise (also added to all other stimulus)
% 1 -> constant bar: stimulation of one single E-PG neuron
% 2 -> jumping bar: stimulation of one E-PG neuron, followed by another non
%      adjacent E-PG neuron
% 3 -> rotating bar: gradual stimulation of all E-PG neurons
% 4 -> competing bars: stimulation of two, non-adjacent E-PG neurons at the
%      same time
% 5 -> periodic: intermittent stimulation of one E-PG neuron
% 6 -> turn: stimulation of one P-EN neuron, followed by all P-EN neurons
%      in one hemisphere (represents a turn in the heading direction)
%
% The angle of the stimulus can be modified in "a":
% - 'dark' and 'turn' -> no angle since there is no visual stimulus
% - 'constant bar' and 'periodic bar' -> delivered only at one angle
% - 'competing bars' -> delivered at two angles
% - 'jumping bar' -> delivered at many angles
% - 'rotating bar' -> only the starting angle can be set
% 
% The duration of the stimulus "T" is also defined. It may be different for
% each stimulus. Setting the parameter "NoS" (number of simulations) allows
% to perform each simulation several times. Additional variables are set:
% - comb=1 -> the stimulus will be provided as a combination of different
%             types, each of them with their own duration
% - vmis=1 -> the stimulus is provided to several neurons following a von
%             Mises distribution. Otherwise, only one neuron receives it
% - indi=1 -> a prompt in the command window allows to choose one single
%             simulation and plot its results individually
% - rast=1 -> if an individual simulation is plotted, an additional figure
%             will appear with the raster plots of the noise + stimulus
%             and the network activity
% - tunc=1 -> the stimulus is adapted to allow the plotting of tuning
%             curves (several trials at each angle are performed to improve
%             accuracy)
name={'dark','constant bar','jumping bar','rotating bar',...
    'competing bars','periodic','turn'}; %type 0 1 2 3 4 5 6
a={[],0,[-135 45],0,[45 -135],0,[]}; %angles of each stimulus type (deg)
stim={0 1 2 3 4 5 6}'; %stimulus to study
T={.5 1 1 2 2 1 2}'; %duration of each stimulus (s)
NoS=3;  %number of simulations of each stimulus
comb=0; %(0) default, (1) combine types of stimulus in one simulation
vmis=0; %(0) default, (1) stimulus provided with von Mises distribution
indi=0; %(0) default, (1) plot an individual simulation
rast=0; %(0) default, (1) plot raster
tunc=0; %(0) default, (1) plot tuning curves

if tunc==1;stim={1};T={.5};NoS=30;comb=0;vmis=1;
    trial=20; %number of trials
    Tx=linspace(-180,180,NoS+1);Tx(end)=[]; %stimulus angles for tuning curves (deg)
else;trial=1;end
if comb==1
    stim={[0 2 3 4 6]};  %combination of stimuli (several independent
    T={[.5 1 1 1 2]};end %combinations can be defined in different rows)

%% Perform simulation
% Variable "simu" contains all simulations. A loop runs through all of them
% and through all trials, storing all the relevant information in the
% matrix-structure "data". If tuning curves are being plotted, then each
% simulation uses a different stimulus angle, taken from "Tx".
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
% A figure shows the RMSD of all simulations. Simulations may represent one
% single type of stimulus or a combination of types, depending on "comb".
% If a simulation contains 'dark' or 'turn' stimulus, the RMSD of those
% sections will be ignored (equal to "nan"). The colors of the bars put
% together the simulations with same stimulus since they are repeated "NoS"
% times. For each group, the mean RMSD "mu" is computed and displayed in
% the command window. Additionally, the tuning curves may be plotted.
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
% Run this section with indi=1 to choose and plot the results of individual
% simulations. The results are presented in three plots:
% - The activity of all E-PG neurons (shows the frequency and position of
%   APs along time)
% - The stimulus angle and the head direction
% - The simulation error, when available
% The x axis in all these graphs is the time (in s). The y axis in second
% and third plots is the angle (in degrees).
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
% This function "goSIM" receives the information about a stimulus and
% simulates the response of the insect CX. All relevant data is saved in
% the variable "data".

%% Connectivity matrix, parameters and templates
% The number of neurons "N" in the network and the connectivity matrix "M"
% are loaded from a file called "CM.mat", as well as the colors for the
% raster plot "c" and the indices of all neurons. Then, the necessary
% parameters are defined and the templates for the action potential "AP"
% and the post-synaptic current "PC" are constructed calling their
% respective functions.
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
AP=ActPt(Vthr,Vmax,Vmin,Tap,h); %action potential template
PC=PstCr(Ipsc,Tpsc,h);          %post-synaptic current template

%% Noise, stimulus and other inputs
% The angle "A", noise "O" and stimulus "S" matrices are constructed. These
% matrices are size N*L (row = neuron, column = time instant). "A" contains
% the angle represented by each neuron while "O" and "S" contain 0 when the
% noise/stimulus doesn't exist, and 1 (or other value >0) when it does. The
% variable "STM" contains the angle of the stimulus in degrees for eac
% time instant. There is only one angle for all types of stimulus except
% for 'dark' and 'turn' (there is no visual stimulus) and 'competing bars'
% (there are two stimulus). Whenever there is no stimulus, the value of
% "STM" is "nan".
% 
% The equation for the membrane potential receives three input currents:
% - Iin:  produced by noise & stimulus, can be defined before the main loop
% - Iout: produced by other neurons in the CX, defined during main loop
% - Iect: ectopic current, constant
A=repmat([ANG+11.25;ANG-11.25]*pi/180,[1,L]); %angle matrix (rad)
O=ones(N,L); %noise matrix (to all neurons)
[S,STM,turn]=StiMA(a,vmis,stim,N,epg,T*1000,h,ANG); %stimulus matrix and angle (deg)
[Vin,Iin]=InpVI(AP,PC,V0,Fn,Fs,h,O,S); %input voltage and current (noise+stimulus)
Iout=zeros(N,L); %post-synaptic (output) current from other neurons
Iect=0; %ectopic current

%% Membrane potential
% The main loop moves along all time instants and calculates the membrane
% potential "V", which may be defined either by the equation of a graded
% potential (using Euler's method) or by an action potential (inserting the
% template). In each iteration, neurons are classified as 'firing AP' or
% 'not firing AP'. When the loop ends, the matrix "V" contains the value of
% membrane potential for all neurons at all time instants.

% Index explanation
% t: time index (from 1 to L)
% i: neuron index (from 1 to N)
% j: determines if neurons are currently firing an AP (j=1) or not (j=0).
%    All neurons start with j=0 (graded potential) and they trigger an AP
%    when their voltage reaches "Vthr". When the AP ends, they return to
%    graded potential.
% k: for neurons firing an AP, tracks the progress of the neuron within the
%    AP vector. All neurons start with k=1. When a neuron starts firing an
%    AP, k starts counting as time moves forward. When k reaches the AP
%    length (k=AP.t), it means that the whole AP template has been inserted
%    and therefore the neuron must return to graded potential.
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
% The head direction "HD" is encoded by an 'activity bump' in E-PG neurons
% (from 33 to 50 within "M"). APs occur when V=Vmax. For a given time
% instant, there may be several neurons firing AP. To obtain one single
% angle (heading direction), a vector summation is performed where the AP
% position is the vector angle (theta) and the frequency of APs is the
% vector magnitude (rho). The angular position (theta) of all APs is
% obtained from the angle matrix "A". The AP frequency (rho) of each neuron
% along time is computed with a convolution. Later, vector summation of the
% APs at each time instant is done in cartesian coord. to obtain one single
% head direction per instant.
Vap=V==Vmax; %APs in all neurons
th=Vap(epg,:).*A;th(Vap(epg,:)==0)=nan; %theta: vector angle (AP position, rad)
rh=zeros(N,L); %rho: vector magnitude (AP frequency, Hz)
boxcar=1000*ones(1,le/h)/le;
for i=1:N;rh(i,:)=conv(Vap(i,:),boxcar,'same');end
[x,y]=pol2cart(th,rh(epg,:));
[HD,~]=cart2pol(sum(x,'omitnan'),sum(y,'omitnan')); %head direction (rad)

% Removal of zeros and smoothing
% Zeros are removed since they are due to the abscence of AP rather than AP
% firing at an angle of 0. HD signal is smoothen with a second convolution.
time2=time(HD~=0);HD2=HD(HD~=0); %time and HD without zeros
boxcar2=ones(1,le2)/le2;
[x,y]=pol2cart(HD2,1);
[HD2,~]=cart2pol(conv(x,boxcar2,'same'),conv(y,boxcar2,'same'));
HD2=HD2*180/pi; %head direction (deg)

% Plot line separation
% For plotting purposes, when head direction jumps from -180 to 180, "nan"
% is added to the vector so the plot line appears separated in the graph
t=find(abs(diff(HD2))>340);
time3=time2;HD3=HD2; %time and HD with "nan"
for i=1:length(t)
    time3=[time3(1:t(i)) nan time3(t(i)+1:end)];
    HD3=[HD3(1:t(i)) nan HD3(t(i)+1:end)];t=t+1;end

%% Compute simulation error
% The simulation error "E" is computed as the difference between the HD and
% the stimulus (STM-HD). The RMSD of the whole simulation is also computed
% and it gives an idea of the overall error of the simulation. These values
% are not computable when the stimulus type is 'dark' or 'turn' because
% there is no stimulus to which compare the HD. The type 'periodic' also
% has no stimulus at some time instants. The value of "E" when it is not
% computable is "nan". There are two stimulus at the same time in
% 'competing bars' and only the nearest bar (with smaller error) is used to
% compare. For other types of stimulus, the second row is ignored ("nan").
E=reshape(STM([HD;HD]~=0),[2,length(HD2)])-[HD2;HD2]; %error w.r.t. both bars
[~,m]=min(abs(E)); %determines the bar with smaller error
E=E(sub2ind([2,length(m)],m,1:length(m))); %error w.r.t. nearest bar

E=wrapTo180(E); %simulation error (deg)
RMSD=sqrt(mean(E.^2,'omitnan')); %root-mean-square deviation (deg)

data=struct('time',time,'time2',time2,'time3',time3,'STM',STM,...
    'turn',turn,'rh',rh(epg,:),'HD3',HD3,'E',E,'RMSD',RMSD);
data.rast=struct('N',N,'c',c,'epg',epg,'stim',stim,'T',T,...
    'Vmax',Vmax,'time',time,'Vin',Vin,'V',V);
data.tunc=struct('pen',pen,'peg',peg,'epg',epg,'dt7',dt7,'rh',rh);
end

%% Stimulus - matrix and angle
function[S_,STM_,turn_]=StiMA(a_,vmis,stim,N,epg,T,h,ANG)
% This function "StiMA" builds the stimulus matrix "S" and angle "STM". It
% receives as inputs the type of stimulus "stim", matrix size "N", time
% limit "T" and angular positions "ANG" represented by the neurons in
% degrees. As a general process to build the stimulus, the time dimension
% is divided in intervals where different neurons are set to be firing. The
% variable "STM" has two rows containing one angular value for each time
% instant, except for 'dark' and 'turn' (there is no visual stimulus) and
% 'competing bars' (there are two stimulus). Whenever there is no stimulus,
% the value of "STM" is set to "nan". The second row is only used when the
% stimulus is 'competing bars', otherwise its value is "nan". In the case
% of the 'turn' stimulus, an additional variable "turn" is created with
% the times of start and end of the stimulus, in order to indicate it in
% the final plot.

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
% The variable "STM" (or "STM_") determines which angle will receive the
% most intense activity. There are two modes:
% - vmis=0: only one neuron receives the stimulus. If the stimulus
%   coincides with the angular position of a neuron, then that neuron will
%   receive the stimulus. If the stimulus is between two angular positions,
%   then the nearest neuron will receive the stimulus. If the stimulus is
%   exactly between two angular positions, then both neurons will receive
%   the stimulus.
% - vmis=1: several neurons receive the stimulus with different intensities
%   following a von Mises distribution. "STM" determines the peak of the
%   distribution. If the peak is between two angular positions, there's no
%   need of approximation. Instead, the intensities will be defined by the
%   distribution "vM". When 'competing bars', two stimuli lead to two
%   distributions that are summed together. The distribution is the same
%   for the two hemispheres.
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
% This function "InpVI" returns the input voltage and current caused by the
% noise and stimulus. It receives as inputs the AP and PSC templates, the
% initial voltage "V0", the noise "Fn" and stimulus "Fs" frequencies, and
% the noise "O" and stimulus "S" matrices. Action potentials are provided
% as a Poisson process with rates Fn and Fs. The matrix "event" combines
% noise and stimulus to provide a matrix with 1 when APs are fired and 0
% otherwise. Then, the overlapping APs are removed and the templates are
% inserted in the input voltage "Vin" and input current "Iin" matrices. If
% a template exceeds the time dimension, the exceeding part is removed.

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
% This function "ActPt" builds the action potential template. The input
% variables are the threshold "Vthr", maximum "Vmax" and minimum "Vmin"
% voltages, the AP length "Tap" and the step size "h". The template
% consists of an increasing interval "SIN" followed by a decreasing
% interval "EXP". The output variable contains the template "AP.V" and its
% length "AP.t".

t=0:h:Tap; %timespan (ms)
PDF=normpdf(-1+2*t(t<=Tap/2)/Tap);
SIN=sin((t(t>=Tap/2)-Tap/2)*2*pi/Tap+pi/2);
V1=Vthr+(Vmax-Vthr)*(PDF-min(PDF))/(max(PDF)-min(PDF));
V2=Vmin+(Vmax-Vmin)*(SIN-min(SIN))/(max(SIN)-min(SIN));
AP.V=[V1(2:end) V2(2:end)];AP.t=length(AP.V); %action potential (mV)
end

%% Post-synaptic current template
function[PC]=PstCr(Ipsc,Tpsc,h)
% This function "PstCr" builds the post-synaptic current template. The
% input variables are the PSC amplitude "Ipsc", half-life of PSC decay
% "Tpsc" and the step size "h". The template consists of an increasing
% interval "SIN" followed by a decreasing interval "EXP". The output
% variable contains the template "PC.I" and its length "PC.t".

t=0:h:2+7*Tpsc; %timespan (ms)
SIN=sin(t(t<=2)*pi/2-pi/2);
EXP=2.^(-(t(t>=2)-2)/Tpsc);
I1=Ipsc*(SIN-min(SIN))/(max(SIN)-min(SIN));
I2=Ipsc*(EXP-min(EXP))/(max(EXP)-min(EXP));
PC.I=[I1(2:end) I2(2:end)];PC.t=length(PC.I); %post-synaptic current (nA)
end

%% Raster plot
function[]=RasPl(sim,data)
% This function "RasPl" displays the raster plot of the stimulus to
% the E-PG neurons (and P-EN in 'turn') and of the whole network
% activity, visualizing the APs of all neurons along time.

disp('Raster plot is being computed. This may take a few moments...')
N=data.N;       time=data.time;
c=data.c;       Vmax=data.Vmax;
epg=data.epg;   Vin=data.Vin;
stim=data.stim; V=data.V;
T=sum(data.T);

[Rx,Ry]=find((Vin==Vmax)'); %input (stimulus+noise) raster
Rx=time(Rx');Ry=-Ry'+N+1;
Rx=[Rx;Rx];Ry=[Ry+.5;Ry-.5];

figure;subplot(211)
plot(Rx,Ry,'-k','linewidth',.7)
title('Stimulus');xlabel('Time (s)')
xlim([0 T]);ylim(abs([epg(end) epg(1)]-N-[.5 1.5]))
if any(stim==6);ylim([10.5 60.5]);end
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
% This function "TunCv" displays the tuning curves of all neuron types. To
% do so, the average activity of one neuron of each type "I" is obtained
% for several stimulus angles (defined by "Tx"). The input variable "data"
% contains the simulations for all those angles (rows) and for several
% trials (columns). Therefore, for a given angle, the average of all trials
% is computed.

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
    % Therefore, the plot will always be centered, no matter which neurons
    % are defined in "I".
    center=Tx(find(Ty==max(Ty),1));
    if i==4;center=Tx(find(Ty==min(Ty),1));end
    Tx=wrapTo180(Tx-center);
    [Tx,j]=sort(Tx);Ty=Ty(j);
    subplot(2,2,i);errorbar(Tx,Ty,err)
    title("Tuning curve of "+name{i}+" neurons")
    xlabel('Stimulus angle (degrees)');ylabel('Firing rate (Hz)')
    ylim([0 400]);end
end
