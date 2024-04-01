%% HANYANG UNIVERSITY - INTERNSHIP
% Josue Perez Sabater - 9087720216

close all;clear;clc

%% Connectivity matrix, constants and templates
simulation_name='sim_RC_2';
load CM          %connectivity matrix M
type=[2 0 3];    %type of stimulus:
%                0: dark
%                1: jumping bar
%                2: rotating bar
%                3: competing bar
%                4: periodic
%                5: turn
T=[500 100 500]; %time span of each type (ms)
Rm=[1 5 10];     %membrane resistance (MOhm)
Cm=[2 .4 .2];    %membrane capacitance (uF)
Rm=8;Cm=2;
V0=-52;          %resting potential (mV)
Vthr=-45;        %threshold potential (mV)
Vmax=20;         %spike peak (mV)
Vmin=-72;        %spike undershoot (mV)
Tap=2;           %AP length (ms)
Ipsc=5;          %PSC amplitude (nA)
Tpsc=5;          %half-life of PSC decay (ms)
Fn=5;            %noise frequency (Hz)
Fs=200;          %sensory input frequency (Hz)
h=.1;            %step size (ms)
time=0:h:sum(T); %time vector (ms)
L=length(time);R=length(Rm);C=length(Cm);

AP=ActPt(Vthr,Vmax,Vmin,Tap,h); %action potential template
OC=OutCr(Ipsc,Tpsc,h);          %output current template

% Noise, stimulus and other inputs
O=zeros(N,L);O(33:50,:)=1; %noise matrix
if any(type==5);O(1:16,:)=1;end
S=STM(type,N,T,h);         %stimulus matrix
[Vin,Iin]=IN(AP,OC,V0,Fn,Fs,h,O,S); %input voltage and current (noise+stimulus)
Iout=zeros(N,L);           %output current from other neurons
Iect=0;                    %ectopic current

% Membrane potential and plot results
V=cell(R,C);
figure('visible','off');
for r=1:length(Rm)
    for c=1:length(Cm)
        [V{r,c}]=param(M,N,color,Rm(r),Cm(c),r,c,R,C,V0,Vthr,Vmax,T,h,time,L,AP,OC,Iin,Iout,Iect);end;end
disp("All simulations finished. Preparing plots...")
fg=findobj('type','fig');set(fg,'color','w')
savefig(figure(1),simulation_name)

%% Main loop
function[V]=param(M,N,color,Rm,Cm,r,c,R,C,V0,Vthr,Vmax,T,h,time,L,AP,OC,Iin,Iout,Iect)
str=[num2str((r-1)*C+c) '/' num2str(R*C)]; %current/total simulations

% Membrane potential
% t: time index
% i: neuron index
% j: 0 if membrane potential, 1 if action potential
% k: index within action potential vector
V(1:N,1)=V0; %initial condition
i=1:N;j=zeros(1,N);k=ones(1,N);
for t=1:L-1
    %clc;disp("Simulation "+str+" in progress. Total time: "+sum(T)+" ms. Time remaining: "+(sum(T)-time(t)-h))
    i0=i(j==0); %neurons without AP
    i1=i(j==1); %neurons with AP
    
    dVdt=@(V) ((V0-V)/Rm+Iin(i0,t)+(Iout(:,t)'*M(:,i0))'+Iect)/Cm;
    V(i0,t+1)=V(i0,t)+dVdt(V(i0,t))*h+normrnd(0,3e-7,size(i0))'; %membrane potential (mV)
    V(i1,t+1)=AP.V(k(i1)); %action potential (mV)
    k(i1)=k(i1)+1;

    APon=i0(V(i0,t+1)>=Vthr); %trigger action potential & output current
    V(APon,t+1)=Vthr;
    m=min(OC.t,L-t);
    Iout(APon,t+1:t+m)=Iout(APon,t+1:t+m)+OC.I(1:m);
    j(APon)=1;

    APoff=i1(k(i1)>AP.t); %return to membrane potential
    j(APoff)=0;k(APoff)=1;
end
clc;disp("Simulation "+str+" finished. Running next simulation...")

% Plot network activity
[Rx,Ry]=find((V==Vmax)');
clr=color(Ry);
Rx=time(Rx');Ry=-Ry'+N+1;
Rx=[Rx;Rx];Ry=[Ry+.5;Ry-.5];

subplot(R,C,(r-1)*C+c)
for i=1:size(Rx,2)
    plot(Rx(:,i),Ry(:,i),clr(i),'linewidth',.7);hold on;end
title("R="+Rm+", C="+Cm);xlim([0 sum(T)]);ylim([.5 N+.5]);
set(gca,'ticklength',[0 0])
end

%% Stimulus - type
function[S]=STM(type,N,T,h)
S=[];

for i=1:length(type)
    time=0:h:T(i);L=length(time);
    Si=zeros(N,L);
    if type(i)==1 %jumping bar
        t=round(1:(L-1)/2:L);
        Si([34 43],t(1):t(2))=1;
        Si([38 47],t(2):t(3))=1;
    
    elseif type(i)==2 %rotating bar
        t=round(1:(L-1)/8:L);
        for j=1:8; Si([32+j 41+j],t(j):t(j+1))=1;end
    
    elseif type(i)==3 %competing bar
        Si([34 43],:)=1;
        Si([38 47],:)=1;
    
    elseif type(i)==4 %periodic
        P=50;D=50;I=[37 46]; %periodicity (ms), duration (ms) and neurons
        PD=repmat([ones(1,D/h) zeros(1,(P-D)/h)],[1 ceil(T/P)]);
        Si(I,:)=repmat([0 PD(1:L-1)],[length(I) 1]);
    
    elseif type(i)==5 %turn
        t=round(1:(L-1)/4:L);
        Si([3 11],t(1):t(2))=1;
        Si(1:8,t(3):t(4))=1;end
    S=[S Si(:,2:end)];end
S=[S(:,1) S];
end

%% Stimulus/Noise - input voltage and current
function[Vin,Iin]=IN(AP,OC,V0,Fn,Fs,h,O,S)
[N,L]=size(O);
Vin=V0*ones(N,L);
Iin=zeros(N,L);

% Action potentials are provided as a Poisson process with rates Fn and Fs
event=O.*(rand(N,L)<Fn*h/1000)+S.*(rand(N,L)<Fs*h/1000); %0 = no AP / 1 = AP

for i=1:N
    ind=find(event(i,:)); %index of APs
    if isempty(ind);else
        ind=ind(logical([1 diff(ind)>AP.t]));end %remove APs that overlap
    
    for t=ind
        m1=min(AP.t,L-t);
        m2=min(OC.t,L-t);
        Vin(i,t:t+m1-1)=AP.V(1:m1);
        Iin(i,t:t+m2-1)=Iin(i,t:t+m2-1)+OC.I(1:m2);
    end;end
end

%% Action potential template
function[AP]=ActPt(Vthr,Vmax,Vmin,Tap,h)
t=0:h:Tap; %timespan (ms)

PDF=normpdf(-1+2*t(t<=Tap/2)/Tap); PDF=(PDF-min(PDF))/(max(PDF)-min(PDF));
SIN=sin((t(t>=Tap/2)-Tap/2)*2*pi/Tap+pi/2); SIN=(SIN-min(SIN))/(max(SIN)-min(SIN));

V1=Vthr+(Vmax-Vthr)*PDF;
V2=Vmin+(Vmax-Vmin)*SIN;
AP.V=[V1(2:end) V2(2:end)];AP.t=length(AP.V); %action potential (mV)
end

%% Output current template
function[OC]=OutCr(Ipsc,Tpsc,h)
t=0:h:2+7*Tpsc; %timespan (ms)

SIN=sin(t(t<=2)*pi/2-pi/2); SIN=(SIN-min(SIN))/(max(SIN)-min(SIN));
EXP=2.^(-(t(t>=2)-2)/Tpsc); EXP=(EXP-min(EXP))/(max(EXP)-min(EXP));

I1=Ipsc*SIN;
I2=Ipsc*EXP;
OC.I=[I1(2:end) I2(2:end)];OC.t=length(OC.I); %output current (nA)
end
