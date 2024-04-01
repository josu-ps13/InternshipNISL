%% HANYANG UNIVERSITY - INTERNSHIP
% Josue Perez Sabater - 9087720216

close all;clear;clc
color='krbg'; %colors to use in raster plot for each type of neuron
plotf=1;      %(0) default, (1) plot templates and matrix
insect=1;     %type of insect
% (1) fruit fly - Kakaria
% (2) fruit fly - Pisokas
% (3) locust - Pisokas
% (4) bee - Thomas Stone

%% Connectivity matrix
% row is pre-synaptic neuron, column is post-synaptic neuron

if insect==1 %fruit fly - Kakaria
    N=60;
    Pen_Pen=zeros(16);
    Pen_Peg=zeros(16);
    A=eye(8);A(1,9)=1;A=[flip(fliplr(A));A];Pen_Epg=[A A];
    Pen_Pin=zeros(16,10);

    Peg_Pen=zeros(16);
    Peg_Peg=zeros(16);
    Peg_Epg=Pen_Epg;
    Peg_Pin=zeros(16,10);

    A=eye(8);A(11:18,9:16)=eye(8);Epg_Pen=A;
    A=[zeros(1,16);eye(16);zeros(1,16)];Epg_Peg=A;
    Epg_Epg=zeros(18);
    A=ones(18,10);A(5:end,1)=0;A(1:end-4,end)=0;Epg_Pin=A;

    A=eye(8);A([2,10])=A([10,2]);A(9:10,2)=[1;0];Pin_Pen=[A flip(fliplr(A))];
    A=zeros(10,8);A(9,1)=1;A(2,8)=1;A(3:8,2:7)=eye(6);Pin_Peg=[A A];
    Pin_Epg=zeros(10,18);
    A=ones(10)-eye(10);A(5:end-1,1)=0;A(2:end-4,end)=0;Pin_Pin=A;

    M=[    Pen_Pen     Pen_Peg 5*Pen_Epg     Pen_Pin; %weights * connectivity matrix
           Peg_Pen     Peg_Peg 5*Peg_Epg     Peg_Pin;
        20*Epg_Pen  20*Epg_Peg   Epg_Epg  20*Epg_Pin;
       -20*Pin_Pen -20*Pin_Peg   Pin_Epg -20*Pin_Pin];

    pen=1:16;peg=17:32;epg=33:50;dt7=51:60;
    c=repelem(color,[16 16 18 10]);

elseif insect==2 %fruit fly - Pisokas
    N=60;
    Pen_Pen=zeros(16);
    Pen_Peg=zeros(16,18);
    A=eye(8);A(1,9)=1;A=[flip(fliplr(A));A];Pen_Epg=[A A];
    Pen_Dt7=zeros(16,8);

    Peg_Pen=zeros(18,16);
    Peg_Peg=zeros(18);
    A1=eye(9);A1(1,end)=1;A2=A1;A2([end 1],[1 1])=[1 1;0 0];A=[A1;A2];Peg_Epg=[A flip(fliplr(A))];
    Peg_Dt7=zeros(18,8);

    A=eye(8);A(11:18,9:16)=eye(8);Epg_Pen=A;
    Epg_Peg=eye(18);
    Epg_Epg=zeros(18);
    A=~eye(8);Epg_Dt7=[A;A;A(1:2,:)];

    A=eye(8);Dt7_Pen=[A [A(7:8,:);A(1:6,:)]];
    A=eye(8);A=[A A A];Dt7_Peg=A(:,1:18);
    Dt7_Epg=zeros(8,18);
    Dt7_Dt7=~eye(8);

    M= [    Pen_Pen     Pen_Peg 5*Pen_Epg     Pen_Dt7; %weights * connectivity matrix
            Peg_Pen     Peg_Peg 5*Peg_Epg     Peg_Dt7;
         30*Epg_Pen  30*Epg_Peg   Epg_Epg  30*Epg_Dt7;
        -30*Dt7_Pen -30*Dt7_Peg   Dt7_Epg -10*Dt7_Dt7];

    pen=1:16;peg=17:34;epg=35:52;dt7=53:60;
    c=repelem(color,[16 18 18 8]);

elseif insect==3 %locust - Pisokas
elseif insect==4 %bee - Thomas Stone
    TL_CL1=eye(16);

    A=eye(8);CL1_TB1=[A;A];

    TB1_TB1
    A=eye(8);TB1_CP4=[A A];
    A=eye(8);A=[A(:,2:end) A(:,1)];TB1_CP1=[A flip(fliplr(A))];

    CP4_Pon=eye(16);
    A=eye(8);CP4_CP1=[zeros(8) A;A zeros(8)];

    A=[ones(1,8);zeros(1,8)];TN2_CP4=[A;flip(A)];

    A=eye(8);A=[A(5:8,:);A(1:4,:)];Pon_CP1=[zeros(8) A;A zeros(8)];

end

save CM N M c pen peg epg dt7
disp('New data has been saved.')

%% Parameters
Vthr=-45; %threshold potential (mV)
Vmax=20;  %spike peak (mV)
Vmin=-72; %spike undershoot (mV)
Ipsc=5;   %PSC amplitude (nA)
Tap=2;    %AP length (ms)
Tpsc=5;   %half-life of PSC decay (ms)
h=.01;    %step size (ms)

%% Action potential
t=0:h:Tap; %timespan (ms)

PDF=normpdf(-1+2*t(t<=Tap/2)/Tap);
SIN=sin((t(t>=Tap/2)-Tap/2)*2*pi/Tap+pi/2);

V1=Vthr+(Vmax-Vthr)*(PDF-min(PDF))/(max(PDF)-min(PDF));
V2=Vmin+(Vmax-Vmin)*(SIN-min(SIN))/(max(SIN)-min(SIN));
AP.t=t;AP.V=[V1 V2(2:end)]; %action potential (mV)

%% Post-synaptic current
t=0:h:2+7*Tpsc; %timespan (ms)

SIN=sin(t(t<=2)*pi/2-pi/2);
EXP=2.^(-(t(t>=2)-2)/Tpsc);

I1=Ipsc*(SIN-min(SIN))/(max(SIN)-min(SIN));
I2=Ipsc*(EXP-min(EXP))/(max(EXP)-min(EXP));
PC.t=t;PC.I=[I1 I2(2:end)]; %post-synaptic current (nA)

%% Plot templates and connectivity matrix
if plotf==1
    figure;heatmap(M,'XDisplayLabel',nan(N,1),'YDisplayLabel',nan(N,1))
    cm=[linspace(.60,1,100) linspace(1,.47,100);
        linspace(.35,1,100) linspace(1,.67,100);
        linspace(.75,1,100) linspace(1,.19,100)]';
    colormap(cm);grid off
    title("Connectivity matrix ("+N+" neurons)")
    xlabel('Post-synaptic neuron');ylabel('Pre-synaptic neuron')

    figure;plot(AP.t,AP.V,'LineWidth',2,'Color',[.6,.35,.75]);axis tight
    yline(Vthr,'--','LabelHorizontalAlignment','center','Color',[.65,.65,.65])
    annotation('textbox', [0, 0.5, 0, 0], 'string','Threshold potential (-45 mV)',...
        'EdgeColor','none','Position',[.36,.41,.33,0])
    title('Action potential template')
    xlabel('Time (ms)');ylabel('Voltage (mV)')

    figure;plot(PC.t,PC.I,'LineWidth',2,'Color',[.6,.35,.75]);axis tight
    yline(Ipsc/2,'--','Color',[.65,.65,.65])
    xline([2 6.944],'--','Color',[.65,.65,.65])
    annotation('doublearrow',[0.174 0.275],[.4 .4],'color',[.65,.65,.65],...
        'Head1Length',7,'Head1Width',7,'Head2Length',7,'Head2Width',7)
    annotation('textbox', [0, 0.5, 0, 0], 'string', 'Half-life (5 ms)',...
        'EdgeColor','none','BackgroundColor','w','Position',[.15,.3,.18,.07])
    title('Post-synaptic current template')
    xlabel('Time (ms)');ylabel('Current (nA)')
    set(findobj('type','axe'),'TickDir','out')
    set(findobj('type','fig'),'color','w');end
