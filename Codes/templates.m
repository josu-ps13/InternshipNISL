%% HANYANG UNIVERSITY - INTERNSHIP
% Josue Perez Sabater - 9087720216

close all;clear;clc
color='krbg'; %colors to use in raster plot for each type of neuron
plotf=0;      %(0) default, (1) plot templates and matrix
insect=1;     %type of insect
% (1) fruit fly - Kakaria
% (2) fruit fly - Pisokas
% (3) locust - Pisokas
% (4) bee - Thomas Stone

%% Connectivity matrix
% row is pre-synaptic neuron, column is post-synaptic neuron

if insect==1 %fruit fly - Kakaria
    N=60;
    pen_pen=zeros(16);
    pen_peg=zeros(16);
    A=eye(8);A(1,9)=1;A=[rot90(A,2);A];pen_epg=[A A];
    Pen_Pin=zeros(16,10);

    peg_pen=zeros(16);
    peg_peg=zeros(16);
    peg_epg=pen_epg;
    Peg_Pin=zeros(16,10);

    A=eye(8);A(11:18,9:16)=eye(8);epg_pen=A;
    A=[zeros(1,16);eye(16);zeros(1,16)];epg_peg=A;
    epg_epg=zeros(18);
    A=ones(18,10);A(5:end,1)=0;A(1:end-4,end)=0;Epg_Pin=A;

    A=eye(8);A([2,10])=A([10,2]);A(9:10,2)=[1;0];Pin_Pen=[A rot90(A,2)];
    A=zeros(10,8);A(9,1)=1;A(2,8)=1;A(3:8,2:7)=eye(6);Pin_Peg=[A A];
    Pin_Epg=zeros(10,18);
    A=ones(10)-eye(10);A(5:end-1,1)=0;A(2:end-4,end)=0;Pin_Pin=A;

    M=[    pen_pen     pen_peg 5*pen_epg     Pen_Pin; %weights * connectivity matrix
           peg_pen     peg_peg 5*peg_epg     Peg_Pin;
        20*epg_pen  20*epg_peg   epg_epg  20*Epg_Pin;
       -20*Pin_Pen -20*Pin_Peg   Pin_Epg -20*Pin_Pin];

    pen=1:16;peg=17:32;epg=33:50;dt7=51:60;
    c=repelem(color,[16 16 18 10]);

elseif insect==2 %fruit fly - Pisokas
    N=60;
    pen_pen=zeros(16);
    pen_peg=zeros(16,18);
    A=eye(8);A(1,9)=1;A=[rot90(A,2);A];pen_epg=[A A];
    pen_dt7=zeros(16,8);

    peg_pen=zeros(18,16);
    peg_peg=zeros(18);
    A1=eye(9);A1(1,end)=1;A2=A1;A2([1 9])=[0 1];A=[A1;A2];peg_epg=[A rot90(A,2)];
    peg_dt7=zeros(18,8);

    A=eye(8);A(11:18,9:16)=eye(8);epg_pen=A;
    epg_peg=eye(18);
    epg_epg=zeros(18);
    A=~eye(8);epg_dt7=[A;A;A(1:2,:)];

    A=eye(8);dt7_pen=[A circshift(A,2)];
    A=eye(8);dt7_peg=[A A A(:,1:2)];
    dt7_epg=zeros(8,18);
    dt7_dt7=~eye(8);

    M= [    pen_pen     pen_peg 5*pen_epg     pen_dt7; %weights * connectivity matrix
            peg_pen     peg_peg 5*peg_epg     peg_dt7;
         30*epg_pen  30*epg_peg   epg_epg  30*epg_dt7;
        -30*dt7_pen -30*dt7_peg   dt7_epg -10*dt7_dt7];

    pen=1:16;peg=17:34;epg=35:52;dt7=53:60;
    c=repelem(color,[16 18 18 8]);

elseif insect==3 %locust - Pisokas
    N=56;
    v=rescale(normpdf(linspace(-3,3,7),0,1));
    pen_pen=zeros(16);
    pen_peg=zeros(16);
    pen_epg=eye(16)+diag(ones(1,7),9)+diag(ones(1,7),-9);
    pen_dt7=zeros(16,8);

    peg_pen=zeros(16);
    peg_peg=zeros(16);
    peg_epg=eye(16);
    peg_dt7=zeros(16,8);

    A=eye(16);A([8 9],[9 8])=1;epg_pen=A;
    epg_peg=A;
    epg_epg=zeros(16);
    epg_dt7=zeros(8);%epg_dt7 defined in the loop

    dt7_pen=[eye(8) eye(8)];
    dt7_peg=[eye(8) eye(8)];
    dt7_epg=zeros(8,16);
    dt7_dt7=zeros(8);%dt7_dt7 defined in the loop

    for i=1:7
        epg_dt7=epg_dt7+diag(repelem(v(i),abs(i-7)+1),i);
        epg_dt7=epg_dt7+diag(repelem(v(i),abs(i-7)+1),-i);
        dt7_dt7=dt7_dt7+diag(repelem(v(i),abs(i-7)+1),i);
        dt7_dt7=dt7_dt7+diag(repelem(v(i),abs(i-7)+1),-i);end
    epg_dt7=[epg_dt7;epg_dt7];

    M= [    pen_pen     pen_peg 20*pen_epg     pen_dt7; %weights * connectivity matrix
            peg_pen     peg_peg 20*peg_epg     peg_dt7;
         10*epg_pen  10*epg_peg    epg_epg  70*epg_dt7;
        -80*dt7_pen -80*dt7_peg    dt7_epg -20*dt7_dt7];

    pen=1:16;peg=17:32;epg=33:48;dt7=49:56;
    c=repelem(color,[16 16 16 8]);

elseif insect==4 %bee - Thomas Stone
    TL_CL1=eye(16);

    A=eye(8);CL1_TB1=[A;A];

    TB1_TB1
    A=eye(8);TB1_CP4=[A A];
    A=eye(8);A=[A(:,2:end) A(:,1)];TB1_CP1=[A rot90(A,2)];

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
    figure;heatmap(M,'XDisplayLabel',nan(N,1),'YDisplayLabel',nan(N,1),'colorlimits',[-100 100])
    c1=[102 21 148]/256;c2=[84 145 3]/256; % negative / positive colors
    cm=[linspace(c1(1),1,100) 1 linspace(1,c2(1),100);
        linspace(c1(2),1,100) 1 linspace(1,c2(2),100);
        linspace(c1(3),1,100) 1 linspace(1,c2(3),100)]';
    colormap(cm);grid off
    title("Connectivity matrix ("+N+" neurons)")
    xlabel('Post-synaptic neuron');ylabel('Pre-synaptic neuron')

    figure;plot(AP.t,AP.V)
    yline(Vthr,'--')
    title('Action potential template')
    xlabel('Time (ms)');ylabel('Voltage (mV)')

    figure;plot(PC.t,PC.I)
    yline(Ipsc/2,'--');xline([2 6.944],'--')
    title('Post-synaptic current template')
    xlabel('Time (ms)');ylabel('Current (nA)')
    set(findobj('type','axe'),'TickDir','out')
    set(findobj('type','fig'),'color','w');end
