
%% Assignment 1: 1D steady conduction heat transfer case in cylindrical walls
clc 
clear all 
close all 
%% Input Data
%% Physical Data
Rint=5; %Internal radius (m) 
Rext=10; %External Radius (m)
H=1; %Height (m)
TA=100+273.15; %temperature of the left environment (K) 
TB=25+273.15; %temperature on the right environment (K)  
alpha_A= 2000; % Boiling water insdie the cylinder (W/m2K)
alpha_B= 10; % Air with free convection outside the cylinder (W/m2K)
qv=10; % internal heat generation (W/m3)
lemda=51.9; %Thermal conductivity of Carbon steel(W/mK)
e_t=Rext-Rint; %thickness (m) 
%% Numerical Data 
maxiter=1e6; %Limit for number of iterations 
tolerance=1e-9; %Tolerance for Iterative process
N_e= 200; %Number of CVs 
d=e_t/N_e; %Thickness of each elements (m)
%% Previous Calculations and vector definitions 
%% Mesh Definition
loop_interval=1;
N_list=1:loop_interval:N_e; %Changing the number of control volumes
Error_change=zeros(1,length(N_list)); %Array for storing chagning error wrt CVs
cas=1;
for N_e=N_list %Loop for changing CVs
R_cv=linspace (Rint,Rext,N_e+1); %Positions of Xcv(m)
V=zeros(1,N_e+2); %Volume of control volume for each node
S=2*pi*R_cv*H; %Surface Area
R_nodes=[]; %Position of the nodes

%Positions of Xp and Volume calculations
for i=1:1:N_e+2
    if i==1
        R_nodes(i)=Rint;
        V(i)=0;
    elseif i<N_e+2
        R_nodes(i)=(R_cv(i-1)+R_cv(i))/2;
        V(i)=pi*(R_cv(i)^2-(R_cv(i-1)^2))*H;
    else
        R_nodes(i)=Rext;
        V(i)=0; 
    end 
end 
%% Discretization coefficients 
ap=ones(1,N_e+2);
bp=ones(1,N_e+2);
aw=ones(1,N_e+2);
ae=ones(1,N_e+2);
%Node i=1
ap(1)=(lemda/(R_nodes(2)-R_nodes(1)))+alpha_A; %first node on left edge
aw(1)=0; %first node on left edge
ae(1)=(lemda/(R_nodes(2)-R_nodes(1))); %first node on left edge
bp(1)=alpha_A*TA; %first node on left edge
%Node i=2:N+1 (Central Nodes)
for i=2:1:N_e+1    
    aw(i)=(lemda*S(i-1))/(R_nodes(i)-R_nodes(i-1)); %central nodes 
    ae(i)=(lemda*S(i))/(R_nodes(i+1)-R_nodes(i)); %central nodes 
    ap(i)=aw(i)+ae(i); %central nodes
    bp(i)=qv*V(i); %central nodes 
end 
%Node i=N+2
ap(N_e+2)=(lemda/(R_nodes(N_e+2)-R_nodes(N_e+1)))+alpha_B;%last node on right edge
aw(N_e+2)=(lemda/(R_nodes(N_e+2)-R_nodes(N_e+1)));%last node on right edge
ae(N_e+2)=0;%last node on right edge
bp(N_e+2)=alpha_B*TB;%last node on right edge
%% Gauss-Seidal Solver 
%% Initial Guess and Iterative Temperature vectors 
T_i=zeros(1,N_e+2); %Numerical Temperature array
T_prev=ones(1,N_e+2); %Guess Temperature array
iter=1; 
while iter<maxiter
    for i=1:1:N_e+2
        if i==1 %For left-most node 
            T_i(i)=(ae(i)*T_i(i+1)+bp(i))/ap(i);
        elseif i<N_e+2 %For central nodes 
            T_i(i)=(aw(i)*T_i(i-1)+ae(i)*T_i(i+1)+bp(i))/ap(i);
        else %For right-most node 
            T_i(i)=(aw(i)*T_i(i-1)+bp(i))/ap(i); 
        end
    end 
    
    residual=max(abs(T_i-T_prev));
    if residual<tolerance
        break; 
    else 
        T_prev=T_i;  
        iter=iter+1;
    end     
end 
%% Global Energy Balance 
Energy_Balance=S(1)*alpha_A*(TA-T_i(1))-S(end)*alpha_B*(T_i(end)-TB)+qv*sum(V);
%% Analytical Solution Procedure 
syms Tw1 Tw2 C1 C2
%% Applying wall and specific boundary conditions 
equa01=C1==((Tw1-Tw2)-(qv/(4*lemda))*(Rext^2-Rint^2))/(log(Rint/Rext)); 
equa02=C2==Tw1+((qv*(Rint^2))/(4*lemda))-C1*log(Rint);
equa03=(alpha_A*(TA-Tw1))-lemda*(((qv*Rint)/(2*lemda))-(C1/Rint));
equa04=(alpha_B*(Tw2-TB))-lemda*(((qv*Rext)/(2*lemda))-(C1/Rext));
%% Determine wall temperatures 
Sol_Ana=solve([equa01,equa02,equa03,equa04],[Tw1 Tw2 C1 C2]);
T_w1=double(Sol_Ana.Tw1); %Left boundary temperature
T_w2=double(Sol_Ana.Tw2); %Right boundary temperature 
C1=double(Sol_Ana.C1); 
C2=double(Sol_Ana.C2);
% fprintf('The temperature of the left wall: Tw1= %f\nThe temperature of the right wall:Tw2=%f\n',T_w1, T_w2);
%% Solving 
T=[]; %Values of temperature 
Q=[]; %Values of heat flux 
R_ana=R_nodes;
i=1;
for r=1:1:N_e+2
    T(i)=T_w2+(qv/(4*lemda))*(Rext^2-R_nodes(i)^2)-((T_w2-T_w1)+(qv/(4*lemda))*(Rext^2-Rint^2))*((log(Rext/R_nodes(i))/(log(Rext./Rint))));
    Q(i)=((qv*R_nodes(i))/2)-lemda*(((T_w1-T_w2)+(qv/(4*lemda))*(Rext^2-Rint^2))/(R_nodes(i)*log(Rint/Rext)));
    i=i+1;
end 

Error_change(N_e)=max(abs(T-T_i));

cas=cas+1;
end
%% Error between Analytical and Numerical
Final_Error=(abs((T-T_i)));
Max_Error=max(Final_Error);
%% plotting 
figure 
plot (R_ana, T,'-o'); %Plotting the analytical temperature distribution 
title ('Temperature Distribution Analytical Solution')
xlabel ('Wall distance (m)')
ylabel ('Temperature (K)')
hold on 
%Plotting the Numerical temperature distribution 
plot (R_nodes, T_i,'.-') 
title ('Temperature Distribution Numerical Solution')
xlabel ('Wall distance (m)')
ylabel ('Temperature (K)')
legend ('Analytical','Numerical') 
hold off
figure %Plotting the error between the Analytical and Numerical 
plot (R_ana, Final_Error,'-o'); 
title ('Error for 200 nodes wrt Analytical Solution')
xlabel ('Wall distance (m)')
ylabel ('Error')
figure %Evolution of Error with increasing number of control volumes 
plot (N_list, Error_change,'-o'); 
title ('Error wrt number of nodes')
xlabel ('Number of nodes')
ylabel ('Error')