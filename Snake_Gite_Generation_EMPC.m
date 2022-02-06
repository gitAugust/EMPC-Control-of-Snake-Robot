clc

% % Please uncomment followingh comment without install casadi package
% currentpath = mfilename('fullpath');
% i=strfind(currentpath,'\');
% currentpath=currentpath(1:i(end));
% casadipath = fullfile(currentpath,'casadi-windows-matlabR2016a-v3.5.5');
% addpath(casadipath)

import casadi.*
%%
%Robot parameters 
R.Nl=9;
l=0.14;
R.m=1;
R.ct=1;
R.cn=3;
R.cp=(R.cn-R.ct)/(2*l);
%Control parameters
R.T = 0.05; % Sampling time 
N = 20;% MPC control horizon lenght
R.lambda1=0.5;
R.lambda2=20;
gamma=0.025;
u_max=0.2276;
phi_max=0.052;
V_phi_max=0.109;

%Modeling
%State of robot
phi = SX.sym('Phi',R.Nl-1,1);theta = SX.sym('theta' );Px = SX.sym('Px');Py = SX.sym('Py');
Vphi = SX.sym('Vphi',R.Nl-1,1);Vtheta = SX.sym('Vtheta');Vt = SX.sym('Vt');Vn = SX.sym('Vn');
states = [phi;theta;Px;Py;Vphi;Vtheta;Vt;Vn]; %State vactor 
n_states = size(states,1);%Get the number of state£¬save with formal£¨n_states, 1£©

%Control vactor
u = SX.sym('u',R.Nl-1,1);
n_controls = size(u,1);

%Construction of MPC
U = SX.sym('U',n_controls,N); % Control sequency 
P = SX.sym('P',n_states);%
X = SX.sym('X',n_states,N+1);%System state

%Kinematics Auxiliary Matrix
R.A=zeros(R.Nl-1,R.Nl);
R.D=zeros(R.Nl-1,R.Nl);
for i=1:R.Nl-1
    R.A(i,i)=1;R.A(i,i+1)=1;
    R.D(i,i)=1;R.D(i,i+1)=-1;
end
R.D_D=R.D'/(R.D*R.D');
R.e_e = ones(R.Nl-1,1);

%Set the initial stste
X(:,1) = P;
%Predict the next state
for k = 1:N
    %Predict with the snake model
    st = X(:,k);  con = U(:,k);
    st_next = (st(1:R.Nl-1)+R.T*st(R.Nl+3:2*R.Nl+1));
    st_next = [st_next;st(R.Nl)+R.T*st(2*R.Nl+2)];
    st_next = [st_next;st(R.Nl+1)+R.T*(st(2*R.Nl+3)*cos(st(R.Nl))-st(2*R.Nl+4)*sin(st(R.Nl)))];
    st_next = [st_next;st(R.Nl+2)+R.T*(st(2*R.Nl+3)*sin(st(R.Nl))+st(2*R.Nl+4)*cos(st(R.Nl)))];
    %   st_next = [st_next;st(R.Nl+3:2*R.Nl+1)+R.T*(R.cp/R.m*st(2*R.Nl+3)*R.A*R.D'*st(1:R.Nl-1)+1/R.m*R.D*R.D'*con)];
    st_next = [st_next;st(R.Nl+3:2*R.Nl+1)+R.T*con];
    st_next = [st_next;st(2*R.Nl+2)+R.T*(-R.lambda1*st(2*R.Nl+2)+R.lambda2/(R.Nl-1)*st(2*R.Nl+3)*R.e_e'*st(1:R.Nl-1))];
    st_next = [st_next;st(2*R.Nl+3)+R.T*(-R.ct/R.m*st(2*R.Nl+3)+2*R.cp/(R.Nl*R.m)*st(2*R.Nl+4)*R.e_e'*st(1:R.Nl-1)...
        -R.cp/(R.Nl*R.m)*st(1:R.Nl-1)'*R.A*R.D_D*st(R.Nl+3:2*R.Nl+1))];
    st_next = [st_next;st(2*R.Nl+4)+R.T*(-R.cn/R.m*st(2*R.Nl+4)+2*R.cp/(R.Nl*R.m)*st(2*R.Nl+3)*R.e_e'*st(1:R.Nl-1))];
    X(:,k+1) = st_next;
end

%Get the function ff about input(control sequency and initial state) and output(system state)
ff=Function('ff',{U,P},{X});

%%
obj = 0;
g = [];

for k=1:N
    %Obtain the optimization target expression in N steps
    st = X(:,k);  con = U(:,k);
    obj = obj-st(2*R.Nl+3);
end

for k = 1:N+1
    %Set Restrictions of phi & Vphi for every step
    g = [g ; X(1:R.Nl-1,k)];%state phi
    g = [g ; X(R.Nl+3:2*R.Nl+1,k)];%state Vphi
end

OPT_variables = reshape(U,n_controls*N,1);
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

%set of ipot reference https://coin-or.github.io/Ipopt/OPTIONS.html
opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;%Acceptable Convergence Tolerance
%opts.ipopt.acceptable_obj_change_tol = 1e-6;

%Get the nlp solver
solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

%%
%Set Restrictions
args = struct;
args.lbg=[];
args.ubg=[];
for k = 1:N+1
    args.lbg = [args.lbg ; ones(R.Nl-1,1)* -phi_max];
    args.lbg = [args.lbg ; ones(R.Nl-1,1)*  -V_phi_max];
    args.ubg = [args.ubg ; ones(R.Nl-1,1)* phi_max];
    args.ubg = [args.ubg ; ones(R.Nl-1,1)*  V_phi_max];
end

args.lbx(1:1:n_controls*N,1) = -u_max;
args.ubx(1:1:n_controls*N,1) = u_max;

%%
%Initialization
t0 = 0;%Simulation timer
x0 = [0,0.01,-0.01,0.01,0,0,0.01,-0.01,zeros(1,R.Nl+5,'double')]'; %initial stste
xx(:,1) = x0;
xxlu(:,1) = x0;
u0 = zeros(n_controls,N);
sim_tim = 15; %Simulation time 
mpciter = 0;%iteration counter
xx1 = [];%Store the robot position for each step
u_cl=[];%Stores all calculated control instructions

%%
%Start MPC
main_loop = tic;
while(mpciter < sim_tim / R.T)
    args.p = x0;
    args.x0 = reshape(u0,n_controls*N,1);
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    u = reshape(full(sol.x)',n_controls,N);%Get the optimal control sequency
    
    %Stor the result data
    u_cl= [u_cl ; u(:,1)'];
    [t0, x0, u0] = shiftLYC(t0, x0, u, R); %Get the next state
    
    xx(:,mpciter+2) = x0;
    %mpciter
    mpciter = mpciter + 1;
end
MPC0_loop_time = toc(main_loop)

%%
% %---------------For comparation Set another MPC problem------------------
obj = 0;

for k=1:N
    st = X(:,k);  con = U(:,k);
    obj = obj-st(2*R.Nl+3)+gamma*(con'*con);
end


OPT_variables = reshape(U,n_controls*N,1);
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);
solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

%%
t0 = 0;
x0 = [0,0.01,-0.01,0.01,0,0,0.01,-0.01,zeros(1,R.Nl+5,'double')]';
xxgamma(:,1) = x0;
u0 = zeros(n_controls,N);
mpciter = 0;
xx1 = [];

%%
main_loop = tic;
while(mpciter < sim_tim / R.T)
    args.p = x0;
    args.x0 = reshape(u0,n_controls*N,1);
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    u = reshape(full(sol.x)',n_controls,N);
    

    [t0, x0, u0] = shiftLYC(t0, x0, u, R);
    
    xxgamma(:,mpciter+2) = x0;
    %mpciter
    mpciter = mpciter + 1;
end
MPCgamma_loop_time = toc(main_loop)

%%
%Lateral undulation Gait for comperation
%-----------------------------------LU GAIT-------------------------------------
luiter=0;

A.alpha=0.05;%15*(pi/180);
A.omega=120*(pi/180);
A.delta=40*(pi/180);

A.kphi=20;
A.kvphi=5;

main_loop = tic;
u_cllu=[];
while(luiter < sim_tim / R.T)
    A.t=R.T*luiter;
    st = xxlu(:,luiter+1);
    
    phi_ref=zeros(1,R.Nl-1);
    dphi_ref=zeros(1,R.Nl-1);
    ddphi_ref=zeros(1,R.Nl-1);
    for i=1:R.Nl-1
        phi_ref(i)=A.alpha*sin(A.omega*A.t+(i-1)*A.delta);%(8.17)
        dphi_ref(i)=A.alpha*A.omega*cos(A.omega*A.t+(i-1)*A.delta);
        ddphi_ref(i)=-A.alpha*A.omega*A.omega*sin(A.omega*A.t+(i-1)*A.delta);
    end
    u_u1=ddphi_ref'+A.kvphi*(dphi_ref'-st(R.Nl+3:2*R.Nl+1))+A.kphi*(phi_ref'-st(1:R.Nl-1));%¿ØÖÆÆ÷ P134 (8.20)
    u11=R.m*(R.D*R.D')\(u_u1+R.cn/R.m*st(R.Nl+3:2*R.Nl+1)-R.cp/R.m*st(2*R.Nl+3)*R.A*R.D'*st(1:R.Nl-1));%P133 (8.19)
    u_cllu= [u_cllu ; u11(1)];
    st_next = (st(1:R.Nl-1)+R.T*st(R.Nl+3:2*R.Nl+1));
    st_next = [st_next;st(R.Nl)+R.T*st(2*R.Nl+2)];
    st_next = [st_next;st(R.Nl+1)+R.T*(st(2*R.Nl+3)*cos(st(R.Nl))-st(2*R.Nl+4)*sin(st(R.Nl)))];
    st_next = [st_next;st(R.Nl+2)+R.T*(st(2*R.Nl+3)*sin(st(R.Nl))+st(2*R.Nl+4)*cos(st(R.Nl)))];
    %   st_next = [st_next;st(R.Nl+3:2*R.Nl+1)+R.T*(R.cp/R.m*st(2*R.Nl+3)*R.A*R.D'*st(1:R.Nl-1)+1/R.m*R.D*R.D'*u1)];
    st_next = [st_next;st(R.Nl+3:2*R.Nl+1)+R.T*u_u1];
    st_next = [st_next;st(2*R.Nl+2)+R.T*(-R.lambda1*st(2*R.Nl+2)+R.lambda2/(R.Nl-1)*st(2*R.Nl+3)*R.e_e'*st(1:R.Nl-1))];
    st_next = [st_next;st(2*R.Nl+3)+R.T*(-R.ct/R.m*st(2*R.Nl+3)+2*R.cp/(R.Nl*R.m)*st(2*R.Nl+4)*R.e_e'*st(1:R.Nl-1)...
        -R.cp/(R.Nl*R.m)*st(1:R.Nl-1)'*R.A*R.D_D*st(R.Nl+3:2*R.Nl+1))];
    st_next = [st_next;st(2*R.Nl+4)+R.T*(-R.cn/R.m*st(2*R.Nl+4)+2*R.cp/(R.Nl*R.m)*st(2*R.Nl+3)*R.e_e'*st(1:R.Nl-1))];
    xxlu(:,luiter+2) = st_next;
    %luiter
    luiter = luiter + 1;
end
LU_loop_time = toc(main_loop)


%%
%--------------------------------------Plot------------------------------------
time_step=(1:1:sim_tim /R.T+1);
figure(1)
hold on
a=plot(time_step,xxlu(14,:),'color',[0.3,0.6,1],'linewidth',1.5);
b=plot(time_step,xx(14,:),'r','linewidth',1.5);
x=get(gca,'xlim');
plot(x,[V_phi_max V_phi_max],'--','color',[1,0.5,0],'linewidth',1.5);
plot(x,[-V_phi_max -V_phi_max],'--','color',[1,0.5,0],'linewidth',1.5);
legend([a,b],'LU','EMPC')
axis([0 time_step(end) -0.15 0.15])
xlabel('Time Steps (1/20 s)')
ylabel('Joint V elocity v¦Õ,3(m/s)')
box on
axes('Position',[0.65,0.2,0.2,0.3]);
axis([155 180 0 0.12])
hold on
a=plot(time_step(155:180),xxlu(14,155:180),'color',[0.3,0.6,1],'linewidth',1.5);
b=plot(time_step(155:180),xx(14,155:180),'r','linewidth',1.5);
box on
%%
figure(2)
hold on
plot(time_step,xx(21,:),'r','linewidth',1.5);
plot(time_step,xxgamma(21,:),'color',[1,0.5,0],'linewidth',1.5);
plot(time_step,xxlu(21,:),'color',[0.3,0.6,1],'linewidth',1.5);
ax = gca;
ax.YAxis.Exponent = -2;
box on
% ax.BoxStyle='full';
axis([0 time_step(end) 0 0.07])
xlabel('Time Steps (1/20 s)')
ylabel('Forward V elocity vt(m/s)')
axes('Position',[0.55,0.2,0.3,0.3]);
hold on
plot(time_step(160:210),xx(21,160:210),'r','linewidth',1.5);
plot(time_step(160:210),xxgamma(21,160:210),'color',[1,0.5,0],'linewidth',1.5);
plot(time_step(160:210),xxlu(21,160:210),'color',[0.3,0.6,1],'linewidth',1.5);
ax = gca;
ax.YAxis.Exponent = -2;
box on
