clear all
clc
nbus = 5; % IEEE-5
Y = ybusppg(nbus); % Calling ybusppg.m to get Y-Bus Matrix..
busd = busdatas(nbus);
BMva=100; % Calling busdatas.. % Base MVA..
Pg = busd(:,5)/BMva; % gen erated real power
Qg = busd(:,6)/BMva; % generated reactive power.
Pl = busd(:,7)/BMva; % load real power
Ql = busd(:,8)/BMva; % load reactive power
Qlim1 = busd(:,9)/BMva;
Qlim2 = busd(:,10)/BMva;
P = Pg - Pl; % Pi = PGi - PLi..
Q = Qg - Ql; % Qi = QGi - QLi..
Psp = P ; % P Specified..
Qsp = Q ; % q specified
Qmin=Qlim1(3);
Qmax=Qlim2(3);
G = real(Y) ; % Conductance matrix..
B = imag(Y) ; % Susceptance matrix..
%---------------------PSO PARAMETERS INITIALIZATION --------------%
particle=[];
mn=[];
fr1=[];
it=input('maximum no. of iterations');
p=input('enter initial no. of particles');
pf=input('final no. of particles'); % no of particle
fs=input('frequency for sorting');
r=input('enter r');
% rfv=input('enter rfv');
% rft=input('enter rft');
c12=[];
tic
particle(1)=p;
rp=1;
T=10;
fr=[];
count=1;
deltai=zeros(p,1);
zeta=0;
i2=0;
rg=1;
rf=1;
%f=[];
f=zeros(p,it);
fp=zeros(p,1);
thp=zeros(p,5);
thg=zeros(p,5);
vp=zeros(p,5);
vg=zeros(p,5);
rft=[];
rfv=[];
for j=1:p
    rft(j)=1;
    rfv(j)=1;
end
%v=[];
%th=[];
%vv=[];
%vth=[];
v=zeros(p,it,5);
th=zeros(p,it,5);
vth=zeros(p,it,5);
vv=zeros(p,it,5);
%vtemp=zeros(p,it,5);
%thtemp=zeros(p,it,5);
%vthtemp=zeros(p,it,5);
%vvtemp=zeros(p,it,5);
%ftemp=zeros(p,it);
a=.5;
b=-0.5;
vth(:,1,:)=a+(b-a)*rand(p,5); %initial velocity of theta vector%
a=-.1;
b=0.1;
vv(:,1,:)=a+(b-a)*rand(p,5); %initial velocity of voltage vector%
vth(:,:,1)=0;
vv(:,:,1)=0;
a=0.5;
b=-0.5;
th(:,1,:)=a+(b-a)*rand(p,5);
a=1.1;
b=0.9;
v(:,1,:)=a+(b-a)*rand(p,5);
v(:,:,1)=1.02;
v(:,:,3)=1.04;
th(:,:,1)=0;
vg(:,1)=1.02;
vp(:,1)=1.02;
vg(:,3)=1.04;
vp(:,3)=1.04;
thp(:,1)=0;
thg(:,1)=0;
%-----------------------initial value of objective function-------------%
% Calculate P and Q
PVIND=zeros(p,1);
fev=0;
for j=1:p
    P = zeros(nbus,1);
    Q = zeros(nbus,1);
    MPS=zeros(p,1);
    MQS=zeros(p,1);
    for i = 2:nbus
        for k = 1:nbus
            P(i) = P(i) + v(j,1,i)* v(j,1,k)*(G(i,k)*cos(th(j,1,i)-th(j,1,k)) + B(i,k)*sin(th(j,1,i)-th(j,1,k)));
        end
    end
for i = 2:nbus
    for k = 1:nbus
         Q(i) = Q(i) + v(j,1,i)* v(j,1,k)*(G(i,k)*sin(th(j,1,i)-th(j,1,k)) - B(i,k)*cos(th(j,1,i)-th(j,1,k)));
    end
end
% real power mismatch
MP=P-Psp;
MPS=MP.^2;
%reactive power mismatch in third bus
Qsp(3)=Q(3);
if Q(3)<Qmin;
    Q(3)=Qmin;
    PVIND(j,1)=1;
else 
    PVIND(j,1)=0;
end
if Q(3)>Qmax;
    Q(3)=Qmax;
    PVIND(j,1)=1;
else 
    PVIND(j,1)=0;
end
%reactive power mismatch
MQ=Q-Qsp;
MQ(3)=0;
MQS=MQ.^2;
%objective function value
f(j,1)=sum(MPS)+sum(MQS);
% fr(j,1)=f(j,1);
% fr(j,2)=j;
fev=fev+1;
end
%Initial personal best values
for i=1:p
    for k=2:5
        thp(i,k)=th(i,1,k);
    end
    for k=2:5
        vp(i,k)=v(i,1,k);
    end
end
%for Initial Global best values updation
fmin=min(f(:,1));
for k=1:p
    if f(k,1)==fmin
       gb=k;
    else
    end
end
%Initial global best value
for k=1:p
    for j=2:5
        thg(k,j)=th(gb,1,j);
    end
    for j=2:5
        vg(k,j)=v(gb,1,j);
    end
end
fgm = min(f(:,1));
Q3=zeros(p,it);
for i=1:it
%for inertia weight W
%wmax=.4;
%wmin=.4;
% w=wmax-((wmax-wmin)*i/it);
%w =0.1+(rand()/2);
%velocity update
%position update
      w=.7;
      for j=1:p
          for k=2:5
                vth(j,(i+1),k) = w*vth(j,i,k) + rp*rand()*(thp(j,k)-th(j,i,k)) + rg*rand()*(thg(j,k)-th(j,i,k));
% if vth(j,(i+1),k)<-0.1
% vth(j,(i+1),k)=-0.1;
%end
%if vth(j,(i+1),k)>0.1
% vth(j,(i+1),k)=0.1;
%end
                th(j,(i+1),k) = th(j,i,k) + rft(j)*vth(j,(i+1),k);
          end
          for q=2:5
                vv(j,(i+1),q) = w*vv(j,i,q) + rp*rand()*(vp(j,q)-v(j,i,q)) + rg*rand()*(vg(j,q)-v(j,i,q));
                v(j,(i+1),q) = v(j,i,q) + rfv(j)*vv(j,(i+1),q);
          end
          for q=3
                if PVIND(j,1)==0
                    v(j,(i+1),q)=1.04;
                end
          end
       end
%th(:,i,5)
%objective function value
      for j=1:p
          P = zeros(nbus,1);
          Q = zeros(nbus,1);
          MPS=zeros(p,1);
          MQS=zeros(p,1);
          for m = 2:nbus
              for k = 1:nbus
                    P(m) = P(m) + v(j,(i+1),m)* v(j,(i+1),k)*(G(m,k)*cos(th(j,(i+1),m)-th(j,(i+1),k)) + B(m,k)*sin(th(j,(i+1),m)-th(j,(i+1),k)));
              end
          end
      for m = 2:5
          for k = 1:nbus
                Q(m) = Q(m) + v(j,(i+1),m)* v(j,(i+1),k)*(G(m,k)*sin(th(j,(i+1),m)-th(j,(i+1),k)) - B(m,k)*cos(th(j,(i+1),m)-th(j,(i+1),k)));
          end
      end
% real power mismatch
MP=P-Psp;
MPS=MP.^2;
%reactive power mismatch in third bus
Qsp(3)=Q(3);
      if Q(3)<Qmin;
          Q(3)=Qmin;
          PVIND(j,1)=1;
      else 
          PVIND(j,1)=0;
      end
      if Q(3)>Qmax;
         Q(3)=Qmax;
         PVIND(j,1)=1;
     else 
          PVIND(j,1)=0;
      end
      Q3(j,i)=Q(3);
      %reactive power mismatch
      MQ=Q-Qsp;
      MQ(3)=0;
      MQS=MQ.^2;
      %objective function value
      f(j,(i+1))=sum(MPS)+sum(MQS);
      fr(j,1)=f(j,i+1);
      fr(j,2)=j;
      fev=fev+1;
end
%personal best values updatio
for j=1:p
    P = zeros(nbus,1);
    Q = zeros(nbus,1);
    MPS=zeros(p,1);
    MQS=zeros(p,1);
    for t =2:nbus
          for k = 1:nbus
                P(t) = P(t) + vp(j,t)* vp(j,k)*(G(t,k)*cos(thp(j,t)-thp(j,k)) + B(t,k)*sin(thp(j,t)-thp(j,k)));
          end
    end
    for t = 2:nbus
          for k = 1:nbus
                Q(t) = Q(t) + vp(j,t)* vp(j,k)*(G(t,k)*sin(thp(j,t)-thp(j,k)) - B(t,k)*cos(thp(j,t)-thp(j,k)));
          end 
    end
      % real power mismatch
      MP=P-Psp;
      MPS=MP.^2;
      %reactive power mismatch in third bus
      Qsp(3)=Q(3);
      if Q(3)<Qmin;
           Q(3)=Qmin;
           PVIND(j,1)=1;
      else 
           PVIND(j,1)=0;
      end
      if Q(3)>Qmax;
          Q(3)=Qmax;
          PVIND(j,1)=1;
      else
          PVIND(j,1)=0;
      end
      %reactive power mismatch
      MQ=Q-Qsp;
      MQ(3)=0;
      MQS=MQ.^2;
      %objective function value
      %objective function value
      fp(j)=sum(MPS)+sum(MQS);
end
if(p>pf && mod(i,fs)==0)
  fr=sortrows(fr);
  fr1=fr;
  k1=1;
  k2=0;
  k3=((p/r));
  k4=1;
  for i1=(k3+1):p
      mn(1,k1)=fr(i1,2);
      k1=k1+1;
  end
  v(mn,:,:)=[];
  th(mn,:,:)=[];
  vv(mn,:,:)=[];
  vth(mn,:,:)=[];
  f(mn,:,:)=[];
  fr(mn,:)=[];
  p=p/r;
end
mn=[];
%personal best value updation
  for k=1:p
      for m=2:5
          if f(k,i+1)<fp(k)
              thp(k,m)=th(k,i+1,m);
          else
      end
  end
end
for k=1:p
    for m=2:5
        if f(k,i+1)<fp(k)
            vp(k,m)=v(k,i+1,m);
        else
        end
    end
end
%for Global best values updation
fgm=min(f(:,(i+1)));
for m=2:5
      for k=1:p
          if f(k,i+1)==fgm
             for l=2:p
                  thg(l,m) = th(k,i+1,m); %global best values
              end
          else
          end
      end

end
for m=2:5
        for k=1:p
            if f(k,i+1)==fgm
                for l=2:p
                    vg(l,m) = v(k,i+1,m); %global best values
                end
            else
                end
        end
end
%stopping
gb1=gb;
fgm=min(f(:,(i+1)));
for k=1:p
      if f(k,i+1)==fgm
          gb1=k;
      else
      end
end
if(p==1)
    break
end
if (abs(f(gb1,i+1)-f(gb1,i))<=10^(-T))
      if(abs(f(gb1,i+1))<=10^(-T))
            break
      end
end
if abs(max(vv(gb1,i,:)))<=10^-(2.9) && abs(max(vth(gb1,i,:)))<=10^-(2.9) && i>=10 && abs(max(vv(gb1,i,:)))>=10^-(3.1) && abs(max(vth(gb1,i,:)))>=10^-(3.1)
      for j=1:p
            rfv(j)=(((abs(max(vv(gb1,i,:))/max((vv(j,i,:)))))^-1));
            rft(j)=(((abs(max(vth(gb1,i,:))/max((vth(j,i,:)))))^-1));
      end
end
end
%Q3(:,it)
bus=zeros(5,10) ;
bus(1,3)=v(gb1,i,1);
bus(2,3)=v(gb1,i,2);
bus(3,3)=v(gb1,i,3);
bus(4,3)=v(gb1,i,4);
bus(5,3)=v(gb1,i,5);
%bus angle updation
bus(1,4)=th(gb1,i,1);
bus(2,4)=th(gb1,i,2);
bus(3,4)=th(gb1,i,3);
bus(4,4)=th(gb1,i,4);
bus(5,4)=th(gb1,i,5);
bus(3,6)=Q3(gb1,it)*BMva;
x=bus(3,6);
V = bus(:,3) ; % Specified Voltage..
del = bus(:,4) ; % Voltage Angle..
toc
%load flow function calling
loadflow(nbus,V,del,BMva,x);
disp('function value');
f(gb1,i)
disp('no. of function evaluations =');
fev
