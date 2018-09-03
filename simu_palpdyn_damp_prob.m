clear all
close all

fs=12; %fontsize
kp=[10:50:1550]; %stiffness probe
FF=zeros(size(kp))

for i=1:length(kp)
%     ind = @(t)(tanh(10*(t-0.5))+1)*-2e-3; %indentation step
    indcc=3;
    indfreq=1;
    ind= @(t)(indcc/2*cos(2*pi*indfreq*t)-indcc/2)*1e-3; %ind sinus
    mp=0.3; %masse probe

    eta1=2.0218e5;
    E0=5.8309e4;
    E1=5.3513e4;
    A=0.01*0.01;
    l0=20e-3;

    k0=A*E0/l0; 
    k1=A*E1/l0;
    cp=30;
    c1=A*eta1/l0;      

    % k0=500; 
    % k1=500;
    % kp=1000; %stiffness probe
    % cp=100;
    % c1=900;       

    Ts = 1e-3; %Time step
    t_ini=0; %Initial time
    t_end = 2; %Simulation time
    x0=[0 0 0]'; % Initial state

    [t1,x1,F1]=run_simutissuedyn_damp_prob(ind,mp,kp(i),k0,k1,cp,c1,Ts,t_ini,t_end,x0); % run simulation for safe tissue

figure(1)
%subplot(2,1,1)
hold on
plot(t1,ind(t1),'b',t1,x1(:,1),'r')
xlabel('Time [s]')
ylabel('displacement [m]')
legend('indentation','x1')

figure(2)
%subplot(2,1,2)
hold on
plot(t1,F1)
xlabel('Time [s]')
ylabel('Force [N]')

figure(3)
hold on
plot(kp(i),max(F1)-min(F1),'o')
FF(i)=max(F1)-min(F1);
end
legend(num2str(kp))