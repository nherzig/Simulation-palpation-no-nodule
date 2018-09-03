function [Tout, Xout,Fout] = run_simutissuedyn_damp_prob(ind,mp,kp,k0,k1,cp,c1,Ts,t_ini,t_end,x0)
syms t
% Time-varying parts:
u = [ind(t);diff(ind(t),t);diff(diff(ind(t),t),t)];
u=matlabFunction(u);
% Set up simulation
dx = @(t,x)ode_sys1(x,u(t),mp,kp,k0,k1,cp,c1);
t_sim = t_ini:Ts:t_end;
[ Tout, Xout ] = ode45(dx,t_sim,x0);
uout=u(Tout');
Fout=kp*(ind(t_sim')-Xout(:,1))+cp*(uout(2,:)'+Xout(:,2));
end

% Define system equations
function dx = ode_sys1(x,u,mp,kp,k0,k1,cp,c1)

% Defined time-varying system matrices elementwise
A1=[0 1 0;0 0 1;-k1*(kp+k0)/(mp*c1) -(c1*(kp+k1+k0)+cp*k1)/(c1*mp) -(cp*c1+mp*k1)/(c1*mp)];
B1=[0 0 0; 0 0 0; kp*k1/(mp*c1) (kp*c1+k1*cp)/(c1*mp) cp/mp];

dx = A1*x + B1*u;

end
