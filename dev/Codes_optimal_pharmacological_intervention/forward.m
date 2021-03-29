function [ x] = forward(f,u,OCP)
%Calculates the corresponding state x for the ordinary differential equation dx(t)/dt=f(x(t),u(t)) 
%for a given right hand-side f, see main for details, and given external stimuli u, see main.m Variables for details about u, x 

N=OCP.numNodes;
dt=OCP.timeInterval;
T=OCP.timeHorizon;
Nt=round(T/dt);                                               %Number of time steps
x=zeros(N,Nt+1);                                              %Number of time steps +1 for the inital state
x(:,1)=OCP.initialState;                                      %Initial value of the state at t=0 equals the intial state of the network

for j=1:Nt                                                    %For loop over all time steps                                     
    x(:,j+1)=x(:,j)+dt*cellfun(@(c) c(x(:,j)',u(:,j)'),f)';   %Calculate the value of the state x at time t=t+dt with an explicit Euler method x(t+dtd)=x(t)+dt*f(x(t),u(t)), detailed x_i(t+dt)=x_i(t) + dt*f(x_1(t),...,x_numNodes(t),u_1(t),...,u_numControls(t)) for all i=1,...,numNodes                                             
end                                                           %cellfun(@(c) c(x,u),f) applies the argument (x,u) to each function handle f{i}



end

