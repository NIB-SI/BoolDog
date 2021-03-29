function [u,valueTargetFunctional,count ] = projected_gradient_method( grad, proj, J,f,df_x,cmx,df_u,cmu, u, xd, tol3, max_iter, OCP )
%Output:u the external stimuli optimizing the target functional, 
        %valueTargetFunctional a 1 x count+1 row vector whose i-th column
        %corresponds to the value of the target functional J after the (i-1)-th
        %iteration of the projected gradient method, 
        %count: number of iterations of the projected gradient method
%Input:grad function handle for the gradient of the reduced target funcitonal J(u):=J(x(u),u), proj function handel of the projection, projects u into [0,1] 
        %J function handle for the target functional, f right hand-side of the ordinary differential equation corresponding to the network with dx/dt=f(x(t),u(t))
        %df_x function handle for the derivative of f with respect to x, cmx denotes the nonzero elements of df_x, df_u function handle for the derivative of f with respect to u, cmu denotes the nonzero elements of df_u, see output function createJacobian
        %u inital guess for the external stimuli, can be taken from the combinatorial method, details like format see main Variables 
        %xd desired state for the values x of the corresponding nodes, tol3 stopping criterion, max_iter maximum iteration number of the projected gradient method 

fprintf('\n');
fprintf('Starting projected gradient method...\n');
sigma=0.01;                                         %Parameter for the Armijo step size strategy
beta1=0.1;                                          %Parameter for the Armijo step size strategy
beta2=1.1;                                          %Paramter for enlarging the step size
dt=OCP.timeInterval;
g=grad(u,f,df_x,cmx,df_u,cmu,xd,OCP);               %get the gradient dJ(u)/du of the reduced target functional J(u):=J(x(u),u)
count=0;                                            %Counter for the number of iterations of the projected gradient method
valueTargetFunctional(1,1)=J(u,f,xd,OCP);           %value of the target functional with the inital guess of the external stimuli u
res=sum(sum((u-proj(u-g)).^2))*dt;                  %Residuum of the projected gradient
fprintf('Initial value of J=%d and residuum res=%d\n',valueTargetFunctional(1,1),res);
t=1;                                                %Initial step size for the step size strategy
while(res>tol3 && count<max_iter)                   %Perform the projected gradient method as long as the residuum of the projected gradient is greater than tol2 or the numer of iterations is less than max_iter
    while(J(proj(u-t*g),f,xd,OCP)-valueTargetFunctional(1,count+1)>sigma*dt*reshape(g',[1,OCP.numControls*OCP.timeHorizon/OCP.timeInterval])*reshape((proj(u-t*g)-u)',[OCP.numControls*OCP.timeHorizon/OCP.timeInterval,1])) %Check if Armijo condition is fulfilled
        t=t*beta1;                                  %If Armijo condition is not fulfilled, decrease step size t
    end
    u=proj(u-t*g);                                  %Update the external stimuli u with a projected gradient step with step size t in direction g
    g=grad(u,f,df_x,cmx,df_u,cmu,xd,OCP);           %Calculate new gradient with the new u
    count=count+1;                                  %Increase counter of the iterations
    valueTargetFunctional(1,count+1)=J(u,f,xd,OCP); %Calculate value of the target functional with the new u
    res=sum(sum((u-proj(u-g)).^2))*dt;              %Calculate residuum of the new projected gradient 
    fprintf('Iteration %i, value target functional J=%d with resdiuum res=%d and step size t=%d\n',count,valueTargetFunctional(1,count+1),res,t);
    t=t*beta2;                                      %Increase the step size
end

end

