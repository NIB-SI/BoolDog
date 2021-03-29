function [ u ] = SQH_method( J,f,df_x,cmx,df_u,cmu,tol2,u,xd,max_iter,OCP)
%Output: External stimuli u as a matrix of numControls rows ans Nt columns,
%see main.m for details
%Input: function handle J for the target functional, f right hand-side of
%the ordinary differential equation corresponding to the network with dx/dt=f(x(t),u(t)),
%df_x function handle for the derivative of f with respect to x, cmx denotes the nonzero elements of df_x, df_u function handle for the derivative of f with respect to u, cmu denotes the nonzero elements of df_u, see output function createJacobian
%u inital guess for the external stimuli, can be taken from the combinatorial method, details like format see main.m Variables 
%xd desired state for the values x of the corresponding nodes, tol2
%stopping criterion, max_iter maximum number of updates on the control u

fprintf('\n');
fprintf('Starting SQH method...\n');

dt=OCP.timeInterval;
Nt=round(OCP.timeHorizon/dt);                           %Number of time steps
N=OCP.numNodes;
M=OCP.numControls;
x=forward(f,u,OCP);                                     %Calculates the state corresponding to the initial guess u
p=backward_SQH(df_x,cmx,x,u,xd,OCP );                   %Calculates the corresponding adjoint equation
u_old=u;                                                %Copies u

epsilon=0.1;                                            %Initial guess for epsilon                         
sigma=50;                                               %Paramter for enlarging epsilon
zeta=0.15;                                              %Paramter for diminishing epsilon
eta=10^-12;                                             %Paramter for comparing the values of the target functional for new and old u
valueTargetFunctional(1)=J(u,x,xd,OCP);                 %Calulates value of target functional for the initial guess of u
counter=1;                                              %Counter counts the updates on u where counter-1 corresponds to the number of updates

nonzeroel=find(cmu(1,:));                               %Finds non zero derivatives of f with respect to the external stimuli u(j), j=1,...,m
while(1==1)
    for j=1:Nt                                                                                      %For loop over time steps
        eval_deriv=cmu;                                                                             %Sets the derivatives which are equal zero to zero and the others equal 1
        eval_deriv(nonzeroel)=sparse(cellfun(@(c) c(x(:,j+1)',u(:,j)'),df_u(nonzeroel)));           %Evaluates the Jacobi matrix at (x(:,j),u(:,j)), the values of x and u at time j, cellfun(@(c) c(x,u),df_u) applies the argument (x,u) to each function handle df_u{i,j}
        eval_deriv=sparse(reshape(eval_deriv,[M,N])');                                              %Reshapes vector to a matrix and transfrom it into a sparse one
        u(:,j)=max(0,min(((-(OCP.alpha+eval_deriv(:,:)'*p(:,j+1))/(2*epsilon))+u_old(:,j)),1));     %Update on u with the minimum of the corresponding augmented Hamiltonian, only correct if df_u is constant with respect to u
    end
    du=sum(sum((u-u_old).^2))*dt;                        %Calculating ||u-u_old||
    x_int=forward(f,u,OCP);                              %Calculating the corresponding state x to u
    J_int=J(u,x_int,xd,OCP);
    if (J_int-valueTargetFunctional(counter)>-eta*du)    %Checks if the update on u with corresponding x induces at least -eta*du descent for the target functional
        epsilon=epsilon*sigma;                           %If not enough descent enlarge epsilon
    else                                                 %If enough descent
        counter=counter+1;                          
        epsilon=epsilon*zeta;                            %Dimish epsilon
        u_old=u;                                         %Adopt u
        x=x_int;                                         %Adopt x
        p=backward_SQH(df_x,cmx,x,u,xd,OCP );            %Calculate adjoint
        valueTargetFunctional(counter)=J_int;            %Adopt value of target functional
        fprintf('Update number %i, value target functional J=%d with du=%d and epsilon=%d\n',counter-1,valueTargetFunctional(counter),du,epsilon);
    end
    if(du<tol2 || counter>max_iter)                      %If du is smaller than tol2, stop returning the last u with resulted in an descent
        fprintf('SQH method converged with du=%d\n',du);
        u=u_old;
        break;
    end
end

end

