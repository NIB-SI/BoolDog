function [ grad] = get_gradient(u,f,df_x,cmx,df_u,cmu,xd,OCP )
%Calculates the gradient grad of the Lagrange function dL(x,u,p)/du with respect
%to u at (x,u,p), where x is the state and p the adjoint, 
%which corresponds to the gradient dJ(u)/du of the reduced target functional J(u):=J(x(u),u), 
%see corresponding paper for details.

%Output:
%numControls x time steps matrix grad, row i corresponds to the external stimulus i, column j
%corresponds to the time step j, entry (i,j) corresponds to value of the
%derivative dL(x,u,p)/d(u_i)^j of the Lagrange function with respect to (u_i)^j, 
%(u_i)^j variable corresponding to the value of the external stimulus i at time step j 

%Input: numControls x time steps matrix of external stimuli u, see main Variables for
        %details
        %right hand-side f, numNodes x numNodes matrix df_x of function handles, 
        %the derivatives of f with respect to x, see createJacobian for
        %details, as well as for cmx, which denotes the nonzero elements
        %of df_x
        %numNodes x numControls matrix df_u of function handles,
        %the derivatives of f with respect to u, see createJacobian for
        %details, as well as for cmu, which denotes the nonzero elements of df_u 
        %desired state xd, see main Variables for details.

dt=OCP.timeInterval;
Nt=round(OCP.timeHorizon/dt);                                                       %Number of time steps
N=OCP.numNodes;
M=OCP.numControls;
x=forward(f,u,OCP);                                                                 %Calculating state corresponding to u
p=backward(df_x,cmx,u,x,xd,OCP);                                                    %Calculating the adjoint 
x=x(:,1:Nt);
grad=OCP.alpha*ones(M,Nt);                                                          %Term of the gradient originating from the targe functional

nonzeroel=find(cmu(1,:));
for j=1:Nt                                                                          %Loop over all time steps
    eval_deriv=cmu;                                                                 %Set the derivatives which are equal zero to zero and the other equal 1
    eval_deriv(nonzeroel)=sparse(cellfun(@(c) c(x(:,j)',u(:,j)'),df_u(nonzeroel))); %Evaluate the Jacobi matrix at (x(:,j),u(:,j)), the values of x and u at time j, cellfun(@(c) c(x,u),df_u) applies the argument (x,u) to each function handle df_u{i,j}
    eval_deriv=sparse(reshape(eval_deriv,[M,N])');                                  %Reshape vector to a matrix and transfrom it into a sparse one
    grad(:,j)=grad(:,j)+dt*eval_deriv(:,:)'*p(:,j);                                 %Assembling gradient
end 

end

