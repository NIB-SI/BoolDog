function [ p ] = backward_SQH( df_x,cmx,x,u,xd,OCP )
%Calculates the adjoint p according to the Pontryagin maximum principle framework,
%see corresponding paper for details
%numNodes x time steps matrix p, row i corresponds to the node i, column j
%corresponds to the time step j, entry (i,j) corresponds to the value of p
%for the node i at time step j.
%Input: numNodes x numNodes matrix df_x of function
%handles, the derivatives of f with respect to x, see createJacobian for
%details, as well as for cmx denoting the nonzero elements of df_x, external stimuli u and desired state xd, see main Variables for
%details.

N=OCP.numNodes;
dt=OCP.timeInterval;
Nt=round(OCP.timeHorizon/dt);                                                   %Number of time steps          
numNoi=OCP.DataNoi(:,1);                                                        %Index of the nodes of interest
wNoi=OCP.DataNoi(:,3);                                                          %Weights of the nodes of interest in the target functional
cnumNoi=setdiff(1:N,numNoi);                                                    %Other indices of the nodes which are not considerd in the target functional
p(:,Nt+1)=zeros(N,1);                                                           
   
nonzeroel=find(cmx(1,:));                                                       %Find the nonzero elements of the Jacobian
for j=Nt:-1:1                                                                   %Loop over all time steps backwards starting from Nt-1 because of the final time condition p(T)=x-xd
    eval_deriv=cmx;                                                             %Set the derivatives which are equal zero to zero and the other equal 1
    eval_deriv(nonzeroel)=cellfun(@(c) c(x(:,j+1)',u(:,j)'),df_x(1,nonzeroel)); %Substitute the ones with the value of the corresponding non-zero derivative, evaluate the Jacobi matrix at (x(:,j),u(:,j)), the values of x and u at time j, cellfun(@(c) c(x,u),df_x) applies the argument (x,u) to each function handle df_x{i,j}
    eval_deriv=sparse(reshape(eval_deriv,[N,N])');                              %Reshape vector to a matrix and transfrom it into a sparse one
    p_int=eval_deriv(:,:)'*p(:,j+1);
    p(cnumNoi,j)=p(cnumNoi,j+1)+dt*p_int(cnumNoi);                              %Solve adjoint equation backward in time with an explicit Euler scheme for the nodes which are not considered in the target functional
    p(numNoi,j)=p(numNoi,j+1)+dt*(wNoi.*(x(numNoi,j+1)-xd(:,j+1))+p_int(numNoi)); %Solve adjoint equation backward in time with an explicit Euler scheme for the nodes of interest  
end

end

