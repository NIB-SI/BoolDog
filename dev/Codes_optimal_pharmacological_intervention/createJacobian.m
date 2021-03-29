function [ df_x,cmx,df_u,cmu ] = createJacobian(f,OCP)
%Creates the Jacobian of the function handle f(x,u) with respect to u and x
%Output: 1xnumNodes*numNodes vector of function handles df_x containing the derivatives of f with respect to x, 1xnumNodes*numNodes matrix of 0 and 1 cmx to note the non-zero elements of the derivatives of df_x,  
%1xnumNodes*numControls vector of function handles df_u containing the derivatives of f with respect to u, 1xnumNodes*numControls matrix of 0 and 1 cmu to note the non-zero elements of the derivatives of df_u,  

fprintf('\n');
fprintf('Creating Jacobian...\n');
N=OCP.numNodes;
M=OCP.numControls;
x=sym('x',[1,N]);                                           %Creating symbolic variables
u=sym('u',[1,M]);

tic;                                                        %Start time to measure time for creating Jacobian

%Derivative of f with respct to x
parfor i=1:N                                                %If no Parallel Computing Toolbox is installed, use for instead of parfor                                   
    for j=1:N
        if(max(ismember(symvar(f{i}(x,u)),x(j)))==1)        %Ceck if f{i} depends on the variable x(j), if yes calculate derivative, else set it to the zero function
            h_x=diff(f{i}(x,u),x(j));                       %Derivate f_i with respct to x_j
            df_x{i,j}=matlabFunction(h_x,'vars',{x,u});     %Create function handle from symbolic expression
            cmx(i,j)=1                                      %Note that the derivative f_i depends on x_j
        else
            df_x{i,j}=@(x,u)0;                              %As f_i does not depend on x_j and thus df_i/dx_j is constant zero, set the corresponding function to zero
            cmx(i,j)=0;                                     %Note that df_i/dx_j is the constant zero function
        end
    end
end

%Derivative of f with respect to u
parfor i=1:N                                                %If no Parallel Computing Toolbox is installed, use for instead of parfor
    for j=1:M
        if(max(ismember(symvar(f{i}(x,u)),u(j)))==1)        %Ceck if f{i} depends on the variable u(j), if yes calculate derivative, else set it to the zero function
            h_u=diff(f{i}(x,u),u(j));                       %Derivate f_i with respct to u_j
            df_u{i,j}=matlabFunction(h_u,'vars',{x,u});     %Create function handle from symbolic expression
            cmu(i,j)=1                                      %Note that the derivative f_i depends on u_j
        else
            df_u{i,j}=@(x,u)0;                              %As f_i does not depend on u_j and thus df_i/du_j is constant zero, set the corresponding function to zero
            cmu(i,j)=0;                                     %Note that df_i/du_j is the constant zero function
        end
        
    end
end

df_x=reshape(df_x',[1,N*N]);                                %Reshape matrices to vectors
cmx=reshape(cmx',[1,N*N]);
df_u=reshape(df_u',[1,M*N]);
cmu=reshape(cmu',[1,N*M]);

fprintf('Jacobian created in %d seconds\n\n', toc);           

end

