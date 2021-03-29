function [  ] = main_comp_therapies()
%Main-file for calculating external stimuli for a regulatory network
%causing the lowest target functional value
%For details see corresponding paper

%Paramters determinded by the network
%numNodes:     Number of nodes of the network, has to correspond to the number of entries of the right hand-side f 
%numControls:  Number of external stimuli of the network, has to correspond to the number of different external stimuli used in the right hand-side f, corresponds to the maximum index of the external stimuli

%Paramters determined by the user
%timeInterval:  Number, time discretization for the ordinary differential equations
%timeHorizon:   Number, time duration, [0, timeHorizon], in which network is simulated
%alpha:         Number, Weights the contribution of the controls to the target functional
%initialState:  1 x numNodes row vector, The steady state the network is in at the beginning, that is time t=0

%Tolerance parameters
tol2=10^-14;            %Tolerance for the stopping criterion for the squential quadratic Hamiltonian method, input in function SQH_method
tol3=10^-4;             %Tolerance for the stopping criterion for the projected gradient method, input in function projected_gradient_method
max_Num=1;              %Intput in function combinatorial_method, Maximum number of external stimuli applied at once to the network in order to find a set of external stimuli that cause the desired switch, max_Num = 1,...,numControls
max_iter=10000;         %Maximum number of updates on the control for the sequential quadratic Hamiltonian method or maximum number of iterations of the projected gradient method

numColu=2;              %Time curves of the active external stimuli are plottet in a window with two columns
numColx=2;              %Time curves of the nodes of interest are plottet in a window with two columns
intv=3;                 %Number of intervals of equal length into which the range of the external stimuli [0,1] is devided

%flags
combi_method=1;                 %If flag equals 1, a combinatorical search with the function combinatorial_method by trial and error search is performed in ordert to determine external stimuli causing the desired switch, if flag equals 0 it is not performed
local_optimization_method=1;    %If flag equals 0, no local optimization method is performed, if local_optimization_method equals 1, then the sequential quadratic Hamiltonian (SQH) method is performed, if local_optimization_method equals 2, then a projected gradient method is performed


%Variables
%x:         numNodes x ((timeHorizon/timeInterval)+1)-matrix, state of the netword, row i corresponds 
%           to the  node i, i=1,...,numNodes, column j corresponds to
%           the time t=(j-1)*timeInterval, j=1,...,((timeHorizon/timeInterval)+1), entry (i,j) corresponds to value x(i) of the i-th
%           node at time t=(j-1)*timeInterval
%xd:        number nodes of interest x ((timeHorizon/timeInterval)+1)-matrix, Desired activity level of the nodes of interest, 
%           row i corresponds to the desired state of the i-th node of interest, column j corresponds to the time 
%           t=(j-1)*timeInterval, j=1,...,((timeHorizon/timeInterval)+1), entry (i,j) value of the i-th
%           desired state xd(i) of the i-th node of interest at time t=(j-1)*timeInterval
%u:         numControls x (timeHorizon/timeInterval)-matrix, external stimuli, row
%           i corresponds to the external stimulus u(i), i=1,...,numControls, column j
%           corresponds to the time t=(j-1)*timeInterval,
%           j=1,...,(timeHorizon/timeInterval), entry (i,j) value of the i-th
%           external stimulus u(i) at time t=(j-1)*timeInterval

A=[8,1,1.5;9,1,1.5;23,0,1;24,0,1;25,0,1];           %Matrix where each row is for a node of interest, first column is for its index which it has in the network, second column is its desired constant activity level, third column is its corresponding weight in the target functional

OCP=struct('numNodes',26,'numControls',1,'timeInterval',0.1,'timeHorizon',20,'alpha',0,'initialState',[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],'DataNoi',A);

xd=get_xd(OCP); %Creates the desired state in which the nodes of interest are expected to be 

%f is the right hand-side of the ordinary differential equation dx(t)/dt=f(x(t),u(t)) corresponding to the network, here as a list of function handles
f={@(x,u)-x(1),...
   @(x,u)((-exp(5)+exp(-10*(((11/10)*10*x(1)/(1+10*x(1)))-0.5)))/((1-exp(5))*(1+exp(-10*(((11/10)*10*x(1)/(1+10*x(1)))-0.5)))))-x(2)+u(1)*(1-x(2)),...
   @(x,u)((-exp(5)+exp(-10*(((11/10)*10*x(2)/(1+10*x(2)))-0.5)))/((1-exp(5))*(1+exp(-10*(((11/10)*10*x(2)/(1+10*x(2)))-0.5)))))-x(3),...
   @(x,u)((-exp(5)+exp(-10*(((13/12)*(10*x(3)+x(12)+x(19))/(1+10*x(3)+x(12)+x(19)))-0.5)))/((1-exp(5))*(1+exp(-10*(((13/12)*(10*x(3)+x(12)+x(19))/(1+10*x(3)+x(12)+x(19)))-0.5)))))-x(4),...
   @(x,u)((-exp(5)+exp(-10*((((11/10)*10*x(4)/(1+10*x(4)))*(1-(1.01/0.01)*0.01*x(14)/(1+0.01*x(14))))-0.5)))/((1-exp(5))*(1+exp(-10*((((11/10)*10*x(4)/(1+10*x(4)))*(1-(1.01/0.01)*0.01*x(14)/(1+0.01*x(14))))-0.5)))))-x(5),...
   @(x,u)((-exp(5)+exp(-10*((2*x(5)/(1+x(5)))-0.5)))/((1-exp(5))*(1+exp(-10*((2*x(5)/(1+x(5)))-0.5)))))-x(6),...
   @(x,u)((-exp(5)+exp(-10*(((1.001/0.001)*0.001*x(6)/(1+0.001*x(6)))-0.5)))/((1-exp(5))*(1+exp(-10*(((1.001/0.001)*0.001*x(6)/(1+0.001*x(6)))-0.5)))))-x(7),...
   @(x,u)((-exp(5)+exp(-10*(((1.01/0.01)*0.01*x(7)/(1+0.01*x(7)))-0.5)))/((1-exp(5))*(1+exp(-10*(((1.01/0.01)*0.01*x(7)/(1+0.01*x(7)))-0.5)))))-x(8),...
   @(x,u)((-exp(5)+exp(-10*(((1.01/0.01)*0.01*x(7)/(1+0.01*x(7)))-0.5)))/((1-exp(5))*(1+exp(-10*(((1.01/0.01)*0.01*x(7)/(1+0.01*x(7)))-0.5)))))-x(9),...
   @(x,u)-x(10),...
   @(x,u)((-exp(5)+exp(-10*((2*x(10)/(1+x(10)))-0.5)))/((1-exp(5))*(1+exp(-10*((2*x(10)/(1+x(10)))-0.5)))))-x(11),...
   @(x,u)((-exp(5)+exp(-10*(((6/5)*5*x(11)/(1+5*x(11)))-0.5)))/((1-exp(5))*(1+exp(-10*(((6/5)*5*x(11)/(1+5*x(11)))-0.5)))))-x(12),...
   @(x,u)((-exp(5)+exp(-10*(((31/30)*(10*x(3)+10*x(12)+10*x(20))/(1+10*x(3)+10*x(12)+10*x(20)))-0.5)))/((1-exp(5))*(1+exp(-10*(((31/30)*(10*x(3)+10*x(12)+10*x(20))/(1+10*x(3)+10*x(12)+10*x(20)))-0.5)))))-x(13),...
   @(x,u)((-exp(5)+exp(-10*(((2*x(26)/(1+x(26)))*(1-(11/10)*10*x(13)/(1+10*x(13))))-0.5)))/((1-exp(5))*(1+exp(-10*(((2*x(26)/(1+x(26)))*(1-(11/10)*10*x(13)/(1+10*x(13))))-0.5)))))-x(14),...
   @(x,u)((-exp(5)+exp(-10*(((11/10)*10*x(13)/(1+10*x(13)))-0.5)))/((1-exp(5))*(1+exp(-10*(((11/10)*10*x(13)/(1+10*x(13)))-0.5)))))-x(15),...
   @(x,u)((-exp(5)+exp(-10*(((2*x(26)/(1+x(26)))*(1-(31/30)*30*x(15)/(1+30*x(15))))-0.5)))/((1-exp(5))*(1+exp(-10*(((2*x(26)/(1+x(26)))*(1-(31/30)*30*x(15)/(1+30*x(15))))-0.5)))))-x(16),...
   @(x,u)((-exp(5)+exp(-10*(((2*x(26)/(1+x(26)))*(1-(101/100)*100*x(7)/(1+100*x(7))))-0.5)))/((1-exp(5))*(1+exp(-10*(((2*x(26)/(1+x(26)))*(1-(101/100)*100*x(7)/(1+100*x(7))))-0.5)))))-x(17),...
   @(x,u)((-exp(5)+exp(-10*(((11/10)*10*x(10)/(1+10*x(10)))-0.5)))/((1-exp(5))*(1+exp(-10*(((11/10)*10*x(10)/(1+10*x(10)))-0.5)))))-x(18),...
   @(x,u)((-exp(5)+exp(-10*((((11/10)*10*x(18)/(1+10*x(18)))*(1-(1.1/0.1)*0.1*x(16)/(1+0.1*x(16))))-0.5)))/((1-exp(5))*(1+exp(-10*((((11/10)*10*x(18)/(1+10*x(18)))*(1-(1.1/0.1)*0.1*x(16)/(1+0.1*x(16))))-0.5)))))-x(19),...
   @(x,u)((-exp(5)+exp(-10*(((11/10)*10*x(19)/(1+10*x(19)))-0.5)))/((1-exp(5))*(1+exp(-10*(((11/10)*10*x(19)/(1+10*x(19)))-0.5)))))-x(20),...
   @(x,u)((-exp(5)+exp(-10*(((21/20)*(10*x(12)+10*x(19))/(1+10*x(12)+10*x(19)))-0.5)))/((1-exp(5))*(1+exp(-10*(((21/20)*(10*x(12)+10*x(19))/(1+10*x(12)+10*x(19)))-0.5)))))-x(21),...
   @(x,u)((-exp(5)+exp(-10*((((101/100)*100*x(21)/(1+100*x(21)))*(1-2*x(17)/(1+x(17))))-0.5)))/((1-exp(5))*(1+exp(-10*((((101/100)*100*x(21)/(1+100*x(21)))*(1-2*x(17)/(1+x(17))))-0.5)))))-x(22),...
   @(x,u)((-exp(5)+exp(-10*(((9/8)*8*x(22)/(1+8*x(22)))-0.5)))/((1-exp(5))*(1+exp(-10*(((9/8)*8*x(22)/(1+8*x(22)))-0.5)))))-x(23),...
   @(x,u)((-exp(5)+exp(-10*(((8/7)*7*x(22)/(1+7*x(22)))-0.5)))/((1-exp(5))*(1+exp(-10*(((8/7)*7*x(22)/(1+7*x(22)))-0.5)))))-x(24),...
   @(x,u)((-exp(5)+exp(-10*(((7/6)*6*x(22)/(1+6*x(22)))-0.5)))/((1-exp(5))*(1+exp(-10*(((7/6)*6*x(22)/(1+6*x(22)))-0.5)))))-x(25),...
   @(x,u)(1-x(26))};

u=zeros(OCP.numControls,round(OCP.timeHorizon/OCP.timeInterval));   %Initial guess for the controls if no heuristical search is performed before the local optimization framework, any control can be set to any value between 0 and 1 for an initial guess for example ones() instead of zeros() 

if(combi_method==1)                                                 %Block for the combinatorial method
    u=combinatorial_method(f,xd,max_Num,intv,OCP);
end

if(local_optimization_method==1 || local_optimization_method==2)                             %Block for the local optimization method                                                         
    [df_x,cmx,df_u,cmu]=createJacobian(f,OCP);                                               %Creates derivatives of f with respect to x and u, Jacobian of the right hand side f
    if(local_optimization_method==1)
        u=SQH_method( @get_J_SQH,f,df_x,cmx,df_u,cmu,tol2,u,xd,max_iter,OCP); %Sequential quadratic Hamiltonian method as a local optimization scheme, returns u, u the external stimuli optimizing the target functional
                                                                              %Input: @get_J function handle for the target functional, f right hand-side of the ordinary differential equation corresponding to the network with dx/dt=f(x(t),u(t))
                                                                              %df_x function handle for the derivative of f with respect to x, cmx notes the nonzero elements of df_x, df_u function handle for the derivative of f with respect to u, cmu notes the nonzero elements of df_u, see output function createJacobian
                                                                              %u inital guess for the external stimuli, can be taken from the combinatorial method 
                                                                              %xd desired state for the values x of the corresponding nodes, tol2 stopping criterion, max_iter maximum number of updates on the control u of the sequential quadratic Hamiltonian
    end                                                              
                                                                                                                                                            
    if (local_optimization_method==2)
        [u,~,~]=projected_gradient_method( @get_gradient, @projection, @get_J,f,df_x,cmx,df_u,cmu, u, xd, tol3, max_iter, OCP );    %Prjected gradient method as a local optimization scheme, returns [u,J,count], u the external stimuli optimizing the target functional
                                                                                                                                    %Input:@get_gradient function handle for the gradient of the reduced target funcitonal, @projection function handel of the projection, projects u into [0,1] 
                                                                                                                                    %@get_J function handle for the target functional, f right hand-side of the ordinary differential equation corresponding to the network with dx/dt=f(x(t),u(t))
                                                                                                                                    %df_x function handle for the derivative of f with respect to x, cmx notes the nonzero elements of df_x, df_u function handle for the derivative of f with respect to u, cmu notes the nonzero elements of df_u, see output function createJacobian
                                                                                                                                    %u inital guess for the external stimuli, can be taken from the combinatorial method 
                                                                                                                                    %xd desired state for the values x of the corresponding nodes, tol2 stopping criterion, max_iter maximum iteration number of the projected gradient method  
    end
end


drawStimuli( u,numColu,OCP );           %Draws time curves of the of the active external stimuli            

x=forward(f,u,OCP);                     %Calculates the state of the network corresponding to the external stimuli u calculated with the schemes above

drawStates(x,numColx,OCP);              %Draws time curves of the activity level of the nodes of interest

fprintf('\n');
fprintf('Save data to file...\n');
dlmwrite('x.txt',[0:OCP.timeInterval:round(OCP.timeHorizon/OCP.timeInterval)*OCP.timeInterval;x]);               %Writes the state x in a text-file "x.txt" where the first row corresponds to the discrete time steps, separated by commas, row i=2,...,numNode+1 corresponds to state x(i-1), values in the colums, separated by commas, value of x(i) at the corresponding time step  
dlmwrite('u.txt',[0:OCP.timeInterval:(round(OCP.timeHorizon/OCP.timeInterval)-1)*OCP.timeInterval;u]);           %Writes the external stimuli u in a text-file "u.txt" where the first row corresponds to the discrete time steps, separated by commas, row i=2,...,numControls+1 corresponds to external stimlus u(i), values in the colums, separated by commas, value of u(i) at the corresponding time step
fprintf('Done!\n');
end

