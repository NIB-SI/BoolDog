function [ J ] = get_J( u,f,xd,OCP )
%Calculates the value of the target functional
%J=0.5*sum_i(w_i*int((x(i)-xd(i))^2))+alpha*sum_j(int(u(j))), i index of nodes of interest, weights w_i for a given
%u, details see main Variables
%right hand-side f of the ordinary differential equation
%dx(t)/dt=f(x(t),u(t)) and desired state xd, details see main Variables

dt=OCP.timeInterval;
x=forward(f,u,OCP);                                     %Calculate state x
numNoi=OCP.DataNoi(:,1);                                %Indices of the nodes of interest
wNoi=OCP.DataNoi(:,3);                                  %Weights of the nodes of interest
J=0.5*dt*sum(wNoi'.*sum(transpose((x(numNoi,:)-xd).^2)))+dt*OCP.alpha*sum(sum(u));  %Calculate target functional 

end

