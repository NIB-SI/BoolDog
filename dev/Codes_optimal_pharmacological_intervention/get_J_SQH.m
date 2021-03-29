function [ J ] = get_J_SQH( u,x,xd,OCP )
%Calculates the value of the target functional
%J=0.5*sum_i(w_i*int((x(i)-xd(i))^2))+alpha*sum_j(int(u(j))), i indices of the nodes of interest for a given
%u,x, weights w_i and desired state xd, details see main Variables

dt=OCP.timeInterval;
numNoi=OCP.DataNoi(:,1);                %Indices of the nodes of interest
wNoi=OCP.DataNoi(:,3);                  %Weights of the nodes of interest
J=0.5*dt*sum(wNoi'.*sum(transpose((x(numNoi,:)-xd).^2)))+dt*OCP.alpha*sum(sum(u));  %Calculate target functional

end

