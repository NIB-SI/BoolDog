function [ xd ] = get_xd(OCP)
%get_xd creates a matrix in whose rows the desired activity levels of the
%node of interest is given as a constant function over time, see main for details.

Nt=round(OCP.timeHorizon/OCP.timeInterval);         %Number of time steps
numNoi=OCP.DataNoi(:,1);                            %Index of nodes of interest
desVal=OCP.DataNoi(:,2);                            %Desired values of the activity levels of the nodes of interest

for i=1:max(size(numNoi))                           %loop over all nodes of interest
    xd(i,:)=desVal(i)*ones(1,Nt+1);                 
end
end

