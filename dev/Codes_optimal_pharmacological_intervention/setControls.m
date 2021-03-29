function [ u ] = setControls( A,OCP)
%Input: Matrix i X 2 Matrix, i=1;...,numControls, i external Stimuli, first
%column number of the corresponding external stimulus, second column its
%duration of application.
%Output: Matrix u, details see main.m Variables

T=OCP.timeHorizon;
dt=OCP.timeInterval;
Nt=round(T/dt);                                 %Number of time steps
if(isempty(A)==1)                               %Test if A is not emty, if emty return u as a zero matrix, that means no external stimulus is active
    u=zeros(OCP.numControls,Nt);                %Nt time steps instead of Nt+1 as the control is from t=0 to t=timeHorizon-timeInterval, no control at the final time possible
else
    u=zeros(OCP.numControls,Nt);
    for i=1:OCP.numControls                     %Loop over all external stimuli
        [pos,is_in]=find(A(:,1)==i);            %Test if external stimulus i is in A, return its row in pos and is_in=1 if yes.
        if(is_in==1)                            %Construct exeternal stimulus i as a function over time
            Nt1=round((A(pos,2)/T)*Nt);         %Duration of application in time steps
            Nt0=Nt-Nt1;                         
            u(i,:)=[ones(1,Nt1),zeros(1,Nt0)];  %Set external stimulus as a function of time with full application, that means u(i,1:Nt1)=1, for time steps 1 to Nt1 and zeros, that means no application for the time steps Nt1+1 to the final time, u(i,Nt1+1:Nt)=0
        end   
    end
end

end

