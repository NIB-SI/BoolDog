function [ures] = combinatorial_method(f,xd,max_Num,intv,OCP)
%Returns a i x round(OCP.T/OCP.timeInterval) matrix A, i=1,...,max_Num, containing a combination of
%external stimuli, each with the same constant value, causing the smallest target functional value. 
%Textoutput: The index of active external stimuli and their corresponding
%target functional value

valu=1/intv:1/intv:1;                                           %Devide range of the external stimuli [0,1] into number of intervals
T=OCP.timeHorizon;
dt=OCP.timeInterval;
Nt=round(T/dt);
numNoi=OCP.DataNoi(:,1);                                        %Read out the index of the nodes of interest
wNoi=OCP.DataNoi(:,3);                                          %Read out their weight in the target functional
u0=zeros(OCP.numControls,Nt);                                   %Initial guess for the external stimuli by constant zero functions
ures=u0;                                                        %ures: current combination of external stimuli with the smallest target functional value
x=forward(f,u0,OCP);                                            %Calculating the state x corresponding to the external stimuli u
Jmin=0.5*dt*sum(wNoi'.*sum(transpose((x(numNoi,:)-xd).^2)));    %Calculates the correponding target functional value
fprintf('\n');
fprintf('The target functional value with no active stimulus is %d\n\n',Jmin);
for i=1:max_Num                                                 %Going through all possible length of combinations from length i=1 to i=max_Num
    fprintf('Combinations with %i external stimuli are being tested...\n',i)
    tic                                                         %Start time for the calculation of all combinations with i external stimuli
    Combis=nchoosek(1:OCP.numControls,i);                       %Making all the combinations of i different controls out of numControls many external stimuli
    B=size(Combis);                                             %B(1,1) number of different combinations, B(1,2)=i lenght of the combination
    J=zeros(B(1,1),1);                                          %Initialize vector for the targetfunctional value for each combination
    numCon=B(1,2);                                              %Bufferen due to the parallelization with parfor
    for k=1:max(size(valu))                                     %Loop over all values for each combination
        Valu=valu(k);                                           %k-th value of the interval into which the range of the external stimuli is devided into
        parfor j=1:B(1,1)                                       %Loop over all combination; parfor causes parallel calculation of the target functional values of corresponding combination
            B=size(Combis);                                     %B(1,1) number of different combinations, B(1,2)=i lenght of the combination, has to be defined again because of parallel computing
            u=u0;                                               %Intialize u
            u(Combis(j,:),:)=Valu*ones(B(1,2),Nt);              %Activate just the external stimuli corresponding to combination combis(j,:) by setting them to the value Valu for all time steps
            x=forward(f,u,OCP);                                 %Calculating the state x corresponding to the external stimuli u
            J(j)=0.5*dt*sum(wNoi'.*sum(transpose((x(numNoi,:)-xd).^2)))+dt*OCP.alpha*sum(sum(u(Combis(j,:),:))); %Calculates corresponding target functional value
        end
        [valJ,posMin]=min(J);                                   %Determine smallest target functional value and corresponding index of the entry of J where all the corresponding values of the combination of lenght i are stored for Valu
        if (valJ<Jmin)                                          %Check if there is combination with corresponding lower target functional value than the other combination tested before 
            ures=u0;
            ures(Combis(posMin,:),:)=Valu*ones(numCon,Nt);      %If yes, store time curve of corresponding combination in ures
            Jmin=valJ;                                          %New lowest target functional value valJ
        end
    end
    fprintf('All combinations with %i external stimuli have been tested in %d seconds\n\n',i,toc)
end
if (max(max(ures))==0)                                          %Check if at least one external stimulus is different from zero
    fprintf('Result from combinatorial method is no active external stimulus\n');
else
fprintf('Combination with smaller target functional value than inactive external stimuli\n');
    constant_value=0;
     for l=1:OCP.numControls
        value_of_stimulus=ures(l,1);                            %external stimulus l is a constant function. Thus it is alright to take the first entry to figure out the value of the function
        if(value_of_stimulus>0)                                 %If the external stimulus l is an active one, print l out
            fprintf('%i,',l);
            constant_value=value_of_stimulus;
        end
    end
    fprintf('\n');
fprintf('The constant value of all external stimuli is %d\n',constant_value);
fprintf('Combination of exteranl stimuli found by combinatorial method has target functional value J=%d\n',Jmin);
end

end

