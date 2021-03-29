function [  ] = drawStimuli( u,numColu,OCP )
%Input: u, see main for details, numColu number of columns in the output
%window
%Output: plot of the active external stimuli if there is at least one

if(max(max(u))~=0)                          %Print the numbers of active external stimuli if there is an external stimulus being different from a constant zero function
    image=figure('name','Plot of external stimuli','NumberTitle','off');    %Creates figure for the plot
    panel=uipanel('Parent',image,'BorderType','none');
    fprintf('\n');
    fprintf('Active external stimuli, that means being different from constant zero function:\n');                      
    numPlots=0;
    for i=1:OCP.numControls                 %Prints out the numbers which correspond to the active external stimuli, can be copied for the function drawStimuli  
        if(max(u(i,:))~=0)
            fprintf('%i,',i);
            numPlots=numPlots+1;
            v(1,numPlots)=i;
        end
    end
    fprintf('\n');
    numRows=ceil(numPlots/numColu);              %Rounds to the next greater integer to get number of rows for the plot in which the external stimuli are plotted
    timeSteps=0:OCP.timeInterval:(round(OCP.timeHorizon/OCP.timeInterval)-1)*OCP.timeInterval;
    for i=1:numPlots                             %Plots each graph of all external stimuli indexed with v(1,i) in a different subplot
    subplot(numRows,numColu,i,'Parent',panel)
    plot(timeSteps,u(v(1,i),:));                 %Plot of the external stimuli u(v(1,i),:) over timeSteps
    titleplot=sprintf('$u_{%i}$',v(1,i));        %Creates the title of each subplot consisting of u_i where i is an element of v
    title(titleplot,'Interpreter','latex')
    end  
else
    fprintf('No active external stimulus\n');
end

end

