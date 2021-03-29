function [  ] = drawStates(x,numColx,OCP)
%Input: x, see main for details, number of columns in the output window
%Output: plot of the states

timeSteps=0:OCP.timeInterval:round(OCP.timeHorizon/OCP.timeInterval)*OCP.timeInterval;               %time steps
numNoi=OCP.DataNoi(:,1);            %Indices of nodes of interest
numPlots=max(size(numNoi));         %Number of plots
numRows=ceil(numPlots/numColx);     %Rounds to the next greater integer to get number of rows for the plot in which the states are plotted

image=figure('name','Plot of states','NumberTitle','off'); %Creates figure for the plot
panel=uipanel('Parent',image,'BorderType','none');
for i=1:numPlots                                            %Plots each graph of all external stimuli indexed with v(1,i) in a different subplot
    subplot(numRows,numColx,i,'Parent',panel)
    plot(timeSteps,x(numNoi(i),:));                         %Plot of the states x(v(1,i),:) over timeSteps
    titleplot=sprintf('$x_{%i}$',numNoi(i));                %Creates the title of each subplot consisting of x_i where i is an element of v
    title(titleplot,'Interpreter','latex')
end


end

