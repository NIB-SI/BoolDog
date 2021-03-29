function [] = drawStimuliFromFile(v,filename,numCol)
%Input: row vector with numbers corresponding to the external stimuli to be
%plotted from the file filename which has to be given as a string
%'filename', number numCol number of external stimuli in a row in
%the plot
%Output: plot of the external stimuli

A=importdata(filename);     %Imports the data from filename
timeSteps=A(1,:);           %time steps
u=A(2:end,:);               %Time curves of external stimuli      

numPlots=max(size(v));      %Number of plots
numRows=ceil(numPlots/numCol);   %Rounds to the next greater integer to get number of rows for the plot in which the external stimuli are plotted

image=figure('name','Plot of external stimuli','NumberTitle','off');    %Creates figure for the plot
panel=uipanel('Parent',image,'BorderType','none');
for i=1:numPlots                                                        %Plots each graph of all external stimuli indexed with v(1,i) in a different subplot
    subplot(numRows,numCol,i,'Parent',panel)
    plot(timeSteps,u(v(1,i),:));                                        %Plot of the external stimuli u(v(1,i),:) over timeSteps
    titleplot=sprintf('$u_{%i}$',v(1,i));                               %Creates the title of each subplot consisting of u_i where i is an element of v
    title(titleplot,'Interpreter','latex')
end

end

