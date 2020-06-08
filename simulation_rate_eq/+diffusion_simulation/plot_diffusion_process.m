function  output = plot_diffusion_process(allFragPos,allFragLengths,timeInstances)


% ########################################################################
%
% FUNCTION: output =  plot_diffusion_process(allFragPos,allFragLengths,timeInstances)
% 
% Input: allFragPos = cell array with all fragment positions at all times,
% allFragLengths = cell array with all fragment lengths at all times,
% timeInstances = array with all time instances
% 
% Output: plots the obtained barcode-like time evolution of all fragments
% with their lengths 
%
% Description:
%
% ########################################################################



figure; hold on; box on;

for i=1:length(timeInstances)
    
    for j=1:length(allFragPos{i})
        
        plot(timeInstances(i)*ones(1,2),[allFragPos{i}(j)-0.5*allFragLengths{i}(j),...
            allFragPos{i}(j)+0.5*allFragLengths{i}(j)],'k-','LineWidth',2)
        
    end
    
end

xlabel('time')
ylabel('position')



output = 0;



end