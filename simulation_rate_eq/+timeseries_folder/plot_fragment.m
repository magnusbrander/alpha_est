function output = plot_fragment(fragments,kymograph,figureNr)

figure(figureNr);
hold on;
imshow(kymograph);


for i=1:size(fragments,1)
    
    
    for j=1:size(fragments,2)
        
        
        if ~isempty(fragments{i,j})
            
            f = fragments{i,j};
            
            pos = f(:,1);
            time = f(:,7);
            figure(figureNr);
            hold on;
            line(pos,time);
            
        end
        
        
    end
    
 
end



   output = 0;

   
end

