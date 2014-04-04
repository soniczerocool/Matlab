image = importance_map';
importance_map = importance_map(:);
[Nb_points,Nb_labels] = size(image);
   for pixelIndex = 1:Nb_points
       
              for label=1: Nb_labels
              
                  index = pixelIndex *(Nb_labels) + label;
                   Data_cost = k- importance_map(index)
                   pause
                
              end
   end
  