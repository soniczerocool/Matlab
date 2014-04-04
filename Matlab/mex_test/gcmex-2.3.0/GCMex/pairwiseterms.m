function [gradient2,gradient3,gradient4] = pairwiseterms(image,sigma,T1C,z,x,y)
[MXS,MYS,MZS]=size(T1C);
Num_features= size(image,1);
 pixelIndex1 = z*MXS*MYS+x*MYS+y;
 pixelIndex2 = z*MXS*MYS+x*MYS+ y-1;
 pixelIndex3 = z*MXS*MYS+ (x-1)*MYS+y;
 pixelIndex4 = (z-1)*MXS*MYS+x*MYS+y;   
 
 gradient2 = 0;
 gradient3 = 0;
gradient4 = 0;
                for  f=1 : Num_features
                
                     featureIndex1 = pixelIndex1 * Num_features + f;
                    featureIndex2 = pixelIndex2 * Num_features + f;
                    featureIndex3 = pixelIndex3 * Num_features + f;
                    featureIndex4 = pixelIndex4 * Num_features + f;
                    gradient2 =gradient2 + (image(featureIndex1) - image(featureIndex2))^2;
                    gradient3 =gradient3 + (image(featureIndex1) - image(featureIndex3))^2;
                    gradient4 =gradient4 + (image(featureIndex1) - image(featureIndex4))^2;
                end
                    gradient2 = sqrt(gradient2);
                    gradient3 = sqrt(gradient3);   
                    gradient4 = sqrt(gradient4);     
                  %gradient2=  (exp(-gradient2/sigma))  ;    
                  %  gradient3=  (exp(-gradient3/sigma))  ;       
                  %    gradient4=  (exp(-gradient4/sigma)) ;        