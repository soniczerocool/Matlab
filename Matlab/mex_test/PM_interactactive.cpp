MInteractiveGraphCutSegmentation(MImage &mask, float sigma)

        GCoptimizationGeneralGraph gc(MXS*MYS,2);

        /* Connect source and sink with mask pixels 128 and 255 */
        for( int y=0;y<MYS;y++)
                for( int x=0;x<MXS;x++) {

                        int pixelIndex = y*MXS+x;
                        if ((int)(mask.MGetColor(x,y,0))==255){
                                gc.setDataCost(pixelIndex,0,1000000.0);
                        }else if ((int)(mask.MGetColor(x,y,0))==128){
                                gc.setDataCost(pixelIndex,1,1000000.0);
                        }

                }
        /* Connect horizontal neighboring pixels */
        for (int y = 0; y < MYS; y++ ){
                for (int  x = 1; x < MXS; x++ ){

                        float gradient = sqrt(pow(MGetColor(x,y,0)-MGetColor(x-1,y,0),2)+pow(MGetColor(x,y,1)-MGetColor(x-1,y,1),2)+pow(MGetColor(x,y,2)-MGetColor(x-1,y,2),2));
                        gc.setNeighbors(x+y*MXS,x-1+y*MXS,(double)(1*exp(-gradient/sigma)));            }
        }

        /* Connect vertical neighboring pixels */
        for (int y = 1; y < MYS; y++ ){
                for (int  x = 0; x < MXS; x++){

                        float gradient = sqrt(pow(MGetColor(x,y,0)-MGetColor(x,y-1,0),2)+pow(MGetColor(x,y,1)-MGetColor(x,y-1,1),2)+pow(MGetColor(x,y,2)-MGetColor(x,y-1,2),2));
                        gc.setNeighbors(x+y*MXS,x+(y-1)*MXS,(double)(1*exp(-gradient/sigma)));
                }
        }

        /* go graph cut! */
        gc.expansion(2);

        /* save results */
        int i=0;
        for( int y=0;y<MYS;y++)
                for( int x=0;x<MXS;x++,i++)
                        MSetColor(gc.whatLabel(i)*255,x,y,0);

