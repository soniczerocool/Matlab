function avg = write_txt_file_stats(precision,recall,dice,jaccard, f)
if length(precision)~=4
    error('precision must be a vector not a scalor')
end
avg = 0;
 fprintf(f,'Precision: ');
for p=1:size(precision,2)
    fprintf(f,' ');
    fprintf(f,num2str(precision(p)));
   avg = avg + precision(p);

end
 fprintf(f,'\n');
  fprintf(f,'recall: ');
for p=1:size(recall,2)
    fprintf(f,' ');
    fprintf(f,num2str(recall(p)));
      avg = avg + recall(p);


end
 fprintf(f,'\n');
  fprintf(f,'dice: ');
 for p=1:size(dice,2)
       fprintf(f,' ');
    fprintf(f,num2str(dice(p)));
       avg = avg + dice(p);


end
 fprintf(f,'\n');
 fprintf(f,'jaccard: ');
 for p=1:size(jaccard,2)
    fprintf(f,' ');
    fprintf(f,num2str(jaccard(p)));
       avg = avg + jaccard(p);


 end


 fprintf(f,'\n');
 
 avg = avg/16;

 fprintf(f,'average: ');
  fprintf(f,num2str(avg));
 fprintf(f,'\n');


