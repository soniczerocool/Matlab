function make_txt_file_stats(precision,recall,dice,jaccard, result_path,name)
folder = [result_path,'textfiles_stats\'];
mkdir(folder)


save_path = [folder, 'statistics','_',name,'.txt'];


f = fopen(save_path,'w');


 fprintf(f,'Precision: ')
for p=1:size(precision,2)
    fprintf(f,' ');
    fprintf(f,num2str(precision(p)));
   

end
 fprintf(f,'\n');
  fprintf(f,'recall: ')
for p=1:size(recall,2)
    fprintf(f,' ');
    fprintf(f,num2str(recall(p)));
   

end
 fprintf(f,'\n');
  fprintf(f,'dice: ')
 for p=1:size(dice,2)
       fprintf(f,' ');
    fprintf(f,num2str(dice(p)));


end
 fprintf(f,'\n');
 fprintf(f,'jaccard: ')
 for p=1:size(jaccard,2)
    fprintf(f,' ');
    fprintf(f,num2str(jaccard(p)));
    

end
 fprintf(f,'\n');
 
fclose(f);


