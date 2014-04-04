h_idx = mask_idx == 10;
Et_idx = mask_idx == 4;
nEt_idx = mask_idx == 3;
e_idx = mask_idx == 2;

selected_space(1:3,h_idx)

 h=figure(4);
    plot3(selected_space(1,h_idx),selected_space(2,h_idx),selected_space(3,h_idx),'.g');
    hold on
    plot3(selected_space(1,Et_idx),selected_space(2,Et_idx),selected_space(3,Et_idx),'.r');
    plot3(selected_space(1,nEt_idx),selected_space(2,nEt_idx),selected_space(3,nEt_idx),'.black');
    plot3(selected_space(1,e_idx),selected_space(2,e_idx),selected_space(3,e_idx),'.cyan');
    xlabel('T1C')
    ylabel('T2')
    zlabel('FLAIR')
    print('-depsc','-r300','picture2')
    
    
    h=figure(5);
    plot3(space(1,:),space(2,:),space(3,:),'.b');
   
    xlabel('T1C')
    ylabel('T2')
    zlabel('FLAIR')
 print('-depsc','-tiff','-r300','picture2')