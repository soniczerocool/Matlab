I =  imread('lena.bmp');
I = rgb2gray(I);
I= double(I);

a = rand(512,512);
c = demo2(I,a,200);
imshow(c,[1,255])