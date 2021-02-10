function Morphology_Method(f)

se_A = strel('disk',21,0);
A = imopen(f,se_A);
figure, imshow(A),title('Seprate The Largest Circle');


ff= f-A;
se_B = strel('disk',18,0);
B = imopen(ff,se_B);
figure,imshow(B),title('Seprate 2nd Circle');


fff = ff - B;
se_C = strel('square',10);
C = imopen(fff,se_C);
figure,imshow(C),title('Seprate Square');


end

