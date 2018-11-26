I1 = imread('frame_0002.png');
I1_C1 = I1(:,:,1);
I1_C2 = I1(:,:,2);
I1_C3 = I1(:,:,3);
I2 = imread('frame_0003.png');
I2_C1 = I2(:,:,1);
I2_C2 = I2(:,:,2);
I2_C3 = I2(:,:,3);


size_I1 = size(I1);
hh = size_I1(1);
ww = size_I1(2);
cc = size_I1(3);

sum_diff = zeros(hh,ww);
perm = zeros(hh,ww);

theta = 0.005;
alpha = 1.0;

for i = 1 : hh-1
    for j = 1 : ww-1
        for c = 1 : cc
            sum_diff(i,j) = sum_diff(i,j) + (double(I1(i,j,c)) - double(I1(i,j+1,c)))^2;
        end
        perm(i,j) = (1 + (abs(sqrt(sum_diff(i,j) / 3) / (sqrt(3) * theta)))^alpha)^-1;
    end
end

%temp_perm = perm*500;
figure, imshow(perm*500,[]);
%figure, imshow(perm*500);

