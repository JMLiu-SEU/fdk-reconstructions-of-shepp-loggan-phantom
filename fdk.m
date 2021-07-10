%% test dicom 
X= dicomread('001.dcm');
Y=(imshow(X,[]));

%% parameters
size = 1100;
num_angles = 360;
D = 5040; %探测器到坐标原点距离
%% picture read  文件目录转到图片
A=zeros(1100,1100,360);  %1100*1100,360层
for i=1:360
    temp = dicomread(strcat(num2str(i,'%03d'),'.dcm'));
    A(:,:,i)=temp;
end

%%  convolution  文件目录转到程序
%step1
for i=1:size
    for j = 1:size
        A(i,j,:)= A(i,j,:)*D/(D^2+(j-size/2)^2+(size/2-i)^2);
    end
end
%step2 每一行数据 卷积
data_filter = zeros(1100,1100,360); % 卷积后数据
filter = [0:(size/2-1), size/2:-1:1]/size;%ramp
filter_rep = repmat(filter,size,1);
for i=1:num_angles
    row_fft = fft(A(:,:,i),[],2);
    data_filter(:,:,i)=row_fft.*filter_rep;
    data_filter(:,:,i) = real(ifft(data_filter(:,:,i),[],2));
end

%% step3 backprojection
for z = -3:1:3   %-size/2:1:size/2 完整重建
    center = zeros(size); % 初值为0 
    for x = 1:size
        for y = 1:size
            r = ((x-size/2)^2+(y-size/2)^2)^0.5;
            phi = atan2(y-size/2,x-size/2);
            for k = 1:num_angles
                beta = k*pi/180;
                temp1 = D*r*cos(beta-phi)/(D+r*sin(beta-phi));
                temp2 = (D*z)/(D+r*sin(beta-phi));
                scale = D^2/(2*(D+r*sin(beta-phi))^2);
                if(-550<temp1&&temp1<549&&-550<temp2&&temp2<549)
                    center(x,y)=center(x,y)+data_filter(round(temp2+size/2+1),round(temp1+size/2+1),k)*scale;
                end
            end
        end
    end   
    name = strcat(num2str(z+size/2,'%03d'),'.png');
    path = strcat('./reconstruction/',name);
    imwrite(center,path);
end
%imshow(center,[]);


%% normalize
function OutImg = normalize(InImg)
    ymax=255;ymin=0;
    xmax = max(max(InImg)); %求得InImg中的最大值
    xmin = min(min(InImg)); %求得InImg中的最小值
    OutImg = round((ymax-ymin)*(InImg-xmin)/(xmax-xmin) + ymin); %归一化并取整
end














