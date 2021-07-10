%% generate projection data
clc;
I = phantom(256);
size_picture = 256;
D = 500;
dtheta = 0.1;
dtheta_rad = dtheta*pi/180;
fan_increment = 1;
[rawdata,FposArcDeg,Fangles] = fanbeam(I,D,...
                    'FanSensorGeometry','arc',...
                 'FanRotationIncrement',fan_increment,'FanSensorSpacing',dtheta);
             

save arc2.txt rawdata -ascii

%%  convolution
%C++和matlab在多维数组的嵌套上是反的，C++是前面的是高维，matlab后面是高维
%matlab 访问矩阵行坐标在前，列坐标在后，与C++的数据宽（对应列）高（对应行不同）
height = size(rawdata,1);
num_angles = size(rawdata,2);
%step1 系数校正
factor = linspace(cos(-height/2*dtheta_rad)*D,cos(height/2*dtheta_rad)*D,height)';
factor_rep = repmat(factor,1,size(rawdata,2));
rawdata = rawdata.*factor_rep;
%step2 生成滤波器
% h = zeros(height,1);
% for i=-floor(height/2):1:floor(height/2)
%     if i==0
%         h(i+ceil(height/2))=0.125/(dtheta_rad^2);
%     elseif mod(abs(i),2)==0
%         h(i+ceil(height/2))=0;       
%     else
%         h(i+ceil(height/2))=-0.5/(pi*pi*sin(i*dtheta_rad)*sin(i*dtheta_rad));
%     end
% end
% h = fftshift(h);
% h(2:height)=h(1:height-1);
% h(1)=0.125/(dtheta_rad^2);
m_nHeight = height;
dsensor = 0.5*pi/180;
g = zeros(m_nHeight,1);
for i = -floor(m_nHeight / 2.0): -1
    if (mod(i,2)~=0)
        g(i + (floor(m_nHeight / 2.0) + ceil(m_nHeight / 2.0))+1) = -1 / (2 * (pi * sin(i * dsensor))^2);
    end
end
g(0+1) = 1 / (8 * dsensor * dsensor);
for i = 1: floor(m_nHeight / 2.0)-1
    if (mod(i,2)~=0)
        g(i+1) = -1 / (2 * (pi * sin(i * dsensor))^2);
    end
end
h=g;
%step3卷积
temp = zeros(height);
for i = 1:num_angles
    temp=cconv((fliplr(rawdata(:,i)'))',h,height);
    rawdata(:,i)=temp;
end

%% backprojection
recon = zeros(256,256);
for x = 1:size_picture
    for y = 1:size_picture
        r = ((x-size_picture/2)^2+(y-size_picture/2)^2)^0.5;
        phi = atan2(y-size_picture/2,x-size_picture/2);
        for k = 1:num_angles
            beta = k*pi/180+pi;
            gamma = atan2(r*cos(beta-phi),D+r*sin(beta-phi));
            L2 = (D+r*sin(beta-phi))^2+(r*cos(beta-phi))^2;
            gamma = gamma*180/pi/dtheta;
            if(1<=gamma+height/2+1&&gamma+height/2+1<=189)
                index=gamma+height/2+1-floor(gamma+height/2+1);
                recon(x,y)=recon(x,y)+2*pi/num_angles*(rawdata(floor(gamma+height/2+1),k)*(1-index)+...
                    rawdata(ceil(gamma+height/2+1),k)*index)/L2;
            end
        end
    end

end
%%
imshow(recon,[]);
