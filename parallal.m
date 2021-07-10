n=180;%num of views,options 180,36,10
I = phantom(256); 
M = 256;
% 对投影做傅里叶变换
R=radon(I,0:179);
width=length(R);
filter=2*[0:round(width/2-1), width/2:-1:1]'/width;
% 每一列做傅里叶变换
r_fft=fft(R,width); 
% 每一列做傅里叶变换后滤波
r_fft_filter=r_fft.*filter; 
%滤波后反变换(real取实部)
proj_ifft=real(ifft(r_fft_filter));
theta=linspace(0.5*pi,1.5*pi,180);
%反投影到x轴，y轴
fbp = zeros(256); % 假设初值为0
for i = 1:(180/n):180
    rad = theta(i);%弧度
    for x = 1:256
        for y = 1:256
            t = round((x-M/2)*cos(rad)+(y-M/2)*sin(rad));%将每个元素到最接近的整数。
            if t<size(R,1)/2+1 && t>-size(R,1)/2
                fbp(x,y)=fbp(x,y)+proj_ifft(round(t+size(R,1)/2),i);
            end
        end
    end
end
fbp = (fbp*0.5*pi)/n;
figure,imshow(fbp),title('反投影变换后的图像')


