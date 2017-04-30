sizer=100;
MAX=10;
depth=1000;

fileId=fopen('test.txt', 'r');
matSize=[sizer,sizer];

formatSpec='';
for i=1:sizer
    formatSpec=strcat(formatSpec, ' %f');
end
formatSpec=formatSpec(2:end);

x=zeros(sizer, sizer, depth);

for i=1:depth
    temp=fscanf(fileId, formatSpec, matSize);
    x(:,:,i)=temp;
end


figure(1)

a=linspace(0,MAX,sizer);
b=linspace(0,MAX,sizer);

for t=1:depth
    s=surf(a,b,x(:,:,t));
    %zlim([-1 3])
    pause(0.01)
end