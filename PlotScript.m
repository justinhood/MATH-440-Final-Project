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
%mov(1:depth) = struct('cdata',[], 'colormap',[]);
a=linspace(0,MAX,sizer);
b=linspace(0,MAX,sizer);

for t=1:depth
    surf(a,b,x(:,:,t))
    pause(0.01)
    %zlim([-1 3])
    %mov(t)=getframe(1);
end
%close(1)
%movie2avi(mov, 'ExplicitAnimation.avi', 'compression', 'None', 'fps', 50);
