function [ y ] = moving_average( x, win_size )
y1=filter(ones(1,win_size/2+1)/win_size,1,x);
y2=filter(ones(1,win_size/2+1)/win_size,1,fliplr(x));
y=y1+fliplr(y2)-(1/win_size)*x;
end
