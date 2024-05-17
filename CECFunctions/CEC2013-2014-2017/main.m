clear all
% mex cec17_func.cpp -DWINDOWS
func_num=1;
D=30;
Xmin=-100;
Xmax=100;
pop_size=30;
iter_max=3000;
runs=1;

%CEC2013 28
fhd=str2func('cec13_func'); Fnum = 28;
%CEC2014 30
% fhd=str2func('cec14_func'); Fnum = 30;
%CEC2017 29 去掉第二个
% fhd=str2func('cec17_func'); Fnum = 30;
for i=1:Fnum
    func_num=i;
    for j=1:runs
        [gbest,gbestval,FES]= PSO_func(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);
        xbest(j,:)=gbest;
        fbest(i,j)=gbestval;
        fbest(i,j)
    end
    f_mean(i)=mean(fbest(i,:));
end



% for i=1:29
% eval(['load input_data/shift_data_' num2str(i) '.txt']);
% eval(['O=shift_data_' num2str(i) '(1:10);']);
% f(i)=cec14_func(O',i);i,f(i)
% end