clear
%survivor strategy： 新一代种群中,一部分由精英策略选取最优的一些个体，另一部分随机选取个体
%选择算子比例根据 当前解的最佳值和平均值的一些关系 进行自适应条件
%停止准则是迭代次数
Nins=[20,23];
for n_i=1:2
 c11=Nins(n_i);
eval(['load C:\Users\dell\Documents\MATLAB\测试算例\instance',num2str(c11)])
num_run=20;
% num_run=1;
num_iter=2500;
Bestfitness1=zeros(1,num_run);
time1=zeros(1,num_run);
Gbest_val=zeros(num_iter,num_run);
% for c=1:num_run
parfor c=1:num_run
%初始化
%load('C:\Users\dell\Documents\MATLAB\data_4_swt_vs_rch\test_case\s=25,w=20,t=30.mat')
tic;
num_w=capa_w; %  武器能力
num_t=cons_st; % 目标约束
num1=num_t;
q=1;
n1=sum(num_w);
n2=sum(num_t);
n=max(max(n1,n2),max(n2,s));%编码维数%
%% 确定时间限制
% if n<10
%     t_lim=0.6*n;%暂用
% else
%     t_lim=0.2658*n^2 - 6.443*n + 47.166; %运行时间限制 %曲线的对称轴在右半平面，n较小时反而时间过长 
% % t_lim = 0.0019*n^3 + 0.0782*n^2 - 0.9942*n + 7.8805 %新拟合的，t_lim随n的增大而增大 
% end 
u=round((6.9713*n+20.687)/2);%随机生成（或者基于贪婪规则）的解的数目
%u=round((7.9713*n+20.687)/2);%随机生成（或者基于贪婪规则）的解的数目
M=2*u;
%M=2*u+1;
% T=500;%进化代数
T=num_iter;
k=2;%锦标赛竞争规模
pc=0.92;%交叉概率
pm=0.05;%变异概率
num_ls=M;
Mof=round(M*pc);

num_trunc=zeros(1,T); %每一代由精英机制筛选出的个体，初始值为50%

ind_restart=zeros(1,T); %记录是否重启动
div_p=zeros(1,T); %记录每一代种群中不同的个体数目
% num_trunc=zeros(1,T); 
% all_temp=zeros(1,T); 
% num_temp=zeros(1,T); 

ls_st=zeros(num_ls,n);
ls_wt=zeros(num_ls,n);
popfit_ls=zeros(1,num_ls); 

popst_all=zeros(Mof+M,n);
popwt_all=zeros(Mof+M,n);
flag_s=zeros(1,n);
flag_w=zeros(1,n);
flag_t=zeros(1,n);
% var_pop=zeros(1,T);
% std_pop=zeros(1,T);
for i=1:s
    flag_s(i)=i;
end
for i=1:w
    flag_w(i)=i;
end
while(sum(num_t)>0)
   for i=1:t
       for j=1:num1(i)
           if num_t(i)~=0
               flag_t(q)=i;
               num_t(i)=num_t(i)-1;
               q=q+1;
           end
       end
   end
end

t1=clock;
MaxFval=zeros(T,1);
bft=0;
flag_ls=0;
%初始化种群

newpop_st=zeros(M,n);
newpop_wt=zeros(M,n);
son_st=zeros(Mof,n);
son_wt=zeros(Mof,n);
popfit_son=zeros(1,Mof);

%随机初始化种群
%   [pop_st,pop_wt,popfit]=rand_ini(M,P1,P2,v,n,t,flag_s,flag_w,flag_t);
% 
%随机+启发式构造
  [pop_st,pop_wt,popfit]=multi_ini(s,w,t,M,P1,P2,v,n,flag_s,flag_w,flag_t,cons_st,cons_wt);


  bft=max(popfit);
%{
%25,20,30 MRBCH构造的解
pr1=[26,24,27,28,10,29,16,30,14,31,32,9,12,33,6,19,3,34,35,1,13,36,37,38,7,8,39,5,40,41,20,42,43,44,25,2,45,23,11,18,21,46,17,47,48,15,49,50,4,22];
pr2=[21,8,22,23,4,24,1,25,10,26,27,3,9,28,15,29,7,30,31,5,32,33,34,35,19,36,37,20,38,39,2,40,41,42,18,14,43,12,16,44,11,45,6,46,47,17,48,49,13,50]; 
pop_st(u,:)=pr1;
pop_wt(u,:)=pr2;
popfit(M)=value_benefit1(pr1,pr2,P1,P2,v,n,t,flag_s,flag_w,flag_t);
%}

  
%迭代过程

for m=1:T
     %选择参与交叉变异的父代个体,记录父本的适应值
    [newpop_st,newpop_wt,popfit_1]=tournament1(pop_st,pop_wt,popfit,k); %规模为k
    
    %交叉 产生M*pc个后代个体
    for i=1:Mof
         a1=floor(rand*M)+1;
         b1=floor(rand*M)+1; 
         if a1==b1
              b1=floor(rand*M)+1; 
         end
         son_st(i,:)=cross_cx(newpop_st(a1,:),newpop_st(b1,:));  
         a2=floor(rand*M)+1;
         b2=floor(rand*M)+1; 
         if a2==b2
              b2=floor(rand*M)+1; 
         end
         son_wt(i,:)=cross_cx(newpop_wt(a2,:),newpop_wt(b2,:));
        
    end
   %变异
   for i=1:round(M*pm)
      s_m=randi(Mof);
      pr1=son_st(s_m,:);
      pr2=son_wt(s_m,:);
      son_st(s_m,:)=multi_mutate(pr1);
      son_wt(s_m,:)=multi_mutate(pr2);
   end
    %计算子代个体适应值 
    for ii=1:Mof
        pr1=son_st(ii,:);
        pr2=son_wt(ii,:);
        popfit_son(ii)=value_benefit1(pr1,pr2,P1,P2,v,n,t,flag_s,flag_w,flag_t);
    end
    
   %survivor selection
   popfit_all=[popfit,popfit_son];
   popst_all=[pop_st;son_st];
   popwt_all=[pop_wt;son_wt];
   [~,ind]=sort(-popfit_all);  
%    num_trunc(m)=fix(max(0.02*M,p_trunc(m)*M));
    num_trunc(m)=max(round(0.01*M),M-round(0.00075*M*max(popfit)/(max(popfit)-mean(popfit))));
    
    
    %
%    popfit_1=zeros(1,M);
      for i=1:M
          if i<=num_trunc(m)      
           %truncation selection
             newpop_st(i,:)=popst_all(ind(i),:);
             newpop_wt(i,:)=popwt_all(ind(i),:);
             popfit_1(i)=popfit_all(ind(i));
          else
               %random selection
               rp1=randi(M);
               newpop_st(i,:)=popst_all(rp1,:);
               newpop_wt(i,:)=popwt_all(rp1,:);
               popfit_1(i)=popfit_all(rp1);
%              [newpop_st(i,:),newpop_wt(i,:),popfit_1(i)]=tournament11(pop_st,pop_wt,son_st,son_wt,popfit,popfit_son);
          end
      end
     %}
      pop_st=newpop_st;
      pop_wt=newpop_wt;
      popfit=popfit_1;
      

%       div_p(m)=length(unique(popfit)); %种群中不同fitness的个数
%       p_trunc(m+1)=(120.^(-div_p(m)/M)-120.^(-1))/(1-120.^(-1));    
%        p_trunc(m+1)=(exp(-div_p(m)/M)-exp(-1))/(1-exp(-1));    

    [MaxFit,MaxInd]=max(popfit);
    if MaxFit>bft
        bft=MaxFit;
    else
        flag_ls=1;
    end
    MaxFval(m)=bft;
   
    %% Local Search
    if flag_ls==1  
     for j=1:num_ls
        [temp_st,temp_wt]=LS(pop_st(MaxInd,:),pop_wt(MaxInd,:));
        pr1=temp_st;
        pr2=temp_wt;
        temp_fit=value_benefit1(pr1,pr2,P1,P2,v,n,t,flag_s,flag_w,flag_t);
        [Minfit,MinInd]=min(popfit);
         if temp_fit>Minfit
             pop_st(MinInd,:)=temp_r1;
             pop_wt(MinInd,:)=temp_r2;
             popfit(MinInd)=temp_fit; 
         end  
    end		 
     flag_ls=0;
    end
 
%{
    if flag_ls==1  
     for j=1:num_ls
       l=1;
       [temp_st,temp_wt]=LS(pop_st(MaxInd,:),pop_wt(MaxInd,:));
        pr1=temp_st;
        pr2=temp_wt;
        temp_fit=value_benefit1(pr1,pr2,P1,P2,v,n,t,flag_s,flag_w,flag_t);
        if temp_fit>popfit(MaxInd)
           ls_st(l,:)=temp_st;
           ls_wt(l,:)=temp_wt;
           popfit_ls(l)=temp_fit; 
           l=l+1;
        end        
     end
      
     flag_ls=0;
     %更新
      [fit_max,pos_ls]=max(popfit_ls);
      if fit_max>bft
         pop_st(MaxInd,:)=ls_st(pos_ls,:);
         pop_wt(MaxInd,:)=ls_wt(pos_ls,:);
         popfit(MaxInd)=fit_max; 
      end
    end 
    %}
    [MaxFit_LS,~]=max(popfit);
    if MaxFit_LS>bft
        bft=MaxFit_LS;
    end
    MaxFval(m)=bft;
%    var_pop(m)=var(popfit,0);
%     std_pop(m)=std(popfit*10^10,0);
%      fit_temp=popfit*10^10;
%      std_temp=std(fit_temp,0)/(10^10);
%      std_pop(m)=std_temp; %标准差
%      std_pop(m)=std_temp/mean(popfit);%标准差变异系数
     div_p(m)=length(unique(popfit)); %种群中不同fitness的个数
     %重启动
   
%     if std_pop(m)<1*10^(-3) 
      if div_p(m)<2
%       if div_p(m)<0.01*M
          %随机初始化种群
%          [pop_st,pop_wt,popfit]=rand_ini(M,P1,P2,v,n,t,flag_s,flag_w,flag_t);
         ind_restart(m)=1;
      %随机+启发式构造
          [pop_st,pop_wt,popfit]=multi_ini(s,w,t,M,P1,P2,v,n,flag_s,flag_w,flag_t,cons_st,cons_wt);

      end  
end

[x1,y1]=sort(popfit);
Bestfitness1(c)=bft;
Gbest_val(:,c)=MaxFval;
%figure(c);
% plot(MaxFval(1:m-1),'-o');
Xopt=pop_st(y1(M),:);
Yopt=pop_wt(y1(M),:);
t2=clock;
time1(c)=etime(t2,t1);
% c=c+1;
end
fave=mean(Bestfitness1);
fmax=max(Bestfitness1);

%% 数据保存
% serial running
%  eval(['save C:\Users\dell\Desktop\memetic4s-wta\data\memetic0821tlim\instance',num2str(c11)]);
 % parallel running
 eval(['save C:\Users\dell\Desktop\memetic4s-wta\绘图用数据\MA_ins',num2str(c11)]);
%  eval(['save C:\Users\dell\Desktop\memetic4s-wta\data\memetic0821tlim\par_instance',num2str(c11)]);
end
