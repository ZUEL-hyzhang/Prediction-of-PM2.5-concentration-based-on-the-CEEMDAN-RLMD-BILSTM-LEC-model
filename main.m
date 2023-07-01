%--------------------------------------------------------------------------
%                         Copula���ۼ�Ӧ��ʵ��
%--------------------------------------------------------------------------

%******************************��ȡ����*************************************
% ���ļ�hushi.xls�ж�ȡ����
hushi = xlsread('mode.csv');
X = hushi(:,1);
% ���ļ�shenshi.xls�ж�ȡ����
shenshi = xlsread('mode.csv');
Y = shenshi(:,6);


%****************************����Ƶ��ֱ��ͼ*********************************
[fx, xc] = ecdf(X);
figure;
ecdfhist(fx, xc, 30);
xlabel('x');  % ΪX��ӱ�ǩ
ylabel('f(x)');  % ΪY��ӱ�ǩ
[fy, yc] = ecdf(Y);
figure;
ecdfhist(fy, yc, 30);
xlabel('y');  % ΪX��ӱ�ǩ
ylabel('f(y)');  % ΪY��ӱ�ǩ


%****************************����ƫ�Ⱥͷ��*********************************
% ����X��Y��ƫ��
xs = skewness(X)
ys = skewness(Y)

% ����X��Y�ķ��
kx = kurtosis(X)
ky = kurtosis(Y)


%******************************��̬�Լ���***********************************
% �ֱ����jbtest��kstest��lillietest������X������̬�Լ���
[h,p] = jbtest(X)  % Jarque-Bera����
[h,p] = kstest(X,[X,normcdf(X,mean(X),std(X))])  % Kolmogorov-Smirnov����
[h, p] = lillietest(X)  % Lilliefors����

% �ֱ����jbtest��kstest��lillietest������Y������̬�Լ���
[h,p] = jbtest(Y)  % Jarque-Bera����
[h,p] = kstest(Y,[Y,normcdf(Y,mean(Y),std(Y))])  % Kolmogorov-Smirnov����
[h, p] = lillietest(Y)  % Lilliefors����


%****************************����ֲ�����ֵ*******************************
% ����ecdf������X��Y�ľ���ֲ�����
[fx, Xsort] = ecdf(X);
[fy, Ysort] = ecdf(Y);
% ����spline����������������ֵ����ԭʼ�����㴦�ľ���ֲ�����ֵ
U1 = spline(Xsort(2:end),fx(2:end),X);
V1 = spline(Ysort(2:end),fy(2:end),Y);

% ����ecdf������X��Y�ľ���ֲ�����
[fx, Xsort] = ecdf(X);
[fy, Ysort] = ecdf(Y);
% ��ȡfx��fy�ĵ�2�������һ��Ԫ�أ�������������㴦�ľ���ֲ�����ֵ
fx = fx(2:end);
fy = fy(2:end);

% ͨ������ͷ�����ָ�ԭʼ�����㴦�ľ���ֲ�����ֵU1��V1
[Xsort,id] = sort(X);
[idsort,id] = sort(id);
U1 = fx(id);
[Ysort,id] = sort(Y);
[idsort,id] = sort(id);
V1 = fy(id);


%*******************************�˷ֲ�����**********************************
% ����ksdensity�����ֱ����ԭʼ����X��Y���ĺ˷ֲ�����ֵ
U2 = ksdensity(X,X,'function','cdf');
V2 = ksdensity(Y,Y,'function','cdf');


% **********************���ƾ���ֲ�����ͼ�ͺ˷ֲ�����ͼ**********************
[Xsort,id] = sort(X);  % Ϊ����ͼ����Ҫ����X��������
figure;  % �½�һ��ͼ�δ���
plot(Xsort,U1(id),'c','LineWidth',5);
hold on
plot(Xsort,U2(id),'k-.','LineWidth',2); 
legend('����ֲ�����','�˷ֲ�����', 'Location','NorthWest'); % �ӱ�ע��
xlabel('x');  % ΪX��ӱ�ǩ
ylabel('F(x)');  % ΪY��ӱ�ǩ

[Ysort,id] = sort(Y);  % Ϊ����ͼ����Ҫ����Y��������
figure;  % �½�һ��ͼ�δ���
plot(Ysort,V1(id),'c','LineWidth',5); 
hold on
plot(Ysort,V2(id),'k-.','LineWidth',2);
legend('����ֲ�����','�˷ֲ�����', 'Location','NorthWest'); % �ӱ�ע��
xlabel('x');  % ΪX��ӱ�ǩ
ylabel('F(x)');  % ΪY��ӱ�ǩ


%****************************���ƶ�ԪƵ��ֱ��ͼ*****************************
% ����ksdensity�����ֱ����ԭʼ����X��Y���ĺ˷ֲ�����ֵ
U = ksdensity(X,X,'function','cdf');
V = ksdensity(Y,Y,'function','cdf');
figure;  % �½�һ��ͼ�δ���
% ���Ʊ�Ե�ֲ��Ķ�ԪƵ��ֱ��ͼ��
hist3([U(:) V(:)],[30,30])
xlabel('x');  % ΪX��ӱ�ǩ
ylabel('y');  % ΪY��ӱ�ǩ
zlabel('Ƶ��');  % Ϊz��ӱ�ǩ


%****************************���ƶ�ԪƵ��ֱ��ͼ*****************************
figure;  % �½�һ��ͼ�δ���
% ���Ʊ�Ե�ֲ��Ķ�ԪƵ��ֱ��ͼ��
hist3([U(:) V(:)],[30,30])
h = get(gca, 'Children');  % ��ȡƵ��ֱ��ͼ�ľ��ֵ
cuv = get(h, 'ZData');  % ��ȡƵ��ֱ��ͼ��Z������
set(h,'ZData',cuv*30*30/length(X));  % ��Ƶ��ֱ��ͼ��Z���������任
xlabel('x');  % ΪX��ӱ�ǩ
ylabel('y');  % ΪY��ӱ�ǩ
zlabel('c(u,v)');  % Ϊz��ӱ�ǩ


%***********************��Copula�в����Ĺ���ֵ******************************
% ����copulafit�������ƶ�Ԫ��̬Copula�е�������ز���
rho_norm = copulafit('Gaussian',[U(:), V(:)])
% ����copulafit�������ƶ�Ԫt-Copula�е�������ز��������ɶ�
[rho_t,nuhat,nuci] = copulafit('t',[U(:), V(:)])


%********************����Copula���ܶȺ����ͷֲ�����ͼ************************
[Udata,Vdata] = meshgrid(linspace(0,1,31));  % Ϊ��ͼ��Ҫ�������µ���������
% ����copulapdf��������������ϵĶ�Ԫ��̬Copula�ܶȺ���ֵ
Cpdf_norm = copulapdf('Gaussian',[Udata(:), Vdata(:)],rho_norm);
% ����copulacdf��������������ϵĶ�Ԫ��̬Copula�ֲ�����ֵ
Ccdf_norm = copulacdf('Gaussian',[Udata(:), Vdata(:)],rho_norm);
% ����copulapdf��������������ϵĶ�Ԫt-Copula�ܶȺ���ֵ
Cpdf_t = copulapdf('t',[Udata(:), Vdata(:)],rho_t,nuhat);
% ����copulacdf��������������ϵĶ�Ԫt-Copula�ֲ�����ֵ
Ccdf_t = copulacdf('t',[Udata(:), Vdata(:)],rho_t,nuhat);
% ���ƶ�Ԫ��̬Copula���ܶȺ����ͷֲ�����ͼ
figure;  % �½�ͼ�δ���
surf(Udata,Vdata,reshape(Cpdf_norm,size(Udata)));  % ���ƶ�Ԫ��̬Copula�ܶȺ���ͼ
xlabel('U');  % ΪX��ӱ�ǩ
ylabel('V');  % ΪY��ӱ�ǩ
zlabel('c(u,v)');  % Ϊz��ӱ�ǩ
figure;  % �½�ͼ�δ���
surf(Udata,Vdata,reshape(Ccdf_norm,size(Udata)));  % ���ƶ�Ԫ��̬Copula�ֲ�����ͼ
xlabel('U');  % ΪX��ӱ�ǩ
ylabel('V');  % ΪY��ӱ�ǩ
zlabel('C(u,v)');  % Ϊz��ӱ�ǩ

% ���ƶ�Ԫt-Copula���ܶȺ����ͷֲ�����ͼ
figure;  % �½�ͼ�δ���
surf(Udata,Vdata,reshape(Cpdf_t,size(Udata)));  % ���ƶ�Ԫt-Copula�ܶȺ���ͼ
xlabel('U');  % ΪX��ӱ�ǩ
ylabel('V');  % ΪY��ӱ�ǩ
zlabel('c(u,v)');  % Ϊz��ӱ�ǩ
figure;  % �½�ͼ�δ���
surf(Udata,Vdata,reshape(Ccdf_t,size(Udata)));  % ���ƶ�Ԫt-Copula�ֲ�����ͼ
xlabel('U');  % ΪX��ӱ�ǩ
ylabel('V');  % ΪY��ӱ�ǩ
zlabel('C(u,v)');  % Ϊz��ӱ�ǩ


%**************��Kendall�����ϵ����Spearman�����ϵ��***********************
% ����copulastat�������Ԫ��̬Copula��Ӧ��Kendall�����ϵ��
Kendall_norm = copulastat('Gaussian',rho_norm)
% ����copulastat�������Ԫ��̬Copula��Ӧ��Spearman�����ϵ��
Spearman_norm = copulastat('Gaussian',rho_norm,'type','Spearman')
% ����copulastat�������Ԫt-Copula��Ӧ��Kendall�����ϵ��
Kendall_t = copulastat('t',rho_t)
Kendall = corr([X,Y],'type','Kendall')
Spearman = corr([X,Y],'type','Spearman')


%******************************ģ������*************************************
% ����ecdf������X��Y�ľ���ֲ�����
[fx, Xsort] = ecdf(X);
[fy, Ysort] = ecdf(Y);
% ����spline����������������ֵ����ԭʼ�����㴦�ľ���ֲ�����ֵ
U = spline(Xsort(2:end),fx(2:end),X);
V = spline(Ysort(2:end),fy(2:end),Y);
% ���徭��Copula����C(u,v)
C = @(u,v)mean((U <= u).*(V <= v));
% Ϊ��ͼ����Ҫ�������µ���������
[Udata,Vdata] = meshgrid(linspace(0,1,31));
% ͨ��ѭ�����㾭��Copula�������²���������㴦�ĺ���ֵ
for i=1:numel(Udata)
    CopulaEmpirical(i) = C(Udata(i),Vdata(i));
end

figure;  % �½�ͼ�δ���
% ���ƾ���Copula�ֲ�����ͼ��
surf(Udata,Vdata,reshape(CopulaEmpirical,size(Udata)))
xlabel('U');  % ΪX��ӱ�ǩ
ylabel('V');  % ΪY��ӱ�ǩ
zlabel('Empirical Copula C(u,v)');  % Ϊz��ӱ�ǩ

% ͨ��ѭ�����㾭��Copula������ԭʼ�����㴦�ĺ���ֵ
CUV = zeros(size(U(:)));
for i=1:numel(U)
    CUV(i) = C(U(i),V(i));
end

% ����������ز���Ϊ0.9264�Ķ�Ԫ��̬Copula������ԭʼ�����㴦�ĺ���ֵ
rho_norm = 0.9264;
Cgau = copulacdf('Gaussian',[U(:), V(:)],rho_norm);
% ����������ز���Ϊ0.9325�����ɶ�Ϊ4�Ķ�Ԫt-Copula������ԭʼ�����㴦�ĺ���ֵ
rho_t = 0.9325;
k = 4.0089;
Ct = copulacdf('t',[U(:), V(:)],rho_t,k);
% ����ƽ��ŷ�Ͼ���
dgau2 = (CUV-Cgau)'*(CUV-Cgau)
dt2 = (CUV-Ct)'*(CUV-Ct)
