function [eigenvalues,N, size_eigenvalue_problem]=tds_region_roots(tds,region,options)
% tds_region_roots:  compute the characteristic roots of one tds in a
% given rectangle via spectral method. The details of this approach can be referenced in [our paper].

% input:
%---"tds" is the standard structure of a time delay system.
%---"region" is a 1X4 vector to define the rectangle: [ left bound, right bound, lower bound, upper bound]
%---"options" is a structure that can be created with the function "tdsrootsoptions1.m"


% output:
%--- eigenvalues include the approximated roots denoted by l0 and the corrected roots by Newton method denoted by l1; 
%--- number of discretization points;
%---size of the eigenvalue problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tds_check_valid(tds);
tds_assert_properties(tds,'retarded');
tds=tds_normalize(tds);
tds=tds_sort(tds);
tds=tds_compress(tds);

ns=length(tds.A{1});
if tds.hA(1)>0
    tds.hA=[0, tds.hA];
    WK={};
    WK{1}=zeros(ns,ns);
    for i=1:1:length(tds.A)
    WK{i+1}=tds.A{i};
    end
    tds.A=WK;
end

if tds.hA==0
    fprintf('Warning: ****this is an ODE****')
    eigenvalues.l0=eig(tds.A{1});
    eigenvalues.l1=eigenvalues.l0;
    N=[];
    size_eigenvalue_problem=size(tds.A{1});
else    
A=tds.A;    
tau=tds.hA;

N_max=floor(options.max_size_eigenvalue_problem/ns-1);

rmini=region(1);
rmaxi=region(2);
imini=region(3);
imaxi=region(4);

max_it=options.newton_max_iterations;
tol=options.root_accuracy;

n=length(A{1});
m=length(tau);

% scale the dde to the system of maximum delay =1
tau_s=tau/tau(m);
K=A;
for i=1:1:length(A)
    K{i}=tau(m)*A{i};
end

realmini=rmini*tau(m);
realmaxi=rmaxi*tau(m);
imagmini=imini*tau(m);
imagmaxi=imaxi*tau(m);
dh=(realmaxi-realmini)/100;
realpart=[realmini:dh:realmaxi];
origi=realpart+1i*((imagmini+imagmaxi)/2);
a_theta=[-0.2462 -0.2266 -0.2053 -0.2302 -0.2326 -0.2335 -0.2362 -0.2421 -0.2463 -0.2604 -0.2656 -0.2749 -0.2919 -0.3030 -0.3140 -0.3265 -0.3491 -0.3664 -0.3892 -0.4204 -0.4548 -0.4855 -0.5339 -0.5872 -0.6491 -0.7354 -0.8478 -0.9930 -1.1800 -1.4448 -1.9414 -2.7149 -3.0292 -2.7272 -2.2329  -1.9425 -1.7025 -1.5274 -1.3828 -1.2578 -1.1550 -1.0718 -0.9977 -0.9340 -0.8771 -0.8243 -0.7852 -0.7499 -0.7131 -0.6879 -0.6574 -0.6386 -0.6108 -0.5970 -0.5816 -0.5640 -0.5540 -0.5413 -0.5341 -0.5273 -0.5187 -0.5147 -0.5130 -0.5133 -0.5123];
b_theta=[0.9124 0.9123 0.9136 0.9165 0.9195 0.9234 0.9285 0.9345 0.9416 0.9501 0.9592 0.9698 0.9818 0.9947 1.0090 1.0249 1.0427 1.0620 1.0833 1.1069 1.1331 1.1614 1.1936 1.2289 1.2685 1.3132 1.3642 1.4231 1.4913 1.5731 1.6783 1.7867 1.8183 1.7507 1.6451 1.5566 1.4801 1.4158 1.3602 1.3108 1.2670 1.2281 1.1935 1.1621 1.1333 1.1073 1.0841 1.0630 1.0433 1.0260 1.0098 0.9958 0.9823 0.9709 0.9603 0.9508 0.9427 0.9356 0.9295 0.9246 0.9202 0.9170 0.9148 0.9136 0.9133];
theta=[0:pi/64:pi];

theta_points1=angle((realmaxi+1i*imagmaxi)-origi);
theta_points2=angle((realmini+1i*imagmaxi)-origi);
pp1_a=spline(theta,a_theta,theta_points1);
pp1_b=spline(theta,b_theta,theta_points1);
r1_points=abs((realmaxi+1i*imagmaxi)-origi);
N1=(r1_points-pp1_a)./pp1_b;
pp2_a=spline(theta,a_theta,theta_points2);
pp2_b=spline(theta,b_theta,theta_points2);
r2_points=abs((realmini+1i*imagmaxi)-origi);
N2=(r2_points-pp2_a)./pp2_b;
N_int=[];
for i=1:1:length(N1)  
  N_int=[N_int,max(N1(i),N2(i))]; 
end
[N,index]=min(N_int);
N=ceil(N);
origi_fix=origi(index);



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if N==0,
    N=min(20,N_max-1);
elseif N>= N_max,
    disp(['Recommended number of discretization points exceeded maximum value.']) 
    afm=ns*(N_max+1);
    disp(['Discretization around s=',num2str(origi_fix/tau(m)),' with ',num2str(N_max+1),' points (size of eigenvalue problem: ',num2str(afm),'x',num2str(afm),').'])
    N=N_max;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
QQ=A;
if N==N_max
    QQ=K;
end
if N~=N_max
    % introduce a shift of the origin
    QQ{1}=K{1}-origi_fix*eye(n);
end
if N~=N_max 
    for i=2:1:m
        QQ{i}=K{i}*exp(-(origi_fix)*tau_s(i));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if N==1
    N=N+1;
end

b=zeros(N+1,N+1);
b(1,:)=4*ones(1,N+1);
b(2,[1:3])=[2 0 -1];
for i=3:1:N
    b(i,[i-1:i+1])=[1/(i-1) 0 -1/(i-1)];
end
b(N+1,[N:N+1])=[1/N 0];
b=b/4;
Bpi_N=kron(b,eye(n));
Sigma_N=eye(n*(N+1));
Sigma_N(1:n,1:n)=zeros(n,n);
for k=1:1:N+1
    for j=1:1:m-1
    Sigma_N(1:n,(k-1)*n+1:k*n)=Sigma_N(1:n,(k-1)*n+1:k*n)+QQ{j+1}*chebpoly(k-1, -2*tau_s(j+1)+1);
    end
    Sigma_N(1:n,(k-1)*n+1:k*n)=Sigma_N(1:n,(k-1)*n+1:k*n)+QQ{1};
end
size_eigenvalue_problem=size(Sigma_N);

u=eig(Sigma_N,Bpi_N);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if N==N_max
    uu=u;
else
uu=u-(-origi_fix)*ones(length(u),1);
end
uu=uu/tau(m);

uu=uu(real(uu)<=(rmaxi+0.1));
uu=uu(real(uu)>=(rmini-0.1));
uu=uu(imag(uu)<=(imaxi+0.1));
uu=uu(imag(uu)>=(imini-0.1));


if N==N_max,
    warning off
end
warning off

ur=uu;
%uf=zeros(length(uu9),1);
uf=[];ul0=[];
for j=1:1:length(ur)
    %keyboard
    l=ur(j);
    QM=zeros(n,n);
    for i=1:1:m
        QM=QM+tau(i)*A{i}*exp(-l*tau(i));
    end
M=l*eye(n)-A{1};
for ii=1:1:m-1
    M=M-A{ii+1}*exp(-l*tau(ii+1));
end
%some approximation might make M containing element NaN or Inf, in order to
%avoid this we make the following codes
if isscalar(unique(isinf(M)))==1
    if isinf(M)==1
        continue;
    end
elseif isscalar(unique(isinf(M)))==0
    continue;
end
if isscalar(unique(isnan(M)))==1
    if isnan(M)==1
        continue;
    end
elseif isscalar(unique(isinf(M)))==0
    continue;
end   

[U,S,V] = svd(M);
v_a=V(:,n);
v=v_a;
F=M;

 for k=1:1:max_it
     J=[F (eye(n)+QM)*v;(v_a)' 0];
     a=J\[-F*v;-((v_a)'*v-1)];
     v=v+a(1:n,1);
     l=l+a(n+1,1);
     QM=zeros(n,n);
    for i=1:1:m
        QM=QM+tau(i)*A{i}*exp(-l*tau(i));
    end
     F=l*eye(n)-A{1};
    for i=1:1:m-1
        F=F-A{i+1}*exp(-l*tau(i+1));
    end
   
    if norm(F*v)<=tol
    break;
    end   
    
 end
if norm(F*v)>=tol
    l=inf;
end


 if isnan(real(l))==0 && isinf(real(l))==0
        if real(l)<rmaxi
            if real(l)>rmini
                if imag(l)<imaxi
                    if imag(l)>imini
                        if abs(l-ur(j))<0.01  
                       uf=[uf;l]; ul0=[ul0; ur(j)];
                        end
                    end
                end
            end
        end
 end

end
warning on

if isempty(uf)
    eigenvalues.l1=[];
    eigenvalues.l0=[];
    fprintf('No characteristic root in the given right half plane is expected>>>')
else
eigenvalues.l0=ul0;
eigenvalues.l1=uf;
end
end
return;

