function [eigenvalues,N, size_eigenvalue_problem]=tds_charateristic_roots(tds, options)
% tds_charateristic_roots: compute the characteristic roots of one tds in a given right half plane via spectral method with heuristic. 
%The details of this approach can be referenced in [our paper].

% input:
%--"tds" is a standard structure of one time delay system. see tds_create.m;
%--"options" is a structure that can be created with the function tdsrootsoptions.m

% output:
%---eigenvalues include the approximated roots denoted by eigenvalues.l0 and the corrected roots by Newton method denoted by eigenvalues.l1; 
%--- number of discretization points N;
%---size of the discretized eigenvalue problem.
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
if isempty(options.minimal_real_part)
    options.minimal_real_part=-1/(max(tau));
end
r=options.minimal_real_part;

max_it=options.newton_max_iterations;
tol=options.root_accuracy;
tm=[];
if isempty(options.commensurate_basic_delay)==1
    tm=0;
else
tau_com=tau/options.commensurate_basic_delay;
for i=1:1:length(tau_com)
if abs(tau_com(i)-round(tau_com(i)))<1e-10
    tm=[tm, tau_com(i)];
end
end
end
n=length(A{1});
m=length(tau);
if length(tm)==length(tau)
   type_of_delays='commensurate delays';
elseif m==2
    type_of_delays='one delay';
elseif m==3
    type_of_delays='two delays';
elseif m==4
    type_of_delays='three delays';
else 
    type_of_delays='more than three delays';
end

  if (strcmp(type_of_delays,'commensurate delays')==0) && (isempty(options.commensurate_basic_delay)==0)
     disp(['delays are not commensurate delays with option.commensurate_basic_delay =', num2str(options.commensurate_basic_delay)])
  end
  
% scale the dde to the system of maximum delay =1
tau_s=tau/tau(m);
K=A;
for i=1:1:length(A)
    K{i}=tau(m)*A{i};
end
r=r*tau(m);
% introduce a shift of the origin
B=K{1}+(-r)*eye(n);
C=zeros(n,n,m-1);
for ii=1:1:m-1
    C(:,:,ii)=K{ii+1}*exp((-r)*tau_s(ii+1));
end
% r is too big or too small can make matrix C in bad conditions.

 if ((max(unique(isnan(C)))==1) || (max(unique(isinf(C)))==1)),
    disp(['WARNING: the value of options.minimal_real_part is extreme.'])
    eigenvalues.l0=[];eigenvalues.l1=[];N=[]; size_eigenvalue_problem=[];
    return;
 end 

magnitude=zeros(1,m-1);
 for ii=1:1:m-1;
     magnitude(ii)=norm(C(:,:,ii));
 end

 maxmag=max(magnitude);

if (maxmag> 10^9) || (maxmag< 10^(-9)),
    disp(['WARNING: the value of options.minimal_real_part is extreme.'])
    disp(['Push a key to continue.'])
    pause
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=options.number_grid_points_p/2; h=pi/p; factor=1.01*sin(pi/p); 
a_theta=[-0.2462 -0.2266 -0.2053 -0.2302 -0.2326 -0.2335 -0.2362 -0.2421 -0.2463 -0.2604 -0.2656 -0.2749 -0.2919 -0.3030 -0.3140 -0.3265 -0.3491 -0.3664 -0.3892 -0.4204 -0.4548 -0.4855 -0.5339 -0.5872 -0.6491 -0.7354 -0.8478 -0.9930 -1.1800 -1.4448 -1.9414 -2.7149 -3.0292];
b_theta=[0.9124 0.9123 0.9136 0.9165 0.9195 0.9234 0.9285 0.9345 0.9416 0.9501 0.9592 0.9698 0.9818 0.9947 1.0090 1.0249 1.0427 1.0620 1.0833 1.1069 1.1331 1.1614 1.1936 1.2289 1.2685 1.3132 1.3642 1.4231 1.4913 1.5731 1.6783 1.7867 1.8183];
theta=[0:pi/64:pi/2];
 
switch lower(type_of_delays)
    case 'one delay'
   gk=[];    
for k=1:1:p+1
    W=B;
    W=W+C(:,:,1)*exp(1i*(k-1)*h);
    f=eig(W)';
    gk=[gk,f];
end

    case 'two delays'
        gk=[];
 for k=1:1:p+1
    for j=-(p-1):1:p,
        W=B;
        W=W+C(:,:,1)*exp(1i*(k-1)*h)+C(:,:,2)*exp(1i*j*h);
        f=eig(W)';
        gk=[gk,f];
    end
end
        
    case 'three delays'
        gk=[];
for k=1:1:p+1
    for j=-(p-1):1:p,
        for s=-(p-1):1:p,
        W=B;
        W=W+C(:,:,1)*exp(1i*(k-1)*h)+C(:,:,2)*exp(1i*j*h)+C(:,:,3)*exp(1i*s*h);
        f=eig(W)';
        gk=[gk,f];
        end
    end
end
    case 'commensurate delays'
        gk=[];
        if tm(end)<=options.new_basic_delay_q
 ba=options.commensurate_basic_delay/tau(m);
% n_k=tau_s/ba
jhh=pi/(p*tm(end));
for k=1:1:p*tm(end)+1
     W=B;
       for l=1:1:m-1
           W=W+C(:,:,l)*exp(1i*(k-1)*jhh*tm(l+1));
       end
       f=eig(W)';
       gk=[gk,f];
end
    elseif tm(end)>options.new_basic_delay_q
            fa=tau_s(m)/options.new_basic_delay_q;
            n_c=tau_s/fa;
            n_cc=round(n_c);
              jhhh=pi/(p*n_cc(end));
    for k=1:1:p*n_cc(end)+1
    W=B;
    for l=1:1:m-1
        W=W+C(:,:,l)*exp(1i*(k-1)*jhhh*n_cc(l+1));
    end
    f=eig(W)';
    gk=[gk,f];
     end
        end
       
    case 'more than three delays'
gk=[];
fa=tau_s(m)/options.new_basic_delay_q;
n_c=tau_s/fa;
n_cc=round(n_c);
jhhh=pi/(p*n_cc(end));
for k=1:1:p*n_cc(end)+1
    W=B;
    for l=1:1:m-1
        W=W+C(:,:,l)*exp(1i*(k-1)*jhhh*n_cc(l+1));
    end
    f=eig(W)';
    gk=[gk,f];
end
end     


si=max(real(gk));

if si<=0 
    N=0;
else %begin alternative
  
    points=[];
    for j=1:1:length(gk),
        if (real(gk(j))>=0) && (real(gk(j))<=factor*si),
            points=[ points gk(j)]; 
        end
    end

    
    theta_points=abs(angle(points));
    pp_a=spline(theta,a_theta,theta_points);
    pp_b=spline(theta,b_theta,theta_points);
    r_points=abs(points);
    N_gk=(r_points-pp_a)./pp_b;
    if isempty(N_gk),
       N1=0;
    else
       N1=max(N_gk);
    end
    
switch lower(type_of_delays)
    case 'one delay'
        gk=[];
for k=1:1:p+1
    W=B;
    W=W+C(:,:,1)*(exp(1i*(k-1)*h)*exp(-factor*si*tau_s(2)));
    f=eig(W)';
    gk=[gk,f];
end
 
    case 'two delays'
        gk=[];
        
    for k=1:1:p+1,
        for j=-(p-1):1:p,
            W=B;
            W=W+C(:,:,1)*exp(1i*(k-1)*h)*exp(-factor*si*tau_s(2))+C(:,:,2)*exp(1i*j*h)*exp(-factor*si*tau_s(3));
            f=eig(W)';
            gk=[gk,f];
        end
    end

    case 'three delays'
        gk=[];
   for k=1:1:p+1,
        for j=-(p-1):1:p,
            for s=-(p-1):1:p,
            W=B;
            W=W+C(:,:,1)*exp(1i*(k-1)*h)*exp(-factor*si*tau_s(2))+C(:,:,2)*exp(1i*j*h)*exp(-factor*si*tau_s(3))+C(:,:,3)*exp(1i*s*h)*exp(-factor*si*tau_s(4));
            f=eig(W)';
            gk=[gk,f];
            end
        end
   end
   
    case 'commensurate delays'

   gk=[];
if tm(end)<=options.new_basic_delay_q
for k=1:1:p*tm(end)+1
    W=B;
    for l=1:1:m-1
          W=W+C(:,:,l)*exp(1i*(k-1)*jhh*tm(l+1))*exp(-factor*si*ba*tm(l+1));
    end
    f=eig(W)';
    gk=[gk,f];
end
elseif tm(end)>options.new_basic_delay_q
    for k=1:1:p*n_cc(end)+1
    W=B;
    for l=1:1:m-1
          W=W+C(:,:,l)*exp(1i*(k-1)*jhhh*n_cc(l+1))*exp(-factor*si*fa*n_cc(l+1));
    end
    f=eig(W)';
    gk=[gk,f];
    end 
end
    
    
    
    case 'more than three delays'

    gk=[];
for k=1:1:p*n_cc(end)+1
    W=B;
    for l=1:1:m-1
          W=W+C(:,:,l)*exp(1i*(k-1)*jhhh*n_cc(l+1))*exp(-factor*si*fa*n_cc(l+1));
    end
    f=eig(W)';
    gk=[gk,f];
end 
end

points=[];
    for j=1:1:length(gk),
        if  real(gk(j))>=factor*si,
        points=[ points gk(j)];       
        end
    end

    theta_points=abs(angle(points));
    pp_a=spline(theta,a_theta,theta_points);
    pp_b=spline(theta,b_theta,theta_points);
    r_points=abs(points);
    N_gk=(r_points-pp_a)./pp_b;
    if isempty(N_gk),
        N2=0;
    else
        N2=max(N_gk);
    end
   
 N=ceil(max(N1,N2));

 end % end alternative
 
 
if N==0,
    N=min(10,N_max-1);
elseif N>= N_max,
    disp(['Recommended number of discretization points exceeded maximum value.']) 
    afm=ns*(N_max+1);
    disp(['Discretization around zero with ',num2str(N_max+1),' points (size of eigenvalue problem: ',num2str(afm),'x',num2str(afm),').'])
    N=N_max;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

QQ=A;
if N==N_max
    QQ=K;
end
if N~=N_max
    QQ{1}=B;
end
if N~=N_max
    for i=2:1:m
        QQ{i}=C(:,:,i-1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if N==1
    N=N+4;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%newton correction part 3>=-1.01*fj_r+fj, which is a little bit board than give
%right half plane, doing this is in order to capture all the roots when
%their approximations are on the left of the line, however their
%corrections are on the right of the line.
if N==N_max
    uu=u;
else
uu=u-(-r)*ones(length(u),1);
end
uu=uu/tau(m);
%    
% apply newton for the roots(uu) of original problem
%

if N==N_max,
    warning off
end
warning off
fj=max(real(uu));
fj_r=abs(fj-(r/tau(m)));
uu1=uu(real(uu)>=-1.01*fj_r+fj);
ur=uu1(imag(uu1)>=0); 
% here we still need to order the approximations by decreasing their real
% part, if we don't do this, some smaller real part roots can converges to
% the larger real part corrections, therefore we lose the right
% approximations.
DS=[real(ur) imag(ur)];
DS=sortrows(DS,-1);
ur=DS(:,1)+1i*DS(:,2);

uf=[10000];ul0=[];
for j=1:1:length(ur)
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
if max(unique(isinf(M)))==1
    continue;
end
if max(unique(isinf(M)))==1
    continue;
end

[U,S,V] = svd(M);
v=V(:,n);
F=M;
vr=real(v);
vi=imag(v);
vn=vr'*vr+vi'*vi;
v_c=(v)/vn; 

 for k=1:1:max_it
     J=[F (eye(n)+QM)*v;(v_c)' 0];
     a=J\[-F*v;-((v_c)'*v-1)];
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

    if norm(a)<=tol
    break;
    end   
    
 end
if norm(a)>=tol
    l=inf;
end
% since we only correct the roots with their imag>=0, so if there
% corrections go to imag<=0, we remove them in order to get rid of the
% repeat roots.
if imag(l)<0
    l=inf;
end
%this for loop is to avoid containing same roots. also this can make the
%newton corrections be reasonable comparing with the inital approximations,
%such that two of them not too far from each other

if (isnan(real(l))==0) && (isinf(real(l))==0)
  for i=1:1:length(uf)
        if abs(l-uf(i))<=1e-6     
            l=inf;
            break
        end
  end
end      
if (isnan(real(l))==0) && (isinf(real(l))==0)
    if real(l)>=(r/tau(m))
        if abs(l-ur(j))<0.01   %to remove the roots if their approximations are far from their corrections.
        uf=[uf;l];ul0=[ul0; ur(j)];
        end
    end
end
end
warning on
uf=uf(2:end);
if isempty(uf)
    eigenvalues.l1=[];
    eigenvalues.l0=[];
    fprintf('No characteristic root in the given right half plane is expected>>>')
else
uu0=[];
for i=1:1:length(ul0)
     if imag(ul0(i))==0
      uu0 =[uu0;ul0(i)];
    elseif imag(ul0(i))~=0 
      uu0=[uu0;ul0(i);conj(ul0(i))];
     end
end
eigenvalues.l0=uu0;
    
uhh=[];    
for i=1:1:length(uf)
    if imag(uf(i))==0
        uhh=[uhh;uf(i)];
    elseif imag(uf(i))~=0 
        uhh=[uhh;uf(i);conj(uf(i))];
    end
end

uhh=uhh(real(uhh)>=(r/tau(m)));

eigenvalues.l1=uhh;
end
end

return;






























