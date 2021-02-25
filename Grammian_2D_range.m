

%%%%%%%%%%%%%%%%%%%%%%%% Range Observability Grammian lower bound in SLAM for a circular trajectory around any point %%%%%%%%%%%%%%%%


%% Setup


s=sym('s','real');%time variable
t_0=sym('t_0','real');%initial time variable
phi_0=sym('phi_0','real');%initial angle
T=sym('T','real');%time window variable
l=sym('l',[2,1],'real'); %landmark position
x_c_l=sym('x_c_l',[2,1],'real'); %difference between the centre of the circle and the position of the landmark
r_c=sym('r_c','positive');%radius of the circle
omega=sym('omega','positive');%angular velocity

assume(r_c>0 & omega >0)


x_l=x_c_l+r_c*[cos(omega*(s-t_0)+phi_0);sin(omega*(s-t_0)+phi_0)];% position trajectory centered at the landmark

v=r_c*omega*[-sin(omega*(s-t_0)+phi_0);cos(omega*(s-t_0)+phi_0)];% velocity trajectory




%% Simplified Gradient computations

H=[x_l(1);x_l(2)];%simplified gradient of the observations

H=combine(H,'sincos');

C=children(expand(H*H'));%simplified Grammian integrand

%% Simplified Grammian computation
C_11=C(1,1);
C_22=C(2,2);
C_12=C(1,2);
Grammian=sym(zeros(2,2));


for i=1:size(C_11{1},2)
    Grammian(1,1)=Grammian(1,1)+int(C_11{1}{i},s,t_0,t_0+T);
    Grammian(2,2)=Grammian(2,2)+int(C_22{1}{i},s,t_0,t_0+T);
    fprintf('iteration count: %d/%d \n',i,size(C_11{1},2));
    

end
fprintf('\n');

for i=1:size(C_12{1},2)
    Grammian(1,2)=Grammian(1,2)+int(C_12{1}{i},s,t_0,t_0+T);
    
   fprintf('iteration count: %d/%d \n',i,size(C_12{1},2));
    

end


Grammian(2,1)=Grammian(1,2);
Grammian=simplify(Grammian);

%% Exact eigenvalues
[V_gram,lambda_gram]=eig(Grammian);

%% Eigenvalues of the limit
R=limit(Grammian/T,T,Inf);

[V_R,lambda_R]=eig(simplify(R));

V_R=simplify(V_R);
lambda_R=simplify(lambda_R);


%% Saving the results
% savefile='Grammian_range_circle_any_centre.mat';
% save(savefile);



