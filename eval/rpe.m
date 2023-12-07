function Perr = rpe(idx,Q_g,Q_est,P_g,P_est)

Perr=zeros(size(idx,1),1);

for i=1:size(Perr,1)
    R_g1=quat2rotm(Q_g(idx(i,1),:));
    R_est1=quat2rotm(Q_est(idx(i,1),:));
    R_g2=quat2rotm(Q_g(idx(i,2),:));
    R_est2=quat2rotm(Q_est(idx(i,2),:));
    P_g1=P_g(idx(i,1),:);
    P_est1=P_est(idx(i,1),:);
    P_g2=P_g(idx(i,2),:);
    P_est2=P_est(idx(i,2),:);
    
    T_g1=[R_g1 P_g1'; 0 0 0 1];
    T_g2=[R_g2 P_g2'; 0 0 0 1];  
    T_est1=[R_est1 P_est1'; 0 0 0 1];
    T_est2=[R_est2 P_est2'; 0 0 0 1];   
    
    T=(T_g1\T_g2)\(T_est1\T_est2);
    
    Perr(i)=norm(T(1:3,4));
end

end

