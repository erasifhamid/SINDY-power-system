function showTable_PS(Error,FOMSnaptime,POD_ROMtime,NLMM_ROMtime,PODtime,NLMM_Time,yFOM,n,rdefl,dm)

method={['FOM (n=',num2str(n),')'];['POD(r=',num2str(rdefl),')'];['S0-NLMM(r=',num2str(rdefl),')','DEIM(m=',num2str(dm),')']};
Basis_time={'-';PODtime+FOMSnaptime;NLMM_Time};
Online_time={'-';POD_ROMtime;NLMM_ROMtime};
Total_CPU_time={FOMSnaptime;PODtime+FOMSnaptime+POD_ROMtime;NLMM_Time+NLMM_ROMtime};
L1_error={'-';Error.err1PODout1Abs;Error.err1NLMMout1Abs};
L2_error={'-';Error.err2PODout1Abs;Error.err2NLMMout1Abs};
Linf_error={'-';Error.errInftyPODout1Abs;Error.errInftyNLMMOut1Abs};



% method={'FOM (n=110)';'POD (r=10)';'NLMM (r=10, K=28)';'NLMM_DEIM(r=10,k=220)'};
% Basis_time={'-';PODtime+FOMSnaptime;NLMM_Time;FOMSnaptime+DEIM_pre_time};
% Online_time={FOMSnaptime;POD_ROMtime;NLMM_ROMtime;NLMM_DEIM_ROMtime};
% L1_error={'-';Error.err1PODout1Abs;Error.err1NLMMout1Abs;Error.err1DEIMout1Abs};
% L2_error={'-';Error.err2PODout1Abs;Error.err2NLMMout1Abs;Error.err2DEIMout1Abs};
% Linf_error={'-';Error.errInftyPODout1Abs;Error.errInftyNLMMOut1Abs;Error.errInftyDEIMout1Abs};

if length(yFOM(:,1)) >=2
    L1errorout2={'-';Error.err1PODout2Abs;Error.err1NLMMout2Abs;Error.err1NLRKout2Abs};
    L2errorout2={'-';Error.err2PODout2Abs;Error.err2NLMMout2Abs;Error.err2NLRKout2Abs};
    Linferrorout2={'-';Error.errInftyPODout2Abs;Error.errInftyNLMMOut2Abs;Error.errInftyNLRKout2Abs};
    T=table(method,Basis_time,Online_time,Total_CPU_time,L1_error,L2_error,Linf_error,L1errorout2,L2errorout2,Linferrorout2)
else
    T=table(method,Basis_time,Online_time,Total_CPU_time,L1_error,L2_error,Linf_error)
end

end

% method={'FOM (n=110)';'POD (r=10)';'NLMM_DEIM'};
% Basis_time={'-';PODtime;NLMM_Time};
% Online_time={FOMSnaptime;POD_ROMtime;NLMM_DEIM_ROMtime};
% L1_error={'-';Error.err1PODout1Abs;Error.err1DEIMout1Abs};
% L2_error={'-';Error.err2PODout1Abs;Error.err2DEIMout1Abs};
% Linf_error={'-';Error.errInftyPODout1Abs;Error.errInftyDEIMout1Abs};
% 
% if length(yFOM(:,1)) >=2
%     L1errorout2={'-';Error.err1PODout2Abs;Error.err1NLMMout2Abs;Error.err1NLRKout2Abs};
%     L2errorout2={'-';Error.err2PODout2Abs;Error.err2NLMMout2Abs;Error.err2NLRKout2Abs};
%     Linferrorout2={'-';Error.errInftyPODout2Abs;Error.errInftyNLMMOut2Abs;Error.errInftyNLRKout2Abs};
%     T=table(method,Basis_time,Online_time,L1_error,L2_error,Linf_error,L1errorout2,L2errorout2,Linferrorout2)
% else
%     T=table(method,Basis_time,Online_time,L1_error,L2_error,Linf_error)
% end

%end