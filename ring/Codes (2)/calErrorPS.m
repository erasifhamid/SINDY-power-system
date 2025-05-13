function [Error] = calErrorPS(yFOM,yPOD,yNLMM,yDEIM)
%% description
% gives the Error (absolute and relative) using L1norm, L2norm and Linfnorm
% for four methods: 
% Simulation FOM, POD, NLMM and NLRK
%
%% Syntax: 
% Error.err<norm><Method>Out<numberofOutput><absolute or relative>
%
% Example: Error.errinfNLMMOut2Rel
% relative Error of Output 2 of method NLMM calculated with infinity norm 
%
%% prerequisites
%Error.errorNLMMOut1 = abs(yFOM-yNLMM);
%Error.errorNLRKout1 = abs(yFOM-yNLrk);
Error.errorPODout1 = abs(yFOM-yPOD);
%Error.errorDEIMout1 = abs(yFOM-yDEIM);

%% L1-norm
%absolute Error
Error.err1NLMMout1Abs = norm(Error.errorNLMMOut1,1);
Error.err1PODout1Abs = norm(Error.errorPODout1,1);
%Error.err1NLRKout1Abs = norm(Error.errorNLRKout1,1);
%Error.err1DEIMout1Abs = norm(Error.errorDEIMout1,1);

%relative Error
Error.err1NLMMout1Rel = norm(Error.errorNLMMOut1,1)./norm(yFOM(1,:),1); % err1 = sum(abs(error))
Error.err1PODout1Rel = norm(Error.errorPODout1,1)./norm(yFOM(1,:),1);
%Error.err1NLRKout1Rel = norm(Error.errorNLRKout1,1)./norm(yFOM(1,:),1);
%Error.err1DEIMout1Rel = norm(Error.errorDEIMout1,1)./norm(yFOM(1,:),1);

%% L2-norm
% absolute Error
Error.err2NLMMout1Abs = norm(Error.errorNLMMOut1,2); 
Error.err2PODout1Abs = norm(Error.errorPODout1,2);
%Error.err2NLRKout1Abs = norm(Error.errorNLRKout1,2);
%Error.err2DEIMout1Abs = norm(Error.errorDEIMout1,2);

% relative Error
Error.err2out1 = norm(Error.errorNLMMOut1,2)./norm(yFOM(1,:),2); % err2 = sqrt(error'*error)
Error.err2PODout1 = norm(Error.errorPODout1,2)./norm(yFOM(1,:),2);
%Error.err2NLRKout1 = norm(Error.errorNLRKout1,2)./norm(yFOM(1,:),2);
%Error.err2DEIMout1 = norm(Error.errorDEIMout1,2)./norm(yFOM(1,:),2);

%% Linf-norm
% absolute Error
Error.errInftyNLMMOut1Abs = norm(Error.errorNLMMOut1,Inf);
Error.errInftyPODout1Abs = norm(Error.errorPODout1,Inf);
%Error.errInftyNLRKout1Abs = norm(Error.errorNLRKout1,Inf);
%Error.errInftyDEIMout1Abs = norm(Error.errorDEIMout1,Inf);

% relative Error 
Error.errInftyNLMMOut1 = norm(Error.errorNLMMOut1,Inf)./norm(yFOM(1,:),Inf); % errInfty = max(abs(error))
Error.errInftyPODout1 = norm(Error.errorPODout1,Inf)./norm(yFOM(1,:),Inf);
%Error.errInftyNLRKout1 = norm(Error.errorNLRKout1,Inf)./norm(yFOM(1,:),Inf);
%Error.errInftyDEIMout1 = norm(Error.errorDEIMout1,Inf)./norm(yFOM(1,:),Inf);
end