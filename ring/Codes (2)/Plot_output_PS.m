function Plot_output_PS(tFOM,tPOD,tNLMM,yFOM,yPOD,yNLMM)

h=figure
hold on
grid on

plot(tFOM, yFOM(1,:),'color',[1,0.78,0.054],'linewidth',6)  %yellow
plot(tPOD,yPOD(1,:),':','color',[0.1215,0.145,0.46],'linewidth',6)  %blue
plot(tNLMM, yNLMM(1,:),'--','color',[1,0.45,0],'linewidth',6)  %orange
%plot(tFOM, yFOM(1,:),'r','Linewidth',4)       %'y','Linewidth',6
%plot(tPOD, yPOD(1,:),'g-.','Linewidth',4)     %'k--','Linewidth',6
%plot(tNLMM, yNLMM(1,:),'b-.','Linewidth',4)  
%plot(tNLrk, yNLrk(1,:),'b-.','Linewidth',2)
%plot(tNLMM, yDEIM(1,:),'b--','Linewidth',4)     %'r-.','Linewidth',2
xlabel('$t$','FontSize',25,'Interpreter','latex');
ylabel('$y(t)$','FontSize',25,'Interpreter','latex');
grid on; box on;
legend('FOM','POD','SO-NLMM-DEIM','FontName', 'cmr12','Interpreter','latex')
set(gca,'FontSize',25,'TickLabelInterpreter','latex','XMinorGrid','on',...
    'YMinorTick','on');
if length(yFOM(:,1)) >=2
    plot(tFOM,yFOM(2,:),'b-+')
    plot(tPOD,yPOD(2,:),'y--o')
    plot(tNLMM,yNLMM(2,:),'r--')
  %  plot(tNLrk,yNLrk(2,:),'g-.')
    ylabel('$y_1(t)$, $y_2(t)$','Interpreter','LaTex')
end

end