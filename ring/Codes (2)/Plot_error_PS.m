function Plot_error_PS(tPOD,tNLMM,Error)

a=figure()
hold on
set(gca, 'YScale', 'log')
plot(tPOD, Error.errorPODout1,'color',[0.1215,0.145,0.46],'Linewidth',4)    %blue
plot(tNLMM, Error.errorNLMMOut1,'color',[1,0.45,0],'Linewidth',4)  %orange
%plot(tNLrk, Error.errorNLRKout1,'g')
xlabel('$t$','FontSize',25,'Interpreter','latex');
ylabel('$e(t)$','FontSize',25,'Interpreter','latex');
set(gca,'FontSize',25,'TickLabelInterpreter','latex','XMinorGrid','on',...
    'YMinorTick','on','YScale','log');
legend1=legend('POD','SO-NLMM-DEIM');
set(legend1,'Interpreter','latex','FontSize',25);
grid on
%saveas(a,'error2_r30.png')

end