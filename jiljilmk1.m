clear
tic
import odesolve1.*
avals = [0.05,0.1,0.3,0.5,1,2,5];
yin = [1;0;0;0];
probsmk1 = zeros(length(avals),length (yin));
type1 = 'nh';
for i=1:length(avals)
    w = odesolve1;
    w.a = avals(i); w.b =0.85; w.eptol= 1e-6;w.Tmin = -50;w.Tmax = 50;
    w.v = 1;w.type =type1;w.y =yin;
    [~,e12] = max(yin);
    execute(w)
    energies_mk1(w)
    best_plot(w)
    plot_energies(w)
    projections(w)
    plot_projs(w)
    all_amps = w.ampz{w.spidx}./w.probz{w.spidx};
    probsmk1(i,:) = all_amps(end,:);
    FigName = char(strcat(string(type1)+string(length(yin))+...
        'levela'+string(i)+'prob'+'.png'));
    print(figure(2),[pwd+"\"+string(length(yin))+"Level"+upper(string(type1))+"new"+string(e12)+"\"+...
        string(FigName)],'-dpng','-r1200')
    %print(figure(2), [pwd+"\3LevelNHnew\"+string(FigName)]);
    FigName = strcat(string(type1)+string(length(yin))+...
        'levela'+string(i)+'energy'+'.png');
    print(figure(3),[pwd+"\"+string(length(yin))+"Level"+upper(string(type1))+"new"+string(e12)+"\"+...
        string(FigName)],'-dpng','-r1200')
    %saveas(figure(3), [pwd+"\3LevelNHnew\"+string(FigName)]);
    FigName = strcat(string(type1)+string(length(yin))+...
        'levela'+string(i)+'ip'+'.png');
    print(figure(4),[pwd+"\"+string(length(yin))+"Level"+upper(string(type1))+"new"+string(e12)+"\"+...
        string(FigName)],'-dpng','-r1200')
    %saveas(figure(4), [pwd+"\3LevelNHnew\"+string(FigName)]);
    FigName = strcat(string(type1)+string(length(yin))+...
        'levela'+string(i)+'proj'+'.png');
    print(figure(5),[pwd+"\"+string(length(yin))+"Level"+upper(string(type1))+"new"+string(e12)+"\"+...
        string(FigName)],'-dpng','-r1200')
    %saveas(figure(5), [pwd+"\3LevelNHnew\"+string(FigName)]);
    close all
    clear w
end
figure;
hold on
[~,s12] = size(probsmk1);
C = lines(s12);
mks = {'o','+','*','v','.','x','d','s'};
for i=1:s12
    plot(avals,probsmk1(:,i),'color',C(i,:),'LineWidth',1.2,...
    'DisplayName',strcat('\psi_{',string(i),'}'),'Marker',mks{i});
end
legend()
annotation('textbox', [0.25, 0.7, 0.1, 0.1], 'String', ...
               strcat('\Psi(-\infty)= ',mat2str(yin)));
xlabel('\alpha')
ylabel('Probabilties')
FigName = char('alphaplot.png');
print(figure(1),[pwd+"\"+string(length(yin))+"Level"+upper(string(type1))+"new"+string(e12)+"\"+...
        string(FigName)],'-dpng','-r1200')
% w=odesolve1;
% w.a = 5; w.b =0.85; w.eptol= 1e-6;w.Tmin = -50;w.Tmax = 50;
% w.v = 1;w.type ='nh';w.y =[0;1;0];
% execute(w)
% best_plot(w)
% projections(w)
% plot_projs(w)
% energies_mk1(w)
% plot_energies(w)
toc