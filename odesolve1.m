classdef odesolve1<handle
    properties
        a {mustBeNumeric}
        y {mustBeNumeric}
        v {mustBeNumeric}
        b {mustBeNumeric}
        eptol {mustBeNumeric}
        Tmin {mustBeNumeric}
        Tmax {mustBeNumeric}
        ansdatall = cell(5,1);
        timdatall = cell(5,1);
        probz = {}
        rmse = {}
        ampz = {}
        egvalall = {}
        egvecall = {}
        tlimall
        projsall
        count = 0 
        type 
        H1
        V1
        spidx
        tmp2
        tmp1
    end
    methods
        function projections(self)
            tim_max = zeros(1,self.count);
            for i =1:self.count
               tim_max(i)=  max(self.timdatall{i});
            end
            [~,max1rdx] = min(cell2mat(self.rmse));
            [~,max1tdx] = max(tim_max);
            if max1tdx == max1rdx
               [~,max1idx] = min(cell2mat(self.rmse));
            else
                [~,dummy] = sort(cell2mat(self.rmse));
                max1idx =dummy(3);
            end
            tim_data = self.timdatall{max1idx};
            projs = zeros(length(tim_data),length(self.y));            
            evol_data = self.ansdatall{max1idx};
            valall= zeros(length(tim_data),length(self.y));
            for i=1:length(tim_data)
                allHam = self.H1*tim_data(i)+self.V1;
                [eigvec,eigval]=eig(allHam);
                egval = diag(eigval).';
                valall(i,:) = egval;
                edvec1 = evol_data(i,:);
                edvec2 = edvec1/norm(edvec1);
                p1a = conj(eigvec')*edvec2';
                projs(i,:) = p1a';
            end
            self.projsall = projs;
        end
        function V11= spinx(self)
            N=length(self.y);
            sp = [(N-1)/2:-1:-(N-1)/2];
            smax=(N-1)/2;
            V=zeros(N,N);
            for i=1:N
                for j=1:N
                    if j==i+1
                       m1 =  sqrt(smax*(smax+1)-sp(j)*(sp(j)+1)); 
                       V(i,j) = m1/2;
                    elseif  i==j+1
                       m2 =  sqrt(smax*(smax+1)-sp(j)*(sp(j)-1)); 
                       V(i,j) = m2/2;
                    end
                end
            end
            V11 =V;
        end
        function energies_mk1(self)
            tsi = (self.Tmax-self.Tmin)/10000;
            tlim = [self.Tmin:tsi:self.Tmax];
            prodall = zeros(length(tlim),nchoosek(length(self.y),2));
            valall= zeros(length(tlim),length(self.y));
            impidx = zeros(length(tlim));
            for i=1:length(tlim)
                allHam = self.H1*tlim(i)+self.V1;
                [eigvec,eigval]=eig(allHam);
                egval = diag(eigval).';
                if sum(real(egval))==0 && sum(imag(egval))~=0 
                    impidx(i) = tlim(i);
                end
                valall(i,:) = egval;
                dotp1 = triu(eigvec'*eigvec,1);
                At = dotp1.';
                m11  = (1:size(At,1)).' >= (1:size(At,2));
                v11  = At(m11);
                v12 = v11(v11~=0);
                prodall(i,:)=v12';
            end
            impidx(impidx==0)=[];
            if isempty(impidx)
                impidx = 0;
            end
            self.tmp2 = max(impidx);
            self.tmp1 = min(impidx);
            clear impidx;
            self.tlimall = tlim;
            self.egvalall = valall;
            self.egvecall = prodall;
        end
        function dydt = diffeq(self,t1,y1)
            dydt = (-1i)*(self.H1*t1+self.V1)*y1;
        end
        function hammat(self)
            N = length(self.y);
            sp = [(N-1)/2:-1:-(N-1)/2];
            %sp(sp>0) = 1;
            %sp(sp<0) = -1;
            self.H1 = diag(sp)*self.a;
            if regexp(self.type,'^[nN][hH]$')
                V11 = spinx(self);
                self.V1= 1i*self.v*V11;
                %self.V1 =  diag(1i*self.v*ones(1,N-1),1) + ...
                %diag(1i*self.v*ones(1,N-1),-1);
            elseif regexp(self.type,'^[hH]$')
                %self.V1 =  diag(self.v*ones(1,N-1),1) + ...
                %diag(self.v*ones(1,N-1),-1);
                V11 = spinx(self);
                self.V1= self.v*V11;
            end
        end
        function rkhf_mk2(self)
            y1=self.y;
            t1 = self.Tmin;
            T = [self.Tmin self.Tmax];
            h = (T(2)-T(1))/1000;
            hstar = (T(2)-T(1))/(10*self.eptol);
            ansdat = zeros(round(hstar),length(y1));
            timdat = zeros(round(hstar),1);
            iter =0;
            while t1<=T(end)
                k1 = h * diffeq(self, t1, y1) ;
                k2 = h * diffeq(self, t1 + 0.25 * h, y1 + 0.25 * k1) ;
                k3 = h * diffeq(self, t1 + (3/8) * h, y1 + (3/32) * k1 + (9/32) * k2) ;
                k4 = h * diffeq(self, t1 + (12/13) * h, y1 + (1932/2197) * k1 + (-7200/2197) * k2 +...
                    (7296/2197) * k3) ;
                k5 = h * diffeq(self, t1 + 1 * h, y1 + (439/216) * k1 + (-8) * k2 +...
                    (3680/513) * k3-(845/4104) * k4) ;
                k6 = h * diffeq(self, t1 + (1/2) * h, y1 + (-8/27) * k1 + (2) * k2 +...
                    (-3544/2565) * k3+(1859/4104) * k4-(11/40) * k5) ;
                % Update next value of y 
                ynew1 = y1 + (25/216)*(k1 )+ (1408/2565) * k3 + (2197/4104) * k4 + ...
                (-1/5)*k5 ;
                ynew2 = y1 + (16/135)*(k1 )+ (6656/12825) * k3 + (28561/56430) * k4 + ...
                    (-9/50)*k5 +(2/55)*k6 ;
                %error
                ep1 = mean(abs(ynew1-ynew2));
                if ep1>=self.eptol
                    h = self.b*h*(self.eptol/ep1)^(0.2);
                else
                     h = self.b*h*(self.eptol/ep1)^(0.25);
                end
                % Update next value of x
                iter = iter+1;
                ansdat(iter,:)=ynew1';
                y1=ynew1;
                timdat(iter) = t1;
                t1 = t1 + h ;
            end
            ansdat = ansdat(1:iter,:);
            timdat = timdat(1:iter,:);
            self.count = self.count+1;
            self.timdatall{self.count} = timdat;
            self.ansdatall{self.count} = ansdat;
        end
        function runge_kutta_mk2(self)
            y1=self.y;
            t1 = self.Tmin;
            T = [self.Tmin self.Tmax];
            h = (T(2)-T(1))/1000;
            hstar = (T(2)-T(1))/self.eptol;
            ansdat = zeros(round(hstar),length(y1));
            timdat = zeros(round(hstar),1);         
            iter =0;
            while t1<=T(end)
                k1 = h * diffeq(self, t1, y1) ;
                k2 = h * diffeq(self, t1 + 0.5 * h, y1 + 0.5 * k1) ;
                k3 = h * diffeq(self, t1 + 0.5 * h, y1 + 0.5 * k2) ;
                k4 = h * diffeq(self,t1 + h, y1 + k3) ;
                % Update next value of y 
                y1 = y1 + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4) ;
                iter = iter +1;
                ansdat(iter, :) = y1';
                timdat(iter, :) = t1;
                % Update next value of t
                t1 = t1 + h ;
            end
            ansdat = ansdat(1:iter,:);
            timdat = timdat(1:iter,:);
            self.count = self.count+1;
            self.timdatall{self.count} = timdat;
            self.ansdatall{self.count} = ansdat;
        end  
        function rkck_mk2(self)
            y1=self.y;
            t1 = self.Tmin;
            T = [self.Tmin self.Tmax];
            h = (T(2)-T(1))/1000;
            hstar = (T(2)-T(1))/self.eptol;
            ansdat = zeros(round(hstar),length(y1));
            timdat = zeros(round(hstar),1);
            iter =0;
            while t1<=T(end)
                k1 = h * diffeq(self, t1, y1) ;
                k2 = h * diffeq(self, t1 + (1/5) * h, y1 + (1/5) * k1) ;
                k3 = h * diffeq(self, t1 + (3/10) * h, y1 + (3/40) * k1 + (9/40) * k2) ;
                k4 = h * diffeq(self, t1 + (3/5) * h, y1 + (3/10) * k1 + (-9/10) * k2 +...
                    (6/5) * k3) ;
                k5 = h * diffeq(self, t1 + 1 * h, y1 + (-11/54) * k1 + (5/2) * k2 +...
                    (-70/27) * k3+(35/27) * k4) ;
                k6 = h * diffeq(self, t1 + (7/8) * h, y1 + (1631/55296) * k1 + ...
                    (175/512) * k2 + (575/13824) * k3+(44275/110592) * k4+(253/4096) * k5) ;
           
                % Updat1e next1 value of y 
                ynew1 = y1 + (37/378)*(k1 )+ (250/621) * k3 + (125/594) * k4 + (512/1771)*k5 ;
                ynew2 = y1 + (2825/27648)*(k1 )+ (18575/48384) * k3 + (13525/55296) * k4 ... 
                +(277/14336)*k5 +(1/4)*k6 ;
                %error
                ep1 = mean(abs(ynew1-ynew2));
                if ep1>=self.eptol
                    h = self.b*h*(self.eptol/ep1)^(0.2);
                else
                     h = self.b*h*(self.eptol/ep1)^(0.25);
                end
                % Update next value of x 
                iter = iter+1;
                ansdat(iter,:)=ynew2';
                y1=ynew2;
                timdat(iter,:)=t1;
                t1 = t1 + h ;
            end
            ansdat = ansdat(1:iter,:);
            timdat = timdat(1:iter,:);
            self.count = self.count+1;
            self.timdatall{self.count} = timdat;
            self.ansdatall{self.count} = ansdat;
        end
        function rkdp_mk2(self)
            y1=self.y;
            t1 = self.Tmin;
            T = [self.Tmin self.Tmax];
            h = (T(2)-T(1))/1000;
            hstar = (T(2)-T(1))/self.eptol;
            ansdat = zeros(round(hstar),length(y1));
            timdat = zeros(round(hstar),1);
            iter =0;
            while t1<=T(end)
                k1 = h * diffeq(self, t1, y1) ;
                k2 = h * diffeq(self, t1 + (1/5) * h, y1 + (1/5) * k1) ;
                k3 = h * diffeq(self, t1 + (3/10) * h, y1 + (3/40) * k1 + (9/40) * k2) ;
                k4 = h * diffeq(self, t1 + (4/5) * h, y1 + (44/45) * k1 + (-56/15) * k2 +...
                    (32/9) * k3) ;
                k5 = h * diffeq(self, t1 + (8/9)* h, y1 + (19372/6561) * k1 + (-25360/2187) * k2 +...
                    (64448/6561) * k3+(-212/729) * k4) ;
                k6 = h * diffeq(self, t1 + 1 * h, y1 + (9017/3168) * k1 + ...
                    (-355/33) * k2 + (46732/5247) * k3+(49/176) * k4+(-5103/18656) * k5) ;
                k7 = h * diffeq(self, t1 + 1 * h, y1 + (35/384) * k1 + ...
                    (500/1113) * k3 + (125/192) * k4+(-2187/6784) * k5+(11/84) * k6) ; 
                % Updat1e next1 value of y 
                ynew1 = y1 + (35/384) * k1 + ...
                    (500/1113) * k3 + (125/192) * k4+(-2187/6784) * k5+(11/84) * k6;
                ynew2 = y1 + (5179/57600)*(k1 )+ (7571/16695) * k3 + (393/640) * k4 ... 
                +(-92097/339200)*k5 +(187/2100)*k6+(1/40)*k7 ;
                %error
                ep1 = mean(abs(ynew1-ynew2));
                if ep1>=self.eptol
                    h = self.b*h*(self.eptol/ep1)^(0.2);
                else
                     h = self.b*h*(self.eptol/ep1)^(0.25);
                end
                % Update next value of x 
                iter = iter+1;
                ansdat(iter,:)=ynew1';
                y1=ynew1;
                timdat(iter,:)=t1;
                t1 = t1 + h ;
            end
            ansdat = ansdat(1:iter,:);
            timdat = timdat(1:iter,:);
            self.count = self.count+1;
            self.timdatall{self.count} = timdat;
            self.ansdatall{self.count} = ansdat;
        end
        function execute(self)
                hammat(self)
                runge_kutta_mk2(self);
                rkhf_mk2(self);
                rkdp_mk2(self);
                %rkck_mk2(self);
                TSPAN = [self.Tmin,self.Tmax]; % Solve from t=1 to t=5
                IC = self.y; % y(t=0) = 1
                opts = odeset('RelTol',eps*1e+3,'AbsTol',eps);
                [timdat4, ansdat4] = ode45(@(t,y) diffeq(self,t,y), ...
                    TSPAN, IC,opts);
                self.count= self.count+1;
                self.timdatall{self.count} = timdat4;
                self.ansdatall{self.count} = ansdat4;
                probz_error(self)
        end
        function probz_error(self)
            for i=1:self.count
                self.ampz{i} = (abs(self.ansdatall{i}).^2);
                self.probz{i} = (sum(abs(self.ansdatall{i}).^2,2));
                self.rmse{i} = sqrt(mean((self.probz{i} -1).^2));
            end
        end
        function plot_data(self)
            figure;
            %C = {'k','b','r','g','y',[.5 .6 .7],[.2 .2 .6]};
            C = lines(self.count);
            strleg = {'Runge Kutta','RKhf','Rkdp','Ode45'};
            hold on;
            for i=1:self.count
                plot(self.timdatall{i},self.probz{i},'DisplayName',strleg{i},...
                    'LineWidth',1.2);
               % plot(self.timdatall{i},self.rmse{i},'--')
            end
            title('Sum of the probabilities');
            legend();
            xlabel('Time'); ylabel('Sum of the Probabilities');
            figure;
            hold on;
            linS = {'-','--',':','-.'};
            hp = zeros(2*self.count,1);
            l1= zeros(length(self.count),1);
            l2= zeros(length(self.y),1);
            for i=1:self.count
                for j=1:length(self.y)
                hp(i) = plot(self.timdatall{i},self.ampz{i}(:,j)./self.probz{i},'linestyle',linS{length(linS)-rem(j,length(linS))},'LineWidth',1.2,...
                'color',C(i,:));
                %hp(i+4) = plot(self.timdatall{i},self.ampz{i}(:,2)./self.probz{i},'-','LineWidth',1.2,...
                %   'color',C(i,:),'DisplayName', strleg{i});
               % plot(self.timdatall{i},self.rmse{i},'--')
                 if i==1
                     l2(j)=plot([NaN,NaN], 'linestyle',linS{length(linS)-rem(j,length(linS))},'color','k','DisplayName', ...
                         strcat('\psi_{',string(j),'}'));
                 end
                end
                l1(i) = plot([NaN,NaN], 'color', C(i,:),'DisplayName', strleg{i});
            end
            legend(horzcat(l2(1:end)',l1(1:end)));
            %legend(l1(1:end));
            %legend(l2(1:end))
            annotation('textbox', [0.25, 0.7, 0.1, 0.1], 'String', ...
               strcat('\Psi(-\infty)= ',mat2str(self.y')));%,'Interpreter','latex')
            title('Plot of probabilities as a function of time');
            %legend(hp(6:end));
            xlabel('Time'); ylabel('Probability');
        end
        function best_plot(self)
            tim_max = zeros(1,self.count);
            for i =1:self.count
               tim_max(i)=  max(self.timdatall{i});
            end
            [~,max1rdx] = min(cell2mat(self.rmse));
            [~,max1tdx] = max(tim_max);
            if max1tdx == max1rdx
               [~,max1idx] = min(cell2mat(self.rmse));
            else
                [~,dummy] = sort(cell2mat(self.rmse));
                max1idx =dummy(3);
            end
            self.spidx = max1idx;
            figure;
            linS = {'-','--',':','-.'};
            %C = {'k','b','r','g','y',[.5 .6 .7],[.2 .2 .6]};
            C = lines(length(self.y));
            strleg = {'Runge Kutta','RKhf','Rkdp','Ode45'};
            hold on;
            plot(self.timdatall{max1idx},self.probz{max1idx},'DisplayName',strleg{max1idx},...
                    'LineWidth',1.2,'LineStyle',linS{i});
            plot([self.tmp1 self.tmp1],[min(self.probz{max1idx}),max(self.probz{max1idx})],...
                'r--','DisplayName',strcat('t = '+string(self.tmp1)))
            plot([self.tmp2 self.tmp2],[min(self.probz{max1idx}),max(self.probz{max1idx})],...
                'r--','DisplayName',strcat('t = '+string(self.tmp2)))
            %txt = strcat('t =',string(self.tmp1));
            %text(self.tmp1+0.5,(min(self.probz{max1idx})+max(self.probz{max1idx}))/2,txt)
            %txt = strcat('t =',string(self.tmp2));
            %text(self.tmp2-0.5,(min(self.probz{max1idx})+max(self.probz{max1idx}))/2,txt)
            title('Sum of the probabilities');
            xlabel('Sum of the probabilities');
            ylabel('Time');
            figure;
            hold on;
            linS = {'-','--',':','-.'};
            hp = zeros(2*self.count,1);
            for j=1:length(self.y)
                hp(j) = plot(self.timdatall{max1idx},self.ampz{max1idx}(:,j)./self.probz{max1idx},'LineWidth',1.2,...
                'color',C(j,:),'LineStyle',linS{length(linS)-rem(j,length(linS))},'DisplayName', ...
                         strcat('\psi_{',string(j),'}'));
                %hp(j).Color(4) = 0.25;
            end
            plot([self.tmp1 self.tmp1],[0,1],'r--','DisplayName',strcat('t = '+string(self.tmp1)))
            plot([self.tmp2 self.tmp2],[0,1],'r--','DisplayName',strcat('t = '+string(self.tmp2)))
            %txt = strcat('t =',string(self.tmp1));
            %text(self.tmp1+0.5,0.5,txt)
            %txt = strcat('t =',string(self.tmp2));
            %text(self.tmp2-0.5,0.5,txt)
            legend();
            annotation('textbox', [0.25, 0.7, 0.1, 0.1], 'String', ...
               strcat('\Psi(-\infty)= ',mat2str(self.y')));%,'Interpreter','latex')
            title("Plot of probabilities as a function of time, a = "+self.a+",v = "+self.v);
            %legend(hp(6:end));
            xlabel('Time'); ylabel('Probability');
        end
        function plot_energies(self)
            %valall= sortrows(self.egvalall','ComparisonMethod','real')';
            valall = self.egvalall;
            vecall= self.egvecall;
            [~,s12]=size(valall);
            [~,s22]=size(vecall);
            linS = {'-','--',':','-.'};
            figure;
            hold on;
            hp = zeros(s12,1);
            hz = zeros(s12,1);
            C = lines(s12);
            for i=1:s12
                str1 =string("\Re(E_{"+i+"})");
                str2 =string("\Im(E_{"+i+"})");
                hp(i) = plot(self.tlimall,real(valall(:,i)),'LineWidth',1.2,...
                'color',C(i,:),'DisplayName', ...
                         str1);
                hz(i) = plot(self.tlimall,imag(valall(:,i)),'LineStyle','--','LineWidth',1.2,...
                'color',C(i,:),'DisplayName', ...
                         str2);     
            end
            plot([self.tmp1 self.tmp1],[min(min(real(valall))),max(max(real(valall)))],'r--',...
                'DisplayName',strcat('t = '+string(self.tmp1)))
            plot([self.tmp2 self.tmp2],[min(min(real(valall))),max(max(real(valall)))],'r--',...
                'DisplayName',strcat('t = '+string(self.tmp2)))
            %txt = strcat('t =',string(self.tmp1));
            %text(self.tmp1+0.5,(min(min(abs(valall)))+max(max(abs(valall))))/2,txt)
            %txt = strcat('t =',string(self.tmp2));
            %text(self.tmp2-0.5,(min(min(abs(valall)))+max(max(abs(valall))))/2,txt)
            title("a = "+self.a+",v = "+self.v);
            xlabel('Time')
            ylabel('EigenValues')
            legend()
            figure;
            C = lines(s22);
            hold on;
            ve1 = 1:1:length(self.y);
            ve2 = combvec(ve1,ve1);
            ve3 = fliplr(ve2');
            ve4 = unique(sort(ve3,2),'rows');
            ve4(ve4(:,1)==ve4(:,2),:)=[];
            hp1 = zeros(s22,1);
            for i=1:s22
                str1 =string("v"+ve4(i,1)+"-"+"v"+ve4(i,2));
                hp1(i) = plot(self.tlimall,abs(vecall(:,i)),'LineWidth',1.2,...
                'color',C(i,:),'LineStyle',linS{length(linS)-rem(i,length(linS))},'DisplayName', ...
                         str1);
            end
            plot([self.tmp1 self.tmp1],[min(min(abs(vecall))),max(max(abs(vecall)))],'r--',...
                'DisplayName',strcat('t = '+string(self.tmp1)))
            plot([self.tmp2 self.tmp2],[min(min(abs(vecall))),max(max(abs(vecall)))],'r--',...
                'DisplayName',strcat('t = '+string(self.tmp2)))
            %txt = strcat('t =',string(self.tmp1));
            %text(self.tmp1+0.5,(min(min(abs(vecall)))+max(max(abs(vecall))))/2,txt)
            %txt = strcat('t =',string(self.tmp2));
            %text(self.tmp2-0.5,(min(min(abs(vecall)))+max(max(abs(vecall))))/2,txt)
            title("a = "+self.a+",v = "+self.v);
            xlabel('Time')
            ylabel('Inner product of EigenStates')
            legend()
        end
        function plot_projs(self)
            tim_max = zeros(1,self.count);
            for i =1:self.count
               tim_max(i)=  max(self.timdatall{i});
            end
            [~,max1rdx] = min(cell2mat(self.rmse));
            [~,max1tdx] = max(tim_max);
            if max1tdx == max1rdx
               [~,max1idx] = min(cell2mat(self.rmse));
            else
                [~,dummy] = sort(cell2mat(self.rmse));
                max1idx =dummy(3);
            end
            tim_data = self.timdatall{max1idx};
            pnall=sqrt(sum(abs(self.projsall).^2,2));
            pall = self.projsall./pnall;
            [~,s12]=size(pall);
            C = lines(s12);
            hp1 = zeros(s12,1);
            figure;
            linS = {'-','--',':','-.'};
            hold on;
            for i=1:s12
                str1 =string("proj"+i);
                hp1(i) = plot(tim_data,abs(pall(:,i)).^2,'LineWidth',1.2,...
                'color',C(i,:),'LineStyle',linS{length(linS)-rem(i,length(linS))},'DisplayName', ...
                         str1);
            end
            plot([self.tmp1 self.tmp1],[min(min(abs(pall))),max(max(abs(pall)))],'r--'...
                ,'DisplayName',strcat('t = '+string(self.tmp1)))
            plot([self.tmp2 self.tmp2],[min(min(abs(pall))),max(max(abs(pall)))],'r--',...
                'DisplayName',strcat('t = '+string(self.tmp2)))
            title("a = "+self.a+",v = "+self.v);
            annotation('textbox', [0.25, 0.7, 0.1, 0.1], 'String', ...
               strcat('\Psi(-\infty)= ',mat2str(self.y')));
%             txt = strcat('t =',string(self.tmp1));
%             text(self.tmp1+0.5,(min(min(abs(pall)))+max(max(abs(pall))))/2,txt)
%             txt = strcat('t =',string(self.tmp2));
%             text(self.tmp2-0.5,(min(min(abs(pall)))+max(max(abs(pall))))/2,txt)
            xlabel('Time')
            ylabel('Projections on the EigenStates')
            legend()
        end
    end
end