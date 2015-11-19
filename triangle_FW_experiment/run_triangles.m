%% Script to reproduce the experiment of Appendix E in paper:
%
% On the Global Linear Convergence of Frank-Wolfe Optimization Variants, 
% Simon Lacoste-Julien and Martin Jaggi, NIPS 2015.
%
% It reproduces Figure 5 on p.23, and takes about 30 seconds to run on my
% laptop.

% suppose triangle with point at [-1,0],[0,0] and [cos(theta),sin(theta)]
% -> the smaller theta, the worse the width
% % width = sin(theta/2)
% diameter^2 = (1+cos(theta))^2+ sin(theta)^2 = 2*(1+cos(theta))
% or also: diameter = 2 cos(theta/2)
% so constant in rate is:
% width^2/diameter^2 = 1/4*(tan(theta/2))^2 
%
% Note: by looking instead at (<rhat, dhat>)^2 angle (tighter constant)
% we gain a factor of 4: it should just be (sin(theta/2))^2


% consider 1/2 ||x-x_0||^2 as objective with x_0 = [1/2,0];
% note: condition number of 1

SAVEPLOT = 0; % set to 1 to export the flot (it requires export_fig)

xstar = [-0.5; 0];

grad = @(x) (x-xstar);
f = @(x) 0.5*norm(x-xstar)^2;

clear opts

opts.maxK = 2000; 

opts.TOL = 1e-10;
opts.pureFW = 0;

%%

thetas = pi./[4,10,20,50,100,200,500,1000,1500,2000];

seeds = [5:347:(20*347)];

methods = {@AFW, @pFW};

for choice = 1:length(thetas)
    theta = thetas(choice); 

    % theoretical constant:
    constant = 1/4*(tan(theta/2))^2;
    %tighter version: 
    %constant = (sin(theta/2))^2; 
    
    V = [0,0; -1, 0; cos(theta), sin(theta)]'; % vertices for polytope

    for meth = 1:2
        for i_repeat = 1:length(seeds)
            s = seeds(i_repeat);
            opts.seed = s;
            solver = methods{meth};
            r = solver(V, grad, f, opts);
            if r.drop_finished
                rates(choice,meth,i_repeat) = nan;
                continue
            end
            fprintf('theta = %1.2g\n', theta)
            % estimate empirical linear rate from regression:
            y = log(r.primal(10:end));
            x = 1:length(r.primal(10:end));
            p = polyfit(x,y,1);
            rates(choice,meth,i_repeat) = -p(1);
        end
    end
    rates(choice,3,:) = constant; % theoretical rate stored in 3rd column
end

raw_rates = rates;

%%
for i = 1:size(rates,1)
    for j = 1:size(rates,2)
        I = find(~isnan(rates(i,j,:)));
        new_rates(i,j) = median(rates(i,j,I));
        U(i,j) = quantile(rates(i,j,I),0.75);
        L(i,j) = quantile(rates(i,j,I),0.25);
    end
end
rates = new_rates;

U = U-rates; % get difference between max and median...
L = rates - L;

%% PFW vs. theoretical constant:
figure
hd(2) = semilogx(thetas, rates(:,2)./rates(:,3),'xr-','MarkerSize',10);
hold on
eBar = errorbar(thetas, rates(:,2)./rates(:,3), L(:,2)./rates(:,3), U(:,2)./rates(:,3),'r','LineWidth',3);

%axis([thetas(end)/2, thetas(1)*2, 0,15])

set(gca, 'FontSize', 16)
xlabel('\theta (rad)','FontSize',20)
ylabel('ratio','FontSize',20)

%title('Ratio of empirical rate constant for PFW to theoretical constant')

% AFW vs. theoretical constant
%figure
hd(1) = semilogx(thetas, 4*rates(:,1)./rates(:,3),'sk-','MarkerSize',10);
%hold on
% theoretical constant for AFW is 1/4 of the one of PFW:
eBar = errorbar(thetas, 4*rates(:,1)./rates(:,3), 4*L(:,1)./rates(:,3), 4*U(:,1)./rates(:,3),'k','LineWidth',3);

axis([thetas(end)/2, thetas(1)*2, 0,60]);


legend(hd, 'AFW', 'PFW')
grid on

if SAVEPLOT
    export_fig('empirical_vs_theory','-pdf','-transparent');
else
    title('Ratio of empirical rate constant for AFW to theoretical constant')
end

%%
figure
% comparing PFW vs. AFW:
temp = squeeze(raw_rates(:,2,:)./ raw_rates(:,1,:));
for i = 1:size(temp,1)
    I = find(~isnan(temp(i,:)));
    PFW_vs_AFW(i) = median(temp(i,I));
    Ucomparison(i) = quantile(temp(i,I),0.75);
    Lcomparison(i) = quantile(temp(i,I),0.25);
end
Ucomparison = Ucomparison-PFW_vs_AFW;
Lcomparison = PFW_vs_AFW-Lcomparison;

semilogx(thetas, PFW_vs_AFW,'s-b','LineWidth', 3,'MarkerSize',10)
hold on
eBar = errorbar(thetas, PFW_vs_AFW, Lcomparison, Ucomparison,'b','LineWidth',3);
plot(thetas, ones(size(thetas)),'k--','LineWidth', 3)

set(gca, 'FontSize', 16)
xlabel('\theta (rad)','FontSize',20)
ylabel('ratio of PFW vs. AFW','FontSize',20)

grid on

if SAVEPLOT
    export_fig('ratio_PFW_vs_AFW','-pdf','-transparent');
else
    title('Ratio of rate for PFW over AFW')
end


