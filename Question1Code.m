%20.430 PSET 2 
%Authors: Laura Schwendeman


%Question 1
%% Question 1.1 

%calculate Fbound
kon = 0.05; %s
koff = 0.45; %s 

Fbound = kon/(kon+koff)

%simulate diffusion in 2D
Dfree = 1; %um^2/s

tau = 0.01; %s
nsteps = 1000; 
tspan = 0:tau:nsteps*tau; 


%initialize a vector with all the x,y positions
x = zeros(nsteps,1);
y = zeros(nsteps,1); 

%get the parameters of the random gaussian
meanStep = 0; 
sigmaStep = sqrt(2*Dfree*tau);

%get all the displacement steps for ea direction
xsteps = normrnd(meanStep, sigmaStep, nsteps-1, 1);
ysteps = normrnd(meanStep, sigmaStep, nsteps-1, 1);

%get all the possitions by cummulatively summing each step
x(2:end) = cumsum(xsteps); 
y(2:end) = cumsum(ysteps);

%plot the diffusion
figure; 
plot(x,y, 'o-'); 
xlabel('x (um)');
ylabel('y (um)');


%% Question 1.2

error_vec = zeros(1,nsteps);
meanerror_vec = zeros(1,nsteps);

for step = 2:length(error_vec)
    
    error_vec(step) = x(step)^2 + y(step)^2; 
    meanerror_vec(step) = error_vec(step); 


end

figure; 
plot(0:tau:(nsteps*tau-tau), meanerror_vec); 
xlabel("Time (s)")
ylabel("MSD um^2")


%% Question 1.3


%% Question 1.4
close all;

 BoundCounter = 0;

 tau = .01; %s
 kon = .05; %s^-1
 koff = .45; %s^-1
 
 nsteps = 10000; 


 Pon = 1-exp(-kon*tau); %%ask about this
 Poff = 1-exp(-koff*tau); 

 xloc = zeros(nsteps, 1);
 yloc = zeros(nsteps, 1);
 BoundLog = zeros(nsteps,1);

 %for diffusion timesteps
 Dfree = 1; %um^
 sigmaStep = sqrt(2*Dfree*tau);


 for i = 2:nsteps
    
     %generate a uniform random number
     rn = rand();

     %calculate Pon and Poff
     % stepsincelastflip = stepsincelastflip +1;
     % Pon = 1-exp(-kon*tau*stepsincelastflip);
     % Poff = 1-exp(-kon*tau*stepsincelastflip);

     if BoundCounter %bound
        
       xloc(i) = xloc(i-1);
       yloc(i) = yloc(i-1);

        if rn < Poff
            BoundCounter = 0;
            
        end
        
     else %free
        xstep = normrnd(0, sigmaStep);
        ystep = normrnd(0, sigmaStep);
       xloc(i) = xloc(i-1)+xstep;
       yloc(i) = yloc(i-1) + ystep;
      

        %check if crossing the wall
        [xloc(i), yloc(i)] = checkWall(xloc(i), yloc(i), xstep, ystep);


        if rn < Pon
            BoundCounter = 1;
            
        end

     end

     BoundLog(i) = BoundCounter; 
   
    
     %plot real time: 
      %plot
%  figure(1)
% plot(xloc(:), yloc(:),'-k')
% plot(xloc(BoundLog == 1), yloc(BoundLog == 1), 'ob');
% plot(xloc(BoundLog == 0), yloc(BoundLog == 0), 'or');
% hold on
% plot([4,4],[-2,2], 'r')
% plot([-4,-4],[-2,2], 'r')
% plot([4,-4],[2,2], 'r')
% plot([4,-4],[-2,-2], 'r')
% pause(0.01)


 end


 %plot
 figure(1)
 hold on;
plot(xloc(:), yloc(:),'-k')
plot(xloc(BoundLog == 0), yloc(BoundLog == 0), 'or');
plot(xloc(BoundLog == 1), yloc(BoundLog == 1), 'ob', 'LineWidth', 6);
hold on
plot([4,4],[-2,2], 'r')
plot([-4,-4],[-2,2], 'r')
plot([4,-4],[2,2], 'r')
plot([4,-4],[-2,-2], 'r')
xlabel('x (um)')
ylabel('y (um)')
legend("path","Free","Bound")

percentBound = sum(BoundLog)/nsteps

numBindingEvents = sum(diff(BoundLog) == -1)


%% 1.5
close all; 

nproteins = 500;
 nsteps = 1000; 

%diff-rxn parameters
tau = .01; %s
 kon = .05; %s^-1
 koff = .45; %s^-1


 %initialize protein array
 xloc = zeros(nsteps, nproteins);
 yloc = zeros(nsteps, nproteins);
 BoundLog = zeros(nsteps, nproteins);

 xloc(1,:) = rand(1,nproteins)*8-4;
 yloc(1,:) = rand(1,nproteins)*4-2;

 BoundLog(1,:) = binornd(1,Fbound,1,nproteins);


 %generate all of the random values
 xsteps = normrnd(0, sigmaStep,nsteps,nproteins);
 ysteps = normrnd(0, sigmaStep,nsteps,nproteins);

 %now run the simulation for each protein
 for i = 2:nsteps
    
    

  for p = 1:nproteins%loop through each protein
       %generate a uniform random number
        rn = rand();

     if BoundLog(i-1,p) %bound
        
       xloc(i,p) = xloc(i-1,p);
       yloc(i,p) = yloc(i-1,p);

        if rn < Poff
            BoundLog(i,p) = 0;
        else
            BoundLog(i,p) = 1;
        end
        
     else %free
        xstep = xsteps(i,p);
        ystep = ysteps(i,p);
       xloc(i,p) = xloc(i-1,p)+xstep;
       yloc(i,p) = yloc(i-1,p) + ystep;
      

        %check if crossing the wall
        [xloc(i,p), yloc(i,p)] = checkWall(xloc(i,p), yloc(i,p), xstep, ystep);


        if rn < Pon
            BoundLog(i,p) = 1;
        else
            BoundLog(i,p) = 0;
            
        end

     end

  end  
  
  if mod(i,10) == 0
      
         figure(2)
     hold off;
    
    plot(xloc(i, BoundLog(i,:) == 0), yloc(i, BoundLog(i,:) == 0), 'or',xloc(i, BoundLog(i,:) == 1), yloc(i, BoundLog(i,:) == 1), 'ob', 'LineWidth', 6);
    % plot([4,4],[-2,2], 'r')
    % plot([-4,-4],[-2,2], 'r')
    % plot([4,-4],[2,2], 'r')
    % plot([4,-4],[-2,-2], 'r')
    xlabel('x (um)')
    ylabel('y (um)')
    legend("Free","Bound")
    title('Question 1.5')
    
    pause(.01)
    

  end

  
 end



%% Question 1.5 -- Jacquin's verion

% variables
k_on = .05;
k_off = .45;
f_bound = k_on/(k_on+k_off);
D_free = 1;
mu = 0;
delta_tau = .01;
n = 1000;
num_proteins = 500;
sigma = sqrt(2*D_free*delta_tau);
time = [0:.01:n*delta_tau-delta_tau];
rand_vect_x = normrnd(mu, sigma, [num_proteins,n]);
rand_vect_y = normrnd(mu, sigma, [num_proteins,n]);
X = zeros(num_proteins,n);
Y = zeros(num_proteins,n);
X(:,1) = -4 + (4-(-4)).*rand(num_proteins,1);
Y(:,1) = -2 + (2-(-2)).*rand(num_proteins,1);
MSD = zeros(num_proteins,n);
fit = zeros(num_proteins,1);
rand_vect_state = rand(num_proteins,n);
state_mat = zeros(num_proteins,n);
state_mat(:,1)= rand(num_proteins,1) < .1;

% diffusion simulation
P_bound = 1 - exp(-k_on*delta_tau);
P_unbound = 1 - exp(-k_off*delta_tau);
for i = 1:num_proteins
   state = state_mat(i,1);
   for iter = 2:n
       if state == 0
           if rand_vect_state(i, iter) < P_bound
               Y(i,iter) = Y(i,iter-1);
               X(i,iter) = X(i,iter-1);
               state = 1;
           else
               Y(i,iter) = Y(i,iter-1)+rand_vect_y(i,iter);
               X(i,iter) = X(i,iter-1)+rand_vect_x(i,iter);
               if Y(i,iter) > 2
                   Y(i,iter) = Y(i,iter-1)-rand_vect_y(i,iter);
               elseif Y(i,iter) < -2
                   Y(i,iter) = Y(i,iter-1)-rand_vect_y(i,iter);
               end
               if X(i,iter) > 4
                   X(i,iter) = X(i,iter-1)-rand_vect_x(i,iter);
               elseif X(i,iter) < -4
                   X(i,iter) = X(i,iter-1)-rand_vect_x(i,iter);
               end
           end
       else
           if rand_vect_state(i,iter) < P_unbound
               Y(i,iter) = Y(i,iter-1)+rand_vect_y(i,iter);
               X(i,iter) = X(i,iter-1)+rand_vect_x(i,iter);
               if Y(i,iter) > 2
                   Y(i,iter) = Y(i,iter-1)-rand_vect_y(i,iter);
               elseif Y(i,iter) < -2
                   Y(i,iter) = Y(i,iter-1)-rand_vect_y(i,iter);
               end
               if X(i,iter) > 4
                   X(i,iter) = X(i,iter-1)-rand_vect_x(i,iter);
               elseif X(iter) < -4
                   X(i,iter) = X(i,iter-1)-rand_vect_x(i,iter);
               end
               state = 0;
           else
               Y(i,iter) = Y(i,iter-1);
               X(i,iter) = X(i,iter-1);
           end
       end
       state_mat(i,iter) = state;
       MSD(i,iter) = (X(i,iter)-X(i,1)).^2+(Y(i,iter)-Y(i,1)).^2; 
   end
end

%% plot movie
bound_mat = zeros(n, 1)
figure('Position', [100 300 1500 500])

counter = 1;

for t = 1:10:n
   clf
   subplot(1,2,1)
   hold on
   plot(X(:,t), Y(:,t), 'r.')
   plot(X(state_mat(:,t)==1, t), Y(state_mat(:,t)==1, t), 'b.')
   plot([4,4],[-2,2], 'k')
   plot([-4,-4],[-2,2], 'k')
   plot([4,-4],[2,2], 'k')
   plot([4,-4],[-2,-2], 'k')
   xlabel('X (um)')
   ylabel('Y (um)')
   title('XY Plane Position of 500 Simulated Protein in Nucleus at t')
   xlim([-4.5 4.5])
   ylim([-2.5 2.5])
   legend('State = free', 'State = bound')
   
   hold off
  
   subplot(1,2,2)
   ylim([0,1])
   ylabel('Fraction')
   xlabel('Free (0) or Bound (1)')
   title('Fraction of Free and Bound Proteins')
   bound_ratio = [1-(sum(state_mat(:,t))/num_proteins), sum(state_mat(:,t))/num_proteins]
   bound_mat(t,1) = sum(state_mat(:,t))/num_proteins;
   hold on
   bar([0,1], bound_ratio)
   xticks([0,1])
   hold off
   F(counter) = getframe(gcf); 
   counter = counter +1; 
   pause(0.001)

 end
avg_p_bound = mean(nonzeros(bound_mat))


%% make Movie file
writerObj = VideoWriter('Question1pt5.mp4', 'MPEG-4'); %Attempt to create an mp4 file
open(writerObj);
for t= 1:counter

    writeVideo(writerObj,F(t))

end

close(writerObj);

%% Question 1.6
  
%new information
  damageBounds = [-.5, .5; -4,4]; 

  kon_outside = .05; %1/s
  kon_inside = 4.5; %1/s
 Dfree = .1; %s^-1
 tau = .01; %s

   D_free_vals = [100, 10, 1, .1, .01, .001, .0001];
for Dfree = D_free_vals
%initialize the problem
%diff-rxn parameters


 koff = .45; %s^-1

 nsteps = round(30/tau); 

 %initialize protein array
 xloc = zeros(nsteps, nproteins);
 yloc = zeros(nsteps, nproteins);
 BoundLog = zeros(nsteps, nproteins);
 damageLocationState = zeros(nsteps, nproteins);

 xloc(1,:) = rand(1,nproteins)*8-4;
 yloc(1,:) = rand(1,nproteins)*4-2;

 BoundLog(1,:) = binornd(1,Fbound,1,nproteins);
 damageLocationState(1,:) = isInDamagedRegion(xloc(1,:), yloc(1,:),damageBounds);

 %generate all of the random values
 xsteps = normrnd(0, sigmaStep,nsteps,nproteins);
 ysteps = normrnd(0, sigmaStep,nsteps,nproteins);


 %to check diffusion's effects
 
 %now simulate
sigmaStep = sqrt(2*Dfree*tau);
   for i = 2:nsteps
     

  for p = 1:nproteins%loop through each protein

       %generate a uniform random number
        rn = rand();

       %determine if in the damaged region
       damaged = damageLocationState(i-1,p);

       %now calculate the probability of being on
       if damaged %if in the damaged region
           kon = kon_inside; 
       else
           kon = kon_outside;
       end

        Pon = 1-exp(-kon*tau);

     if BoundLog(i-1,p) %bound
        
       xloc(i,p) = xloc(i-1,p);
       yloc(i,p) = yloc(i-1,p);

        if rn < Poff
            BoundLog(i,p) = 0;
        else
            BoundLog(i,p) = 1;
        end
        
     else %free
        xstep = xsteps(i,p);
        ystep = ysteps(i,p);
       xloc(i,p) = xloc(i-1,p)+xstep;
       yloc(i,p) = yloc(i-1,p) + ystep;
      

        %check if crossing the wall
        [xloc(i,p), yloc(i,p)] = checkWall(xloc(i,p), yloc(i,p), xstep, ystep);


        if rn < Pon
            BoundLog(i,p) = 1;
        else
            BoundLog(i,p) = 0;
            
        end

     end

     %update whether in a damaged region now
     damageLocationState(i,p) = isInDamagedRegion(xloc(i,p), yloc(i,p),damageBounds);

  end  
  
  % if mod(i,10) == 0
  % 
  %        figure(3)
  %    hold off;
  %   subplot(1,2,1);
  %   plot(xloc(i, BoundLog(i,:) == 0), yloc(i, BoundLog(i,:) == 0), 'or',xloc(i, BoundLog(i,:) == 1), yloc(i, BoundLog(i,:) == 1), 'ob', 'LineWidth', 6);
  %   % plot([4,4],[-2,2], 'r')
  %   % plot([-4,-4],[-2,2], 'r')
  %   % plot([4,-4],[2,2], 'r')
  %   % plot([4,-4],[-2,-2], 'r')
  %   xlabel('x (um)')
  %   ylabel('y (um)')
  %   legend("Free","Bound")
  %   title('Question 1.6')
  % 
  %   subplot(1,2,2); 
  %   histogram(xloc(i
  % ,:),8/.5)
  %   xlabel('x (um)')
  %   ylabel('number of particles')
  %   title(['t = ' num2str(tau*i)]);
  % 
  % 
  %   pause(.01)
  % 

  % end

  
   end

   %plot Dfree comparisons
   figure(8);
   hold on; 
   plot((1:nsteps)*delta_tau, sum(xloc(:,:) < .5 & xloc(:,:) > -.5 & BoundLog(:,:) == 1, 2))

end

   %% 
   title("number of bound proteins in region over time")
    xlabel("time (s)"); 
    ylabel("number of proteins bound to dna")
    legend("D = 100", "D = 10", "D = 1", "D = .1", "D = .01", "D = .001", "D = .0001")
    set(gca, 'FontSize', 12)



 %% 1.6 Jacquin's version

 %variables
k_on = .05;
k_off = .45;
k_damage = 4.5;
f_bound = k_on/(k_on+k_off);
D_free = 1;
mu = 0;
delta_tau = .01;
n = 3000;
num_proteins = 1000;
sigma = sqrt(2*D_free*delta_tau);
time = [0:.01:n*delta_tau-delta_tau];
rand_vect_x = normrnd(mu, sigma, [num_proteins,n]);
rand_vect_y = normrnd(mu, sigma, [num_proteins,n]);
X = zeros(num_proteins,n);
Y = zeros(num_proteins,n);
X(:,1) = -4 + (4-(-4)).*rand(num_proteins,1);
Y(:,1) = -2 + (2-(-2)).*rand(num_proteins,1);
MSD = zeros(num_proteins,n);
fit = zeros(num_proteins,1);
rand_vect_state = rand(num_proteins,n);
state_mat = zeros(num_proteins,n);
state_mat(:,1)= rand(num_proteins,1) < f_bound;

% diffusion simulation
P_bound = 1 - exp(-k_on*delta_tau);
P_unbound = 1 - exp(-k_off*delta_tau);
P_damage = 1 - exp(-k_damage*delta_tau);
for i = 1:num_proteins
   state = state_mat(i,1);
   for iter = 2:n
       if state == 0
           if X(i, iter-1) >= -.5 && X(i, iter-1) <= .5
               P_on = P_damage;
           else
               P_on = P_bound;
           end
           if rand_vect_state(i, iter) < P_on
               Y(i,iter) = Y(i,iter-1);
               X(i,iter) = X(i,iter-1);
               state = 1;
           else
               Y(i,iter) = Y(i,iter-1)+rand_vect_y(i,iter);
               X(i,iter) = X(i,iter-1)+rand_vect_x(i,iter);
               if Y(i,iter) > 2
                   Y(i,iter) = Y(i,iter-1)-rand_vect_y(i,iter);
               elseif Y(i,iter) < -2
                   Y(i,iter) = Y(i,iter-1)-rand_vect_y(i,iter);
               end
               if X(i,iter) > 4
                   X(i,iter) = X(i,iter-1)-rand_vect_x(i,iter);
               elseif X(i,iter) < -4
                   X(i,iter) = X(i,iter-1)-rand_vect_x(i,iter);
               end
           end
       else
           if rand_vect_state(i,iter) < P_unbound
               Y(i,iter) = Y(i,iter-1)+rand_vect_y(i,iter);
               X(i,iter) = X(i,iter-1)+rand_vect_x(i,iter);
               if Y(i,iter) > 2
                   Y(i,iter) = Y(i,iter-1)-rand_vect_y(i,iter);
               elseif Y(i,iter) < -2
                   Y(i,iter) = Y(i,iter-1)-rand_vect_y(i,iter);
               end
               if X(i,iter) > 4
                   X(i,iter) = X(i,iter-1)-rand_vect_x(i,iter);
               elseif X(iter) < -4
                   X(i,iter) = X(i,iter-1)-rand_vect_x(i,iter);
               end
               state = 0;
           else
               Y(i,iter) = Y(i,iter-1);
               X(i,iter) = X(i,iter-1);
           end
       end
       state_mat(i,iter) = state;
       MSD(i,iter) = (X(i,iter)-X(i,1)).^2+(Y(i,iter)-Y(i,1)).^2; 
   end   
end

% plot movie
bound_mat = zeros(n, 1);
figure('Position', [100 300 1500 500])
counter = 1; 
for t = 1:10:n
   clf
   subplot(1,2,1)
   hold on
   plot(X(:,t), Y(:,t), 'r.')
   plot(X(state_mat(:,t)==1, t), Y(state_mat(:,t)==1, t), 'b.')
   plot([4,4],[-2,2], 'k')
   plot([-4,-4],[-2,2], 'k')
   plot([4,-4],[2,2], 'k')
   plot([4,-4],[-2,-2], 'k')
   xlabel('X (um)')
   ylabel('Y (um)')
   title('XY Plane Position of 1000 Simulated Protein in Nucleus at t')
   xlim([-4.5 4.5])
   ylim([-2.5 2.5])
   legend('State = free', 'State = bound')
   hold off
  
   subplot(1,2,2)
   histogram(X(:,t), 16)
   ylim([0,250])
   xlim([-4,4])
   ylabel('Number of Partices')
   xlabel('Free (0) or Bound (1)')
   title('Particle Density of Proteins Along X-Axis')
   F(counter) = getframe(gcf);
   pause(0.001)
 end
avg_p_bound = mean(nonzeros(bound_mat));

%% make Movie file
writerObj = VideoWriter('Question1pt6.mp4', 'MPEG-4'); %Attempt to create an mp4 file
open(writerObj);
for t= 1:counter

    writeVideo(writerObj,F(t))

end

close(writerObj);



 %% Functions 
 
 %set the new x and y positions if you cross the wall boundary
function [xnew,ynew] = checkWall(x,y,xstep, ystep)

    if y >2
        ynew = y - 2*ystep;
    elseif y <-2 
        ynew = y - 2*ystep; 
    else 
        ynew = y;
    end

    if x > 4
        xnew = x - 2*xstep;
    elseif x < -4
        xnew = x - 2*xstep; 
    else 
        xnew = x; 
    end

end


%to return if a state is damaged or not
function [damageBin] = isInDamagedRegion(x,y,bounds)

    damageBin = all((x < bounds(1,2) & x > bounds(1,1)) & (y < bounds(2,2) & y>bounds(2,1))); 

end
