%--------------------------------------------------------------------------
clear all
clc

t0 = clock;
%Input location
citys=[37.2101283	126.9748391
37.2087386	126.9720948
37.208337	126.9707966
37.2092428	126.9694662
37.2102169	126.9688869
37.2111397	126.96805
37.2112252	126.9669986
37.2104903	126.9671702
37.2098409	126.9681358
37.2089864	126.9680071
37.2078756	126.9680285
37.2082516	126.9692945
37.2063204	126.9658828
37.2073629	126.9679588
37.2078671	126.9697559
37.2077902	126.9710112
37.2080209	126.972127
37.2070382	126.9728458
37.2068246	126.972127
37.2071835	126.9711936
37.2071236	126.9701743
37.207004	126.9691229
37.206816	126.9681358
37.2064913	126.9673848
37.2057906	126.9681358
37.2062435	126.968801
37.2063204	126.9700348
37.2063033	126.9710433
37.2055343	126.971097
37.2055599	126.9700992
37.2054659	126.9690049
37.2051326	126.9680178
37.2044917	126.9666553
37.2045174	126.9686615
37.2036628	126.9684362
37.2025006	126.967814
37.2009197	126.9675779
37.2017486	126.9691765
37.202663	126.9699597
37.202663	126.9709468
37.2013641	126.969949
37.2005437	126.9702387
37.1993302	126.9702172
37.2018085	126.9715369
37.2027143	126.9688439
37.2037312	126.9697559
37.2047994	126.9694769
37.2042097	126.9708073
37.2035688	126.9708717
37.2053634	126.9721913
37.2057223	126.973393
37.2057223	126.9750345
37.2055684	126.9762683
37.2057137	126.9775021
37.2059017	126.9784355
37.2051839	126.9782531
37.2037056	126.9783926
37.2032356	126.979841
37.2028339	126.9792724
37.2038081	126.9772768
37.2045857	126.977309
37.2032954	126.9765902
37.2019025	126.9773519
37.2019195	126.9761717
37.2023212	126.976676
37.204637	126.9749165
37.2040388	126.9740045
37.2035005	126.9737041
37.2039363	126.972127
37.2045943	126.9727707
37.2045857	126.9735324
37.2051412	126.974498
37.2044319	126.976558
37.2067049	126.9777703
37.2065426	126.9783819
37.2075936	126.9737685
37.2078585	126.9728565
37.2088497	126.9687045
37.2095675	126.9686401
37.2106527	126.9678998
37.2102938	126.9666553
37.2115071	126.9664299
37.2117378	126.9673419
37.2065939	126.9725668
];
%--------------------------------------------------------------------------
%% calculate the distance 
n = size(citys,1);
D = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;     
        end
    end    
end
%--------------------------------------------------------------------------

m = 75;                              
alpha = 1;                           
beta = 5;                           
vol = 0.2;                          
Q = 10;                            
Heu_F = 1./D;                      
Tau = ones(n,n);                  
Table = zeros(m,n);               
iter = 1;                        
iter_max = 100;                    
Route_best = zeros(iter_max,n);        
Length_best = zeros(iter_max,1);    
Length_ave = zeros(iter_max,1);     
Limit_iter = 0;                     
%-------------------------------------------------------------------------

while iter <= iter_max
    
      start = zeros(m,1);
      for i = 1:m
          temp = randperm(n);
          start(i) = temp(1);
      end
      Table(:,1) = start; 
      
      citys_index = 1:n;
     
      for i = 1:m
       
         for j = 2:n
             has_visited = Table(i,1:(j - 1));        
             allow_index = ~ismember(citys_index,has_visited);   
             allow = citys_index(allow_index);  
             P = allow;
            
             for k = 1:length(allow)
                 P(k) = Tau(has_visited(end),allow(k))^alpha * Heu_F(has_visited(end),allow(k))^beta;
             end
             P = P/sum(P);
            
            Pc = cumsum(P);     
            target_index = find(Pc >= rand);
            target = allow(target_index(1));
            Table(i,j) = target;
         end
      end
     
      Length = zeros(m,1);
      for i = 1:m
          Route = Table(i,:);
          for j = 1:(n - 1)
              Length(i) = Length(i) + D(Route(j),Route(j + 1));
          end
          Length(i) = Length(i) + D(Route(n),Route(1));
      end
     
      if iter == 1
          [min_Length,min_index] = min(Length);
          Length_best(iter) = min_Length;  
          Length_ave(iter) = mean(Length);
          Route_best(iter,:) = Table(min_index,:);
          Limit_iter = 1; 
          
      else
          [min_Length,min_index] = min(Length);
          Length_best(iter) = min(Length_best(iter - 1),min_Length);
          Length_ave(iter) = mean(Length);
          if Length_best(iter) == min_Length
              Route_best(iter,:) = Table(min_index,:);
              Limit_iter = iter; 
          else
              Route_best(iter,:) = Route_best((iter-1),:);
          end
      end
    
      Delta_Tau = zeros(n,n);
      
      for i = 1:m
         
          for j = 1:(n - 1)
              Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q/Length(i);
          end
          Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1)) + Q/Length(i);
      end
      Tau = (1-vol) * Tau + Delta_Tau;
  
    iter = iter + 1;
    Table = zeros(m,n);
end
%--------------------------------------------------------------------------

[Shortest_Length,index] = min(Length_best);
Shortest_Route = Route_best(index,:);
Time_Cost=etime(clock,t0);
disp(['shortest distance:' num2str(Shortest_Length)]);
disp(['shortest path:' num2str([Shortest_Route Shortest_Route(1)])]);
disp(['iteration:' num2str(Limit_iter)]);
disp(['run time:' num2str(Time_Cost) 's']);
order=[Shortest_Route Shortest_Route(1)]
path=[];
for i=order
    citys(i,:);
    path=[path;citys(i,:)];
end
path
%--------------------------------------------------------------------------

figure(1)
plot([citys(Shortest_Route,1);citys(Shortest_Route(1),1)],...  
     [citys(Shortest_Route,2);citys(Shortest_Route(1),2)],'o-');
grid on
for i = 1:size(citys,1)
    text(citys(i,1),citys(i,2),['   ' num2str(i)]);
end
text(citys(Shortest_Route(1),1),citys(Shortest_Route(1),2),'       start');
text(citys(Shortest_Route(end),1),citys(Shortest_Route(end),2),'       end');
xlabel('x')
ylabel('y')
title(['shortest path(shortest path:' num2str(Shortest_Length) ')'])
figure(2)
plot(1:iter_max,Length_best,'b')
legend('shortest distance')
xlabel('iteration')
ylabel('distance')
title('Convergence trajectory of algorithm')
%--------------------------------------------------------------------------