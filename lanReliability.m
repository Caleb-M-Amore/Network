%% Lan Reliability Project

%% setup
clear; clc; close all;
% Define the range of K and probability p
K_val = [1, 5, 15, 50, 100];
p_val = 0:0.01:0.99;  % Probability of failure from 0 to 1 with increment of 0.01
N = 1000;  % Number of iterations for each simulation


%% Task 1
simulated = zeros(length(K_val), length(p_val));
calculated = zeros(length(K_val), length(p_val));
    
for i = 1:length(K_val)
    K = K_val(i);
    for j = 1:length(p_val)
        p = p_val(j);
        % Run the simulation
        simulated(i, j) = runSingleLinkSim(K, p, N);
        % Calculate the analytical result
        calculated(i, j) = K/(1 - p);  
    end
end
% Plotting the results
for i = 1:length(K_val)
    K = K_val(i);
    figure;  % Create a new figure for each K value
    semilogy(p_val, calculated(i, :), 'color', 'r', 'LineWidth', 2);
    hold on;
    semilogy(p_val, simulated(i, :), 'ko', 'color', 'b');
    title(sprintf('Task 1: Avg # of Transmissions for K = %d', K));
    xlabel('P(Failure)');
    ylabel('# of Transmissions');
    legend('Calculated', 'Simulated');
    grid; 
    hold off;
end

    figure;  
    semilogy(p_val, simulated(1, :), 'bo');
    hold on;
    semilogy(p_val, simulated(2, :), 'ro');
    semilogy(p_val, simulated(3, :), 'go');
    semilogy(p_val, simulated(4, :), 'co');
    semilogy(p_val, simulated(5, :), 'mo');
    semilogy(p_val, calculated(1, :), 'color', 'k');
    semilogy(p_val, calculated(2, :), 'color', 'k');
    semilogy(p_val, calculated(3, :),'color', 'k');
    semilogy(p_val, calculated(4, :), 'color', 'k');
    semilogy(p_val, calculated(5, :), 'color', 'k');
    legend('K=1', 'K=5', 'K=15', 'K=50', 'K=100');
    title('Task 1: Avg # of Transmissions for all K');
    xlabel('P(Failure)');
    ylabel('# of Transmissions');
    grid; 
    hold off;
%% Task 2
% Preallocate space for results
simulated = zeros(length(K_val), length(p_val));
calculated = zeros(length(K_val), length(p_val));

% Main simulation loop
for i = 1:length(K_val)
    K = K_val(i);
    for j = 1:length(p_val)
        p = p_val(j);
        % Run the simulation
        simulated(i, j) = runTwoSeriesLinkSim(K, p, N);
        % Calculate the analytical result
        % This is a placeholder for the analytical calculation
        % You will need to replace this with the actual formula
        calculated(i, j) = K/((1 - p)^2);  
    end
end
% Plotting the results
for i = 1:length(K_val)
    K = K_val(i);
    figure;  % Create a new figure for each K value
    semilogy(p_val, calculated(i, :), 'color', 'r', 'LineWidth', 2);
    hold on;
    semilogy(p_val, simulated(i, :), 'ko', 'color', 'b');
    title(sprintf('Task 2: Avg # of Transmissions for K = %d', K));
    xlabel('P(Failure)');
    ylabel('# of Transmissions');
    legend('Calculated', 'Simulated');
    grid; 
    hold off;
end

    figure;
    semilogy(p_val, simulated(1, :), 'ko', 'color', 'b');
    hold on;
    semilogy(p_val, simulated(2, :), 'ko', 'color', 'r');
    semilogy(p_val, simulated(3, :), 'ko', 'color', 'g');
    semilogy(p_val, simulated(4, :), 'ko', 'color', 'c');
    semilogy(p_val, simulated(5, :), 'ko', 'color', 'm');
    semilogy(p_val, calculated(1, :), 'color', 'k');
    semilogy(p_val, calculated(2, :), 'color', 'k');
    semilogy(p_val, calculated(3, :),'color', 'k');
    semilogy(p_val, calculated(4, :), 'color', 'k');
    semilogy(p_val, calculated(5, :), 'color', 'k');
    legend('K=1', 'K=5', 'K=15', 'K=50', 'K=100');
    title('Task 2: Avg # of Transmissions for all K');
    xlabel('P(Failure)');
    ylabel('# of Transmissions');
    grid; 
    hold off;
%% Task 3
% Preallocate space for results
    simulated = zeros(length(K_val), length(p_val));

% Main simulation loop
for i = 1:length(K_val)
    K = K_val(i);
    for j = 1:length(p_val)
        p = p_val(j);
        % Run the simulation
        simulated(i, j) = runTwoParallelLinkSim(K, p, N);
    end
end
% Plotting the results
for i = 1:length(K_val)
    K = K_val(i);
    figure;  % Create a new figure for each K value
    semilogy(p_val, simulated(i, :), 'ko', 'color', 'b');
    title(sprintf('Task 3: Avg # of Transmissions for K = %d', K));
    xlabel('P(Failure)');
    ylabel('# of Transmissions');
    grid; 
end

    
    figure;  % Create a new figure for each K value
    semilogy(p_val, simulated(1, :), 'ko', 'color', 'b');
    hold on;
    semilogy(p_val, simulated(2, :), 'ko', 'color', 'r');
    semilogy(p_val, simulated(3, :), 'ko', 'color', 'g');
    semilogy(p_val, simulated(4, :), 'ko', 'color', 'c');
    semilogy(p_val, simulated(5, :), 'ko', 'color', 'm');
    legend('K=1', 'K=5', 'K=15', 'K=50', 'K=100');
    title('Task 3: Avg # of Transmissions for all K');
    xlabel('P(Failure)');
    ylabel('# of Transmissions');
    grid; 
    hold off;

%% Task 4
% Preallocate space for results
    simulated = zeros(length(K_val), length(p_val));

% Main simulation loop
for i = 1:length(K_val)
    K = K_val(i);
    for j = 1:length(p_val)
        p = p_val(j);
        % Run the simulation
        simulated(i, j) = runCompoundNetworkSim(K, p, N);
    end
end
% Plotting the results
for i = 1:length(K_val)
    K = K_val(i);
    figure;  % Create a new figure for each K value
    semilogy(p_val, simulated(i, :), 'ko', 'color', 'b');
    title(sprintf('Task 4: Avg # of Transmissions for K = %d', K));
    xlabel('P(Failure)');
    ylabel('# of Transmissions');
    grid; 
end

    
    figure;  % Create a new figure for each K value
    semilogy(p_val, simulated(1, :), 'ko', 'color', 'b');
    hold on;
    semilogy(p_val, simulated(2, :), 'ko', 'color', 'r');
    semilogy(p_val, simulated(3, :), 'ko', 'color', 'g');
    semilogy(p_val, simulated(4, :), 'ko', 'color', 'c');
    semilogy(p_val, simulated(5, :), 'ko', 'color', 'm');
    legend('K=1', 'K=5', 'K=15', 'K=50', 'K=100');
    title('Task 4: Avg # of Transmissions for all K');
    xlabel('P(Failure)');
    ylabel('# of Transmissions');
    grid; 
    hold off;
 %% Task 5.1
% Define the range of K and probability p
    K_val = [1, 5, 10];
    p1=0.1;
    p2=0.6;
    p3_val = 0:0.01:0.99;  % Probability of failure from 0 to 1 with increment of 0.01
    N = 1000;  % Number of iterations for each simulation
% Preallocate space for results
    simulated = zeros(length(K_val), length(p3_val));
% Main simulation loop
for i = 1:length(K_val)
    K = K_val(i);
    for j = 1:length(p3_val)
        p3 = p3_val(j);
        % Run the simulation
        simulated(i, j) = runCustomCompoundNetworkSim(K, p1, p2, p3, N);
    end
end
% Plotting the results
    figure;  % Create a new figure for each K value
    semilogy(p3_val, simulated(1, :), 'ko', 'color', 'c');
    hold on;
    semilogy(p3_val, simulated(2, :), 'ko', 'color', 'r');
    semilogy(p3_val, simulated(3, :), 'ko', 'color', 'g');
    legend('K=1', 'K=5', 'K=10');
    title('Task 5.1: Avg # of Transmissions for all K');
    xlabel('P(Failure)');
    ylabel('# of Transmissions');
    grid; 
    hold off;
%% Task 5.2
% Define the range of probability p
    p1=0.6;
    p2=0.1;
    p3_val = 0:0.01:0.99;  % Probability of failure from 0 to 1 with increment of 0.01
    N = 1000;  % Number of iterations for each simulation
% Preallocate space for results
    simulated = zeros(length(K_val), length(p3_val));
% Main simulation loop
for i = 1:length(K_val)
    K = K_val(i);
    for j = 1:length(p3_val)
        p3 = p3_val(j);
        % Run the simulation
        simulated(i, j) = runCustomCompoundNetworkSim(K, p1, p2, p3, N);
    end
end
% Plotting the results
    figure;  % Create a new figure for each K value
    semilogy(p3_val, simulated(1, :), 'ko', 'color', 'c');
    hold on;
    semilogy(p3_val, simulated(2, :), 'ko', 'color', 'r');
    semilogy(p3_val, simulated(3, :), 'ko', 'color', 'g');
    legend('K=1', 'K=5', 'K=10');
    title('Task 5.2: Avg # of Transmissions for all K');
    xlabel('P(Failure)');
    ylabel('# of Transmissions');
    grid; 
    hold off;
%% Task 5.3
% Define the range of probability p
    p1=0.1;
    p3=0.6;
    p2_val = 0:0.01:0.99;  % Probability of failure from 0 to 1 with increment of 0.01
    N = 1000;  % Number of iterations for each simulation
% Preallocate space for results
    simulated = zeros(length(K_val), length(p2_val));
% Main simulation loop
for i = 1:length(K_val)
    K = K_val(i);
    for j = 1:length(p2_val)
        p2 = p2_val(j);
        % Run the simulation
        simulated(i, j) = runCustomCompoundNetworkSim(K, p1, p2, p3, N);
    end
end
% Plotting the results
    figure;  % Create a new figure for each K value
    semilogy(p2_val, simulated(1, :), 'ko', 'color', 'c');
    hold on;
    semilogy(p2_val, simulated(2, :), 'ko', 'color', 'r');
    semilogy(p2_val, simulated(3, :), 'ko', 'color', 'g');
    legend('K=1', 'K=5', 'K=10');
    title('Task 5.3: Avg # of Transmissions for all K');
    xlabel('P(Failure)');
    ylabel('# of Transmissions');
    grid; 
    hold off;
%% Task 5.4
% Define the range of probability p
    p1=0.6;
    p3=0.1;
    p2_val = 0:0.01:0.99;  % Probability of failure from 0 to 1 with increment of 0.01
    N = 1000;  % Number of iterations for each simulation
% Preallocate space for results
    simulated = zeros(length(K_val), length(p2_val));
% Main simulation loop
for i = 1:length(K_val)
    K = K_val(i);
    for j = 1:length(p2_val)
        p2 = p2_val(j);
        % Run the simulation
        simulated(i, j) = runCustomCompoundNetworkSim(K, p1, p2, p3, N);
    end
end
% Plotting the results
    figure;  % Create a new figure for each K value
    semilogy(p2_val, simulated(1, :), 'ko', 'color', 'c');
    hold on;
    semilogy(p2_val, simulated(2, :), 'ko', 'color', 'r');
    semilogy(p2_val, simulated(3, :), 'ko', 'color', 'g');
    legend('K=1', 'K=5', 'K=10');
    title('Task 5.4: Avg # of Transmissions for all K');
    xlabel('P(Failure)');
    ylabel('# of Transmissions');
    grid; 
    hold off;
%% Task 5.5
% Define the range of probability p
    p3=0.6;
    p2=0.1;
    p1_val = 0:0.01:0.99;  % Probability of failure from 0 to 1 with increment of 0.01
    N = 1000;  % Number of iterations for each simulation
% Preallocate space for results
    simulated = zeros(length(K_val), length(p1_val));
% Main simulation loop
for i = 1:length(K_val)
    K = K_val(i);
    for j = 1:length(p1_val)
        p1 = p1_val(j);
        % Run the simulation
        simulated(i, j) = runCustomCompoundNetworkSim(K, p1, p2, p3, N);
    end
end
% Plotting the results
    figure;  % Create a new figure for each K value
    semilogy(p1_val, simulated(1, :), 'ko', 'color', 'c');
    hold on;
    semilogy(p1_val, simulated(2, :), 'ko', 'color', 'r');
    semilogy(p1_val, simulated(3, :), 'ko', 'color', 'g');
    legend('K=1', 'K=5', 'K=10');
    title('Task 5.5: Avg # of Transmissions for all K');
    xlabel('P(Failure)');
    ylabel('# of Transmissions');
    grid; 
    hold off;
%% Task 5.6
% Define the range of probability p
    p3=0.1;
    p2=0.6;
    p1_val = 0:0.01:0.99;  % Probability of failure from 0 to 1 with increment of 0.01
    N = 1000;  % Number of iterations for each simulation
% Preallocate space for results
    simulated = zeros(length(K_val), length(p1_val));
% Main simulation loop
for i = 1:length(K_val)
    K = K_val(i);
    for j = 1:length(p1_val)
        p1 = p1_val(j);
        % Run the simulation
        simulated(i, j) = runCustomCompoundNetworkSim(K, p1, p2, p3, N);
    end
end
% Plotting the results
    figure;  % Create a new figure for each K value
    semilogy(p1_val, simulated(1, :), 'ko', 'color', 'c');
    hold on;
    semilogy(p1_val, simulated(2, :), 'ko', 'color', 'r');
    semilogy(p1_val, simulated(3, :), 'ko', 'color', 'g');
    legend('K=1', 'K=5', 'K=10');
    title('Task 5.6: Avg # of Transmissions for all K');
    xlabel('P(Failure)');
    ylabel('# of Transmissions');
    grid; 
    hold off;


%% Function runSingleLinkSim()
% Parameters
%  K - the number of packets in the application message
%  p - the probability of failure 
%  N - the number of simulations to run
%
% Returns: the average numeric result across the total simulations

function result = runSingleLinkSim(K,p,N)

    simResults = ones(1,N); % a place to store the result of each simulation
    
    for i=1:N
        txAttemptCount = 0; % transmission count
        pktSuccessCount = 0; % number of packets that have made it across
    
        while pktSuccessCount < K
            
            r = rand; % generate random number to determine if packet is successful (r > p)
            txAttemptCount = txAttemptCount + 1; % count 1st attempt
        
            % while packet transmissions is not successful (r < p)
            while r < p
                r = rand; % transmit again, generate new success check value r
                txAttemptCount = txAttemptCount + 1; % count additional attempt
            end
        
            pktSuccessCount = pktSuccessCount + 1; % increase success count after success (r > p)
        end
    
        simResults(i) = txAttemptCount; % record total number of attempted transmissions before entire application message (K successful packets) transmitted
    end

    result = mean(simResults);
end


%% Function runTwoSeriesLinkSim()
% Parameters
%  K - the number of packets in the application message
%  p - the probability of failure 
%  N - the number of simulations to run
%
% Returns: the average numeric result across the total simulations

function result = runTwoSeriesLinkSim(K,p,N)

    simResults = ones(1,N); % a place to store the result of each simulation
    
    for i=1:N
        txAttemptCount = 0; % transmission count
        pktSuccessCount = 0; % number of packets that have made it across
    
        while pktSuccessCount < K
            
            r1 = rand; % generate random number to determine if packet is successful (r1 > p)
            r2 = rand; % generate random number to determine if packet is successful (r2 > p)
            txAttemptCount = txAttemptCount + 1; % count 1st attempt
        
            % while packet transmissions is not successful (r1 < p)
            while r1 < p
                r1 = rand; % transmit again, generate new success check value r1
                txAttemptCount = txAttemptCount + 1; % count additional attempt
            end
            % while packet transmissions is not successful (r2 < p)
            while r2 < p
                r2 = rand; % transmit again, generate new success check value r2
                txAttemptCount = txAttemptCount + 1; % count additional attempt
            end
        
            pktSuccessCount = pktSuccessCount + 1; % increase success count after success (r > p)
        end
    
        simResults(i) = txAttemptCount; % record total number of attempted transmissions before entire application message (K successful packets) transmitted
    end

    result = mean(simResults);
end

%% Function runTwoParellelLinkSim()
% Parameters
%  K - the number of packets in the application message
%  p - the probability of failure 
%  N - the number of simulations to run
%
% Returns: the average numeric result across the total simulations

function result = runTwoParallelLinkSim(K,p,N)

    simResults = ones(1,N); % a place to store the result of each simulation
    
    for i=1:N
        txAttemptCount = 0; % transmission count
        pktSuccessCount = 0; % number of packets that have made it across
    
        while pktSuccessCount < K
            
            r1 = rand; % generate random number to determine if packet is successful (r1 > p)
            r2 = rand; % generate random number to determine if packet is successful (r2 > p)
            txAttemptCount = txAttemptCount + 1; % count 1st attempt
        
            % while packet transmissions is not successful (r1 < p) && (r2 <p)
            while ((r1 < p)&&(r2 < p))
                r1 = rand; % transmit again, generate new success check value r1
                r2 = rand;% transmit again, generate new success check value r2
                txAttemptCount = txAttemptCount + 1; % count additional attempt
            end
            pktSuccessCount = pktSuccessCount + 1; % increase success count after success (r > p)
        end
    
        simResults(i) = txAttemptCount; % record total number of attempted transmissions before entire application message (K successful packets) transmitted
    end

    result = mean(simResults);
end


%% Function runCompoundNetworkSim()
% Parameters
%  K - the number of packets in the application message
%  p - the probability of failure 
%  N - the number of simulations to run
%
% Returns: the average numeric result across the total simulations

function result = runCompoundNetworkSim(K,p,N)

    simResults = ones(1,N); % a place to store the result of each simulation
    
    for i=1:N
        txAttemptCount = 0; % transmission count
        pktSuccessCount = 0; % number of packets that have made it across
    
        while pktSuccessCount < K
            
            r1 = rand; % generate random number to determine if packet is successful (r1 > p)
            r2 = rand; % generate random number to determine if packet is successful (r2 > p)
            r3 = rand; % generate random number to determine if packet is successful (r2 > p)
            txAttemptCount = txAttemptCount + 1; % count 1st attempt
        
            % while packet transmissions is not successful (r1 < p)
            while ((r1 < p)&&(r2 < p))
                r1 = rand; % transmit again, generate new success check value r1
                r2 = rand;% transmit again, generate new success check value r2
                txAttemptCount = txAttemptCount + 1; % count additional attempt
            end
            % while packet transmissions is not successful (r3 < p)
            while r3 < p
                r3 = rand; % transmit again, generate new success check value r2
                txAttemptCount = txAttemptCount + 1; % count additional attempt
            end
        
            pktSuccessCount = pktSuccessCount + 1; % increase success count after success (r > p)
        end
    
        simResults(i) = txAttemptCount; % record total number of attempted transmissions before entire application message (K successful packets) transmitted
    end

    result = mean(simResults);
end


%% Function runCustomCompoundNetworkSim()
% Parameters
%  K - the number of packets in the application message
%  p - the probability of failure 
%  N - the number of simulations to run
%
% Returns: the average numeric result across the total simulations

function result = runCustomCompoundNetworkSim(K,p1, p2, p3,N)

    simResults = ones(1,N); % a place to store the result of each simulation
    
    for i=1:N
        txAttemptCount = 0; % transmission count
        pktSuccessCount = 0; % number of packets that have made it across
    
        while pktSuccessCount < K
            
            r1 = rand; % generate random number to determine if packet is successful (r1 > p)
            r2 = rand; % generate random number to determine if packet is successful (r2 > p)
            r3 = rand; % generate random number to determine if packet is successful (r2 > p)
            txAttemptCount = txAttemptCount + 1; % count 1st attempt
        
            % while packet transmissions is not successful (r1 < p)
            while ((r1 < p1)&&(r2 < p2))
                r1 = rand; % transmit again, generate new success check value r1
                r2 = rand;% transmit again, generate new success check value r2
                txAttemptCount = txAttemptCount + 1; % count additional attempt
            end
            % while packet transmissions is not successful (r3 < p)
            while r3 < p3
                r3 = rand; % transmit again, generate new success check value r2
                txAttemptCount = txAttemptCount + 1; % count additional attempt
            end
        
            pktSuccessCount = pktSuccessCount + 1; % increase success count after success (r > p)
        end
    
        simResults(i) = txAttemptCount; % record total number of attempted transmissions before entire application message (K successful packets) transmitted
    end

    result = mean(simResults);
end
