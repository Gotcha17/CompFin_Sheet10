clc; clear; funcprot(0);

function V0 = UpOutPut_BS_MC_Richardson (S0, r, sigma, T, K, B, M, m)
    delta_t = T/(2*m); // Fine mash
    // generate a Matrix of mxM normaly distributed r.v. with mean = 0 and ...
    // sd = sqrt(delta_t)
    delta_W1 = grand(m, M, "nor", 0, sqrt(delta_t));
    delta_W2 = grand(m, M, "nor", 0, sqrt(delta_t));
    // Initialize matrices for no barrier hit information
    nobarrier_hit_fine = ones(1,M);
    nobarrier_hit_coarse = ones(1,M);
    // Set initial stock values for first time t=0
    Euler_fine = S0;  
    Euler_coarse = S0;
    // Since in the BS model dSt = r*St*dt + sigma*St*delta_W
    for i=1:m
        // Calculation matrices of values for fine grid and for coarse grid
        // Since grid of fine calculation is twice as fine, twice as many
        // calculation steps are needed
        Euler_fine = Euler_fine + r*Euler_fine*delta_t +...
                         sigma.*Euler_fine.*delta_W1(i,:);
        nobarrier_hit_fine = nobarrier_hit_fine.*(Euler_fine<B);
        
        Euler_fine = Euler_fine + r*Euler_fine*delta_t +...
                         sigma.*Euler_fine.*delta_W2(i,:);
        nobarrier_hit_fine = nobarrier_hit_fine.*(Euler_fine<B);
        
        // For the coarse grid, stepsize is twice as large: 2*delta_t
        // and also both brownian motions are needed for this stepsize
        Euler_coarse = Euler_coarse + r*Euler_coarse*2*delta_t +...
                         sigma.*Euler_coarse.*(delta_W1(i,:)+delta_W2(i,:));
        nobarrier_hit_coarse = nobarrier_hit_coarse.*(Euler_coarse<B);
    end
    // Define put option function g
    function y = g(x)
        y = max((K - x), 0);
    endfunction
    // 
    V_fine = nobarrier_hit_fine.*(g(Euler_fine)*exp(-r*T));
    V_coarse = nobarrier_hit_coarse.*(g(Euler_coarse)*exp(-r*T));
    // Calculating final result for option price at time t=0 by averaging over...
    // all simulated prices with the MC simulation
    V0 = mean(2*V_fine-V_coarse); 
endfunction

// Function from Exercise 12
function V_0 = UpOutPut_BinMod(S_0, r, sigma, T, K, B, M)
    delta_t = T/M;  //calculation of delta_t
    Beta = (exp(-r*delta_t)+exp((r+sigma^2)*delta_t))/2;  //calculation of beta for CRR
    u=Beta+sqrt((Beta^2)-1); //so u>d is true
    d=u^-1; //because ud=1
    q=(exp(r*delta_t)-d)/(u-d); //calculation of succes probability (u)
    S=zeros(M+1,M+1); //creation stock price matrix
    S(1,1)=S_0; //Setting stock price at t=0 as initial price in the stock matrix
    
    for i=2:M+1 //Initializing algo for computation of stock price 
        for j=1:i
            S(j,i)=S(1,1)*u^(j-1)*d^(i-j); //with j upwards and i-j downwards movements
        end
    end
    
    V=-ones(M+1,M+1); //creating option value matrix
    V(:,M+1)=max((K-S(:,M+1)), 0); //calculation of option values for last column
    
    for i=M:-1:1 //Initializing algo for computation of option price
        //option value is zero if at current point stock price is higher than the barrier
        V(1:i,i)=exp(-r*delta_t)*(q*V(2:i+1,i+1)+(1-q)*V(1:i,i+1)).*(S(1:i,i) < B);
    end
    V_0 = V(1,1); //setting of first element of the option value matrix as option price at time t=0
endfunction

// Set values for the input variables
S0 = 100; S_0 = S0; r = 0.05; sigma = 0.2; T = 1; K = 100; B = 110; M = 10000; m = 250;
// Run functions and display result
V0 = UpOutPut_BS_MC_Richardson (S0, r, sigma, T, K, B, M, m)
M = 1000;
V_0 = UpOutPut_BinMod(S_0, r, sigma, T, K, B, M)
disp("Price of UpOutPut BS_MC_Richardson: "+string(V0))
disp("Price of UpOutPut BinMod: "+string(V_0))
