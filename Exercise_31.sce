clc; clear; funcprot(0);

function [S0, V0] = BS_EuPut_FiDi_Explicit (r, sigma, a, b, m, nu_max, T, K)
    q = 2*r/sigma^2;
    delta_t_tilde = T*sigma^2/(2*nu_max);
    delta_x_tilde = (b-a)/m;
    x_tilde = a + [0:m].*delta_x_tilde;
    lambda = delta_t_tilde/delta_x_tilde^2;
    w = max(exp(0.5.*x_tilde(1:(m-1))*(q-1))-exp(0.5.*x_tilde(1:(m-1))*(q+1)),0);

    // t_tilde loop
    for i=1:nu_max
        nu(1) = (1-2*lambda)*w(1)+lambda*w(2);
        nu(m-1) = lambda*w(m-2) + (1-2*lambda)*w(m-1);
        nu(2:(m-2)) = lambda.*w(1:(m-3)) + (1-2*lambda).*w(2:(m-2))+ lambda.*w(3:(m-1));
        w = nu;
    end
    
    // Computation of option values for t=0
    S0 = K.*exp(x_tilde(1:(m-1)));
    V0 = K.*w'.*exp(-0.5.*x_tilde(1:(m-1))*(q-1) - 0.5*sigma^2*T*((q-1)^2/4 + q));
endfunction

function V0 = BS_Price_Int (r, sigma, S0, T, K)
      
    // Specify integrand for integration formula.
    function y = integrand (x)
        y = 1/sqrt(2*%pi) * g( S0*exp((r-0.5*sigma^2)*T + sigma*sqrt(T)*x) ) * exp(-r*T) * exp(-x^2/2);
    endfunction
    // EuPut function
    function y = g(x)
        y = max(K - x, 0)
    endfunction
    V0 = intg(-10, 10, integrand);
endfunction
// Setting given variables
r = 0.05; sigma = 0.2; a = -0.7; b = 0.4; m = 100; nu_max = 2000; T = 1; K = 100;
[S0, V0] = BS_EuPut_FiDi_Explicit (r, sigma, a, b, m, nu_max, T, K)
scf(0)
clf()
plot(S0, V0)
for i=1:length(S0)
    V0(i) = BS_Price_Int (r, sigma, S0(i), T, K)
end
plot(S0, V0, "r")
// Labeling the axis and adding a legend for better overview
legend(["FiDi Explicit"; "Integration"])
ylabel("EuPut Value at t=0", "fontsize", 4, "color", "blue")
xlabel("Stock Price", "fontsize", 4, "color", "blue")
