//
//  project_examples.cpp
//  
//
//  Created by Céline Nguyen on 14/04/2020.
//

#include <stdio.h>
#include "monte_carlo.hpp"
#include <random>
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <set>

int main(){
    std::mt19937 G(time(NULL));
    {
        std::cout << std::endl << "Option Call Up-and-Out" << std::endl;
        // Parameters:
        int T = 1;
        double S0 = 100;
        double sigma = 0.3;
        double r = 0.05;
        double K = 100;
        int M = 50000; // Monte Carlo
        int N = 1000; // discretisation
        double H = 130;
        Up_and_Out up_and_out(S0, sigma, r, K, T, H, "call");
        
        // Euler scheme estimation:
        Stats MonteCarloD = up_and_out.MC_price(G,M,N, "discrete");
        Stats MonteCarloC = up_and_out.MC_price(G,M,N, "continuous");
        std::cout << "Monte Carlo discrete: " << MonteCarloD.get_MonteCarlo() << std::endl;
        std::cout << "Monte Carlo continuous: " << MonteCarloC.get_MonteCarlo() << std::endl;
        std::cout << "Closed formula : " << up_and_out.exact_price() << std::endl;
        // Comparison with discrete barrier option
        std::cout << "Price of option monitored discretely" << std::endl;
        H = H*exp(sigma*0.5826*sqrt(double(T)/double(N)));
        Up_and_Out up_and_out_discrete(S0, sigma, r, K, T, H, "call");
        std::cout << "Closed formula : " << up_and_out_discrete.exact_price() << std::endl;

        // Export data:
        Export("Up_and_Out_Call_Discrete.dat", MonteCarloD, up_and_out.exact_price(), up_and_out_discrete.exact_price());
        Export("Up_and_Out_Call_Continuous.dat", MonteCarloC, up_and_out.exact_price(), up_and_out_discrete.exact_price());
    
    }
    {
        std::cout << std::endl << "Option Put Up-and-Out" << std::endl;
        // Parameters:
        int T = 1;
        double S0 = 100;
        double sigma = 0.3;
        double r = 0.05;
        double K = 100;
        int M = 10000; // Monte Carlo
        int N = 1000; // discretisation
        double H = 130;
        Up_and_Out up_and_out(S0, sigma, r, K, T, H, "put");
        
        // Euler scheme estimation:
        Stats MonteCarloD = up_and_out.MC_price(G,M,N, "discrete");
        Stats MonteCarloC = up_and_out.MC_price(G,M,N, "continuous");
        std::cout << "Monte Carlo discrete: " << MonteCarloD.get_MonteCarlo() << std::endl;
        std::cout << "Monte Carlo continuous: " << MonteCarloC.get_MonteCarlo() << std::endl;
        std::cout << "Closed formula : " << up_and_out.exact_price() << std::endl;
        // Comparison with discrete barrier option
        std::cout << "Discrete version of option" << std::endl;
        H = H*exp(sigma*0.5826*sqrt(double(T)/double(N)));
        Up_and_Out up_and_out_discrete(S0, sigma, r, K, T, H, "put");
        std::cout << "Closed formula : " << up_and_out_discrete.exact_price() << std::endl;
        // Export data:
        Export("Up_and_Out_Put_Discrete.dat", MonteCarloD, up_and_out.exact_price(), up_and_out_discrete.exact_price());
        Export("Up_and_Out_Put_Continuous.dat", MonteCarloC, up_and_out.exact_price(), up_and_out_discrete.exact_price());
    }
    {
        std::cout << std::endl << "Option Call Up-and-In" << std::endl;
        // Parameters:
        int T = 1;
        double S0 = 100;
        double sigma = 0.3;
        double r = 0.05;
        double K = 100;
        int M = 10000; // Monte Carlo
        int N = 1000; // discretisation
        double H = 130;
        Up_and_In up_and_in(S0, sigma, r, K, T, H, "call");
        
        // Euler scheme estimation:
        Stats MonteCarloD = up_and_in.MC_price(G,M,N, "discrete");
        Stats MonteCarloC = up_and_in.MC_price(G,M,N, "continuous");
        std::cout << "Monte Carlo discrete: " << MonteCarloD.get_MonteCarlo() << std::endl;
        std::cout << "Monte Carlo continuous: " << MonteCarloC.get_MonteCarlo() << std::endl;
        std::cout << "Closed formula : " << up_and_in.exact_price() << std::endl;
        // Comparison with discrete barrier option
        std::cout << "Discrete version of option" << std::endl;
        H = H*exp(sigma*0.5826*sqrt(double(T)/double(N)));
        Up_and_In up_and_in_discrete(S0, sigma, r, K, T, H, "call");
        std::cout << "Closed formula : " << up_and_in_discrete.exact_price() << std::endl;
        // Export data:
        Export("Up_and_In_Call_Discrete.dat", MonteCarloD, up_and_in.exact_price(), up_and_in_discrete.exact_price());
        Export("Up_and_In_Call_Continuous.dat", MonteCarloC, up_and_in.exact_price(), up_and_in_discrete.exact_price());
    }
    {
        std::cout << std::endl << "Option Put Up-and-In" << std::endl;
        // Parameters:
        int T = 1;
        double S0 = 100;
        double sigma = 0.3;
        double r = 0.05;
        double K = 100;
        int M = 10000; // Monte Carlo
        int N = 1000; // discretisation
        double H = 130;
        Up_and_In up_and_in(S0, sigma, r, K, T, H, "put");
        
        // Euler scheme estimation:
        Stats MonteCarloD = up_and_in.MC_price(G,M,N, "discrete");
        Stats MonteCarloC = up_and_in.MC_price(G,M,N, "continuous");
        std::cout << "Monte Carlo discrete: " << MonteCarloD.get_MonteCarlo() << std::endl;
        std::cout << "Monte Carlo continuous: " << MonteCarloC.get_MonteCarlo() << std::endl;
        std::cout << "Closed formula : " << up_and_in.exact_price() << std::endl;
        // Comparison with discrete barrier option
        std::cout << "Discrete version of option" << std::endl;
        H = H*exp(sigma*0.5826*sqrt(double(T)/double(N)));
        Up_and_In up_and_in_discrete(S0, sigma, r, K, T, H, "put");
        std::cout << "Closed formula : " << up_and_in_discrete.exact_price() << std::endl;
        // Export data:
        Export("Up_and_In_Put_Discrete.dat", MonteCarloD, up_and_in.exact_price(), up_and_in_discrete.exact_price());
        Export("Up_and_In_Put_Continuous.dat", MonteCarloC, up_and_in.exact_price(), up_and_in_discrete.exact_price());
    }
    {
        std::cout << std::endl << "Option Call Down-and-Out" << std::endl;
        // Parameters:
        int T = 1;
        double S0 = 100;
        double sigma = 0.3;
        double r = 0.05;
        double K = 100;
        int M = 10000; // Monte Carlo
        int N = 1000; // discretisation
        double H = 90;
        Down_and_Out down_and_out(S0, sigma, r, K, T, H, "call");
        
        // Euler scheme estimation:
        Stats MonteCarloD = down_and_out.MC_price(G,M,N, "discrete");
        Stats MonteCarloC = down_and_out.MC_price(G,M,N, "continuous");
        std::cout << "Monte Carlo discrete: " << MonteCarloD.get_MonteCarlo() << std::endl;
        std::cout << "Monte Carlo continuous: " << MonteCarloC.get_MonteCarlo() << std::endl;
        std::cout << "Closed formula : " << down_and_out.exact_price() << std::endl;
        // Comparison with discrete barrier option
        std::cout << "Discrete version of option" << std::endl;
        H = H*exp(-sigma*0.5826*sqrt(double(T)/double(N)));
        Down_and_Out down_and_out_discrete(S0, sigma, r, K, T, H, "call");
        std::cout << "Closed formula : " << down_and_out_discrete.exact_price() << std::endl;
        // Export data:
        Export("Down_and_Out_Call_Discrete.dat", MonteCarloD, down_and_out.exact_price(), down_and_out_discrete.exact_price());
        Export("Down_and_Out_Call_Continuous.dat", MonteCarloC, down_and_out.exact_price(), down_and_out_discrete.exact_price());
    }
    {
        std::cout << std::endl << "Option Put Down-and-Out" << std::endl;
        // Parameters:
        int T = 1;
        double S0 = 100;
        double sigma = 0.3;
        double r = 0.05;
        double K = 100;
        int M = 10000; // Monte Carlo
        int N = 1000; // discretisation
        double H = 90;
        Down_and_Out down_and_out(S0, sigma, r, K, T, H, "put");
        
        // Euler scheme estimation:
        Stats MonteCarloD = down_and_out.MC_price(G,M,N, "discrete");
        Stats MonteCarloC = down_and_out.MC_price(G,M,N, "continuous");
        std::cout << "Monte Carlo discrete: " << MonteCarloD.get_MonteCarlo() << std::endl;
        std::cout << "Monte Carlo continuous: " << MonteCarloC.get_MonteCarlo() << std::endl;
        std::cout << "Closed formula : " << down_and_out.exact_price() << std::endl;
        // Comparison with discrete barrier option
        std::cout << "Discrete version of option" << std::endl;
        H = H*exp(-sigma*0.5826*sqrt(double(T)/double(N)));
        Down_and_Out down_and_out_discrete(S0, sigma, r, K, T, H, "put");
        std::cout << "Closed formula : " << down_and_out_discrete.exact_price() << std::endl;
        // Export data:
        Export("Down_and_Out_Put_Discrete.dat", MonteCarloD, down_and_out.exact_price(), down_and_out_discrete.exact_price());
        Export("Down_and_Out_Put_Continuous.dat", MonteCarloC, down_and_out.exact_price(), down_and_out_discrete.exact_price());
    }
    {
        std::cout << std::endl << "Option Call Down-and-In" << std::endl;
        // Parameters:
        int T = 1;
        double S0 = 100;
        double sigma = 0.3;
        double r = 0.05;
        double K = 100;
        int M = 10000; // Monte Carlo
        int N = 1000; // discretisation
        double H = 90;
        Down_and_In down_and_in(S0, sigma, r, K, T, H, "call");
        
        // Euler scheme estimation:
        Stats MonteCarloD = down_and_in.MC_price(G,M,N, "discrete");
        Stats MonteCarloC = down_and_in.MC_price(G,M,N, "continuous");
        std::cout << "Monte Carlo discrete: " << MonteCarloD.get_MonteCarlo() << std::endl;
        std::cout << "Monte Carlo continuous: " << MonteCarloC.get_MonteCarlo() << std::endl;
        std::cout << "Closed formula : " << down_and_in.exact_price() << std::endl;
        // Comparison with discrete barrier option
        std::cout << "Discrete version of option" << std::endl;
        H = H*exp(-sigma*0.5826*sqrt(double(T)/double(N)));
        Down_and_In down_and_in_discrete(S0, sigma, r, K, T, H, "call");
        std::cout << "Closed formula : " << down_and_in_discrete.exact_price() << std::endl;
        // Export data:
        Export("Down_and_In_Call_Discrete.dat", MonteCarloD, down_and_in.exact_price(), down_and_in_discrete.exact_price());
        Export("Down_and_In_Call_Continuous.dat", MonteCarloC, down_and_in.exact_price(), down_and_in_discrete.exact_price());
    }
    {
        std::cout << std::endl << "Option Put Down-and-In" << std::endl;
        // Parameters:
        int T = 1;
        double S0 = 100;
        double sigma = 0.3;
        double r = 0.05;
        double K = 100;
        int M = 50000; // Monte Carlo
        int N = 1000; // discretisation
        double H = 90;
        Down_and_In down_and_in(S0, sigma, r, K, T, H, "put");
        
        // Euler scheme estimation:
        Stats MonteCarloD = down_and_in.MC_price(G,M,N, "discrete");
        Stats MonteCarloC = down_and_in.MC_price(G,M,N, "continuous");
        std::cout << "Monte Carlo discrete: " << MonteCarloD.get_MonteCarlo() << std::endl;
        std::cout << "Monte Carlo continuous: " << MonteCarloC.get_MonteCarlo() << std::endl;
        std::cout << "Closed formula : " << down_and_in.exact_price() << std::endl;
        // Comparison with discrete barrier option
        std::cout << "Discrete version of option" << std::endl;
        H = H*exp(-sigma*0.5826*sqrt(double(T)/double(N)));
        Down_and_In down_and_in_discrete(S0, sigma, r, K, T, H, "put");
        std::cout << "Closed formula : " << down_and_in_discrete.exact_price() << std::endl;
        // Export data:
        Export("Down_and_In_Put_Discrete.dat", MonteCarloD, down_and_in.exact_price(), down_and_in_discrete.exact_price());
        Export("Down_and_In_Put_Continuous.dat", MonteCarloC, down_and_in.exact_price(), down_and_in_discrete.exact_price());
    }
    {
        std::cout << std::endl << "Error analysis (Up-and-Out Call)" << std::endl;
        // Parameters:
        int T = 1;
        double S0 = 100;
        double sigma = 0.3;
        double r = 0.05;
        double K = 100;
        int M = 50000; // Monte Carlo
        int N; // discretisation
        double H = 130;
        Up_and_Out up_and_out(S0, sigma, r, K, T, H, "call");
        std::cout << "Closed formula : " << up_and_out.exact_price() << std::endl;
        // Error convergence
        int N_start = 5;
        int N_stop = 205;
        int N_step = 10;
        Error error(up_and_out, N_start, N_stop, N_step, M);
        error.run(G);
        // Export data:
        Export("Error_Up_and_Out_Call.dat", error);
        std::cout << "Exported !" << std::endl;
    }
    

    return 0;
}

