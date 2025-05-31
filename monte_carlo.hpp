//
//  monte_carlo.hpp
//  
//
//  Created by Céline Nguyen on 14/04/2020.
//

#ifndef monte_carlo_hpp
#define monte_carlo_hpp

#include <stdio.h>
#include <random>
#include <algorithm>
#include <cmath>
#include <utility>
#include <string>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <unordered_set>
#include <iostream>

class Stats{
protected:
    std::vector<double> Values;
    std::vector<double> Mean;
    std::vector<double> Cumsum;
    std::vector<double> Var;
public:
    Stats(double x=0.): Values(), Mean(), Cumsum(), Var() {};
    std::vector<double> get_Values() const;
    std::vector<double> get_Mean() const;
    std::vector<double> get_Var() const;
    std::vector<double> get_uconfidence() const;
    std::vector<double> get_lconfidence() const;
    std::vector<double> get_error() const;
    std::vector<double> get_relative_err(double value) const;
    double get_MonteCarlo() const;
    friend Stats & operator+=(Stats &, double);
    friend std::ostream & operator << (std::ostream &, const Stats &);
    friend void Export(std::string, const Stats &, double, double);
};

template <class Statistique, class Measurement, class RNG>
void MonteCarlo(Statistique & res, const Measurement & f, RNG & G, long unsigned int n){
    for (long unsigned i=0;i<n; i++){
        res += f(G);
    }
};


class BarrierOption {
protected:
    double S0;
    double sigma;
    double r;
    double K;
    int T;
    double H;
    std::string option;
    double call_value;
    double put_value;
    double lambda;
    double y;
    double x1;
    double y1;
public:
    BarrierOption(double S0_, double sigma_, double r_, double K_, int T_, double H_, std::string option_);
    std::tuple<double, bool> Euler(std::mt19937 & G, int N, std::string type="discrete") const;
    double BrownianBridge(double dt, double z1, double z2) const;
    virtual double payoff(std::tuple<double, bool> asset) const = 0;
    virtual Stats MC_price(std::mt19937 & G, long unsigned int n, int N, std::string type="discrete") const;
    virtual double exact_price() const =0;
};


class Up_and_In : public BarrierOption {
public:
    Up_and_In(double S0_, double sigma_, double r_, double K_, int T_, double H_, std::string option_) : BarrierOption(S0_, sigma_, r_, K_, T_, H_, option_) {
        try { if (H_<S0) throw std::invalid_argument("Invalid barrier value for this barrier type.");}
        catch (const std::invalid_argument& msg) {
            std::cout << std::endl;
            std::cerr << msg.what() << std::endl;
            std::exit( EXIT_FAILURE );
        }
    };
    double payoff(std::tuple<double, bool> asset) const override;
    double exact_price() const override;
};

class Up_and_Out : public BarrierOption {
public:
    Up_and_Out(double S0_, double sigma_, double r_, double K_, int T_, double H_, std::string option_) : BarrierOption(S0_, sigma_, r_, K_, T_, H_, option_) {
        try { if (H_<S0) throw std::invalid_argument("Invalid barrier value for this barrier type.");}
        catch (const std::invalid_argument& msg) {
            std::cout << std::endl;
            std::cerr << msg.what() << std::endl;
            std::exit( EXIT_FAILURE );
        }
    };
    double payoff(std::tuple<double, bool> asset) const override;
    double exact_price() const override;
};
    
class Down_and_In : public BarrierOption {
public:
    Down_and_In(double S0_, double sigma_, double r_, double K_, int T_, double H_, std::string option_) : BarrierOption(S0_, sigma_, r_, K_, T_, H_, option_) {
        try { if (H_>S0) throw std::invalid_argument("Invalid barrier value for this barrier type.");}
        catch (const std::invalid_argument& msg) {
            std::cout << std::endl;
            std::cerr << msg.what() << std::endl;
            std::exit( EXIT_FAILURE );
        }
    };
    double payoff(std::tuple<double, bool> asset) const override;
    double exact_price() const override;
};

class Down_and_Out : public BarrierOption {
public:
    Down_and_Out(double S0_, double sigma_, double r_, double K_, int T_, double H_, std::string option_) : BarrierOption(S0_, sigma_, r_, K_, T_, H_, option_) {
            try { if (H_>S0) throw std::invalid_argument("Invalid barrier value for this barrier type.");}
            catch (const std::invalid_argument& msg) {
                std::cout << std::endl;
                std::cerr << msg.what() << std::endl;
                std::exit( EXIT_FAILURE );
            }
    };
    double payoff(std::tuple<double, bool> asset) const override;
    double exact_price() const override;
};
    

class Error {
protected:
    BarrierOption & option;
    long unsigned int M;
    std::vector<double> Discrete;
    std::vector<double> Continuous;
    std::vector<int> xAxes;
    int N_start;
    int N_stop;
    int N_step;
public:
    Error(BarrierOption & option_, int N_start_, int N_stop_, int N_step_, long unsigned int M_): option(option_), M(M_) {
        try {
            if (N_start_<1) throw std::invalid_argument("Invalid start value for number of discretization steps.");
            else N_start = N_start_;
            if (N_stop_<= N_start_) throw std::invalid_argument("Invalid stop value for number of discretization steps.");
            else N_stop = N_stop_;
            if (N_step_>(N_stop_-N_start_)) throw std::invalid_argument("Invalid step value for number of discretization steps.");
            else N_step = N_step_;
        }
        catch (const std::invalid_argument& msg) {
            std::cout << std::endl;
            std::cerr << msg.what() << std::endl;
            std::exit( EXIT_FAILURE );
        }
    }
    void run(std::mt19937 & G);
    friend void Export(std::string, const Error &);
};

#endif /* monte_carlo_hpp */
