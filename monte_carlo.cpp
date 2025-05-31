//
//  monte_carlo.cpp
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
#include <math.h>
#include <tuple>

// STATS
std::vector<double> Stats::get_Values() const {return Values;};

std::vector<double> Stats::get_Mean() const {return Mean;};

std::vector<double> Stats::get_Var() const {return Var;};

std::vector<double> Stats::get_uconfidence() const {
    std::vector<double> Uconfidence(Values.size());
    for (int i=0; i<Values.size(); i++){
        Uconfidence[i] = Mean[i]+(1.96*sqrt(Var[i]/(i+1)));
    }
    return Uconfidence;
};

std::vector<double> Stats::get_lconfidence() const {
    std::vector<double> Lconfidence(Values.size());
    for (int i=0; i<Values.size(); i++){
        Lconfidence[i] = Mean[i]-(1.96*sqrt(Var[i]/(i+1)));
    }
    return Lconfidence;
};

std::vector<double> Stats::get_error() const {
    std::vector<double> error(Values.size());
    for (int i=0; i<Values.size(); i++){
        error[i] = sqrt(Var[i]/(i+1));
    }
    return error;
}

std::vector<double> Stats::get_relative_err(double value) const {
    std::vector<double> error(Values.size());
    for (int i=0; i<Values.size(); i++){
        error[i] = abs((value-Mean[i])/value);
    }
    return error;
}


double Stats::get_MonteCarlo() const {return Mean.back();};

Stats & operator += (Stats & stat, double x){
    (stat.Values).push_back(x);
    int n = (stat.Values).size();
    if (n == 1){
        (stat.Mean).push_back(x);
        (stat.Cumsum).push_back(x*x);
        (stat.Var).push_back(0);
    }
    else{
        double newMean = ((stat.Mean).back()*(stat.Mean).size()+x)/n;
        (stat.Mean).push_back(newMean);
        double newCumsum = (stat.Cumsum).back()+x*x;
        (stat.Cumsum).push_back(newCumsum);
        double newVar = (stat.Cumsum).back()/n - (stat.Mean).back()*(stat.Mean).back();
        (stat.Var).push_back(newVar);
    }
    return stat;
}

std::ostream & operator << (std::ostream & flux, const Stats & stat){
    std::vector<double> uconf = stat.get_uconfidence();
    std::vector<double> lconf = stat.get_lconfidence();
    for (unsigned k=0; k<(stat.Mean).size(); k++){
        flux << k+1 << " " << (stat.Values)[k] << " " << (stat.Mean)[k] << " " << uconf[k] << " " << lconf[k] <<  std::endl;
    }
    return flux;
}

void Export(std::string s, const Stats & stat, double price, double price_discrete){
    std::ofstream fichier(s);
    std::vector<double> values = stat.get_Values();
    std::vector<double> mean = stat.get_Mean();
    std::vector<double> uconf = stat.get_uconfidence();
    std::vector<double> lconf = stat.get_lconfidence();
    std::vector<double> error1 = stat.get_error();
    std::vector<double> error2 = stat.get_relative_err(price);
    for (int i=0; i<values.size(); i++){
        fichier << i+1 << " " << values[i] << " " << mean[i] << " " << uconf[i] << " " << lconf[i] << " " << error1[i] << " " << error2[i] << " " << price << " " << price_discrete << std::endl;
    }
    fichier.close();
};


// BARRIER OPTION

double normalCDF(double x) {
    return std::erfc(-x/std::sqrt(2))/2;
}


BarrierOption::BarrierOption(double S0_, double sigma_, double r_, double K_, int T_, double H_, std::string option_): S0(S0_), sigma(sigma_), r(r_), K(K_), T(T_), H(H_), option(option_) {
    double d1 = (log(S0/K)+(r+sigma*sigma/2)*T)/(sigma*sqrt(T));
    double d2 = d1-sigma*sqrt(T);
    call_value = S0*normalCDF(d1)-K*exp(-r*T)*normalCDF(d2);
    put_value = K*exp(-r*T)*normalCDF(-d2)-S0*normalCDF(-d1);
    lambda = (r+sigma*sigma/2)/(sigma*sigma);
    y = (log(H*H/(S0*K)))/(sigma*sqrt(T))+lambda*sigma*sqrt(T);
    x1 = log(S0/H)/(sigma*sqrt(T))+lambda*sigma*sqrt(T);
    y1 = log(H/S0)/(sigma*sqrt(T))+lambda*sigma*sqrt(T);
};

double BarrierOption::BrownianBridge(double dt, double z1, double z2) const {
    double p = 1-exp(-2*(H-z1)*(H-z2)/(sigma*sigma*z1*z1*dt));
    if ((S0<H) && (z1<H) && (z2<H)) return p;
    else {
        if ((S0>H) && (z1>H) && (z2>H)) return p;
        else return 0;
    }
}

std::tuple<double, bool> BarrierOption::Euler(std::mt19937 & G, int N, std::string type) const{
    std::normal_distribution<double> Norm(0,1);
    bool Out_of_D = false;
    double dt = double(T)/double(N);
    double St = S0;
    for (int i=0; i<N; i++){
        double new_St = St + r*St*dt + sigma*St*Norm(G)*sqrt(dt);
        if (!Out_of_D) {
            if (((S0<=H) && (new_St>=H)) || ((S0>=H) && (new_St<=H))) Out_of_D = true;
            else {
                if (type=="continuous") {
                    double p = BrownianBridge(dt, St, new_St);
                    std::bernoulli_distribution Bern(p);
                    if (!Bern(G)) Out_of_D = true;
                }
            }
        }
        St = new_St;
    }
    return std::make_tuple(St, Out_of_D);
};


Stats BarrierOption::MC_price(std::mt19937 & G, long unsigned int n, int N, std::string type) const {
    auto Function_to_evaluate = [=](std::mt19937 & G){return exp(-r*T)*payoff(Euler(G, N, type));};
    Stats stat;
    MonteCarlo(stat, Function_to_evaluate, G, n);
    return stat;
};

/*
double BarrierOption::payoff(std::tuple<double, bool> asset) const {
    return 0;
};

double BarrierOption::exact_price() const {
    return 0;
};
*/

// UP-AND-IN

double Up_and_In::payoff(std::tuple<double, bool> asset) const {
    double ST = std::get<0>(asset);
    bool Out_of_D = std::get<1>(asset);
    double payoff;
    if ((option=="call") && (ST>K)) payoff = ST - K;
    else {
        if ((option=="put") && (ST<K)) payoff = K - ST;
        else payoff = 0;
    }
    if (Out_of_D) return payoff;
    else return 0;
};

double Up_and_In::exact_price() const{
    if (option=="call") {
        if (H>K) return S0*normalCDF(x1)-K*exp(-r*T)*normalCDF(x1-sigma*sqrt(T))-S0*pow(H/S0,2*lambda)*(normalCDF(-y)-normalCDF(-y1))+K*exp(-r*T)*pow(H/S0, 2*lambda-2)*(normalCDF(-y+sigma*sqrt(T))-normalCDF(-y1+sigma*sqrt(T)));
        else return call_value;
    }
    else {
        if (H>K) return -S0*pow(H/S0, 2*lambda)*normalCDF(-y)+K*exp(-r*T)*pow(H/S0, 2*lambda-2)*normalCDF(-y+sigma*sqrt(T));
        else return put_value+S0*normalCDF(-x1)-K*exp(-r*T)*normalCDF(-x1+sigma*sqrt(T))-S0*pow(H/S0, 2*lambda)*normalCDF(-y1)+K*exp(-r*T)*pow(H/S0, 2*lambda-2)*normalCDF(-y1+sigma*sqrt(T));
    }
};




// DOWN-AND-IN

double Down_and_In::payoff(std::tuple<double, bool> asset) const {
    double ST = std::get<0>(asset);
    bool Out_of_D = std::get<1>(asset);
    double payoff;
    if ((option=="call") && (ST>K)) payoff = ST - K;
    else {
        if ((option=="put") && (ST<K)) payoff = K - ST;
        else payoff = 0;
    }
    if (Out_of_D) return payoff;
    else return 0;
};
    
double Down_and_In::exact_price() const{
    if (option=="call") {
        if (H<=K) return S0*pow(H/S0, 2*lambda)*normalCDF(y)-K*exp(-r*T)*pow(H/S0, 2*lambda-2)*normalCDF(y-sigma*sqrt(T));
        else return call_value-S0*normalCDF(x1)+K*exp(-r*T)*normalCDF(x1-sigma*sqrt(T))+S0*pow(H/S0, 2*lambda)*normalCDF(y1)-K*exp(-r*T)*pow(H/S0, 2*lambda-2)*normalCDF(y1-sigma*sqrt(T));
    }
    else {
        if (H<K) return -S0*normalCDF(-x1)+K*exp(-r*T)*normalCDF(-x1+sigma*sqrt(T))+S0*pow(H/S0, 2*lambda)*(normalCDF(y)-normalCDF(y1))-K*exp(-r*T)*pow(H/S0, 2*lambda-2)*(normalCDF(y-sigma*sqrt(T))-normalCDF(y1-sigma*sqrt(T)));
        else return put_value;
    }
};


// UP-AND-OUT

double Up_and_Out::payoff(std::tuple<double, bool> asset) const {
    double ST = std::get<0>(asset);
    bool Out_of_D = std::get<1>(asset);
    double payoff;
    if ((option=="call") && (ST>K)) payoff = ST - K;
    else {
        if ((option=="put") && (ST<K)) payoff = K - ST;
        else payoff = 0;
    }
    if (!Out_of_D) return payoff;
    else return 0;
};

double Up_and_Out::exact_price() const{
    if (option=="call") {
        if (H>K) return call_value-S0*normalCDF(x1)+K*exp(-r*T)*normalCDF(x1-sigma*sqrt(T))+S0*pow(H/S0, 2*lambda)*(normalCDF(-y)-normalCDF(-y1))-K*exp(-r*T)*pow(H/S0, 2*lambda-2)*(normalCDF(-y+sigma*sqrt(T))-normalCDF(-y1+sigma*sqrt(T)));
        else return 0;
    }
    else {
        if (H>=K) return put_value+S0*pow(H/S0, 2*lambda)*normalCDF(-y)-K*exp(-r*T)*pow(H/S0, 2*lambda-2)*normalCDF(-y+sigma*sqrt(T));
        else return -S0*normalCDF(-x1)+K*exp(-r*T)*normalCDF(-x1+sigma*sqrt(T))+S0*pow(H/S0, 2*lambda)*normalCDF(-y1)-K*exp(-r*T)*pow(H/S0, 2*lambda-2)*normalCDF(-y1+sigma*sqrt(T));
    }
};


    
// DOWN-AND-OUT

double Down_and_Out::payoff(std::tuple<double, bool> asset) const {
    double ST = std::get<0>(asset);
    bool Out_of_D = std::get<1>(asset);
    double payoff;
    if ((option=="call") && (ST>K)) payoff = ST - K;
    else {
        if ((option=="put") && (ST<K)) payoff = K - ST;
        else payoff = 0;
    }
    if (!Out_of_D) return payoff;
    else return 0;
};

double Down_and_Out::exact_price() const{
    if (option=="call") {
        if (H<=K) return call_value-S0*pow(H/S0, 2*lambda)*normalCDF(y)+K*exp(-r*T)*pow(H/S0, 2*lambda-2)*normalCDF(y-sigma*sqrt(T));
        else return S0*normalCDF(x1)-K*exp(-r*T)*normalCDF(x1-sigma*sqrt(T))-S0*pow(H/S0, 2*lambda)*normalCDF(y1)+K*exp(-r*T)*pow(H/S0, 2*lambda-2)*normalCDF(y1-sigma*sqrt(T));
    }
    else {
        if (H<K) return put_value+S0*normalCDF(-x1)-K*exp(-r*T)*normalCDF(-x1+sigma*sqrt(T))-S0*pow(H/S0, 2*lambda)*(normalCDF(y)-normalCDF(y1))+K*exp(-r*T)*pow(H/S0, 2*lambda-2)*(normalCDF(y-sigma*sqrt(T))-normalCDF(y1-sigma*sqrt(T)));
        else return 0;
    }
};

// ERROR

void Error::run(std::mt19937 & G){
    for (int N=N_start; N<N_stop; N+=N_step){
        Stats MC_d = option.MC_price(G, M, N, "discrete");
        Stats MC_c = option.MC_price(G, M, N, "continuous");
        double price = option.exact_price();
        Discrete.push_back(MC_d.get_MonteCarlo()-price);
        Continuous.push_back(MC_c.get_MonteCarlo()-price);
        xAxes.push_back(N);
    }
};

void Export(std::string s , const Error & error){
    std::ofstream fichier(s);
    for (int i=0; i<error.Discrete.size(); i++){
        fichier << error.xAxes[i] << " " << error.Discrete[i] << " " << error.Continuous[i] << std::endl;
    }
    fichier.close();
};
