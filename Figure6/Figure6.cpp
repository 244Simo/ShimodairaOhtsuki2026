# include <stdio.h>
# include <iostream>
# include <random>
# include <cmath>
# include <fstream>
# include <vector>
using namespace std;

////Function prototype
double get_w_resident();
double get_w_mutant();
double get_z_resident();
////Function
double get_w_resident(double x, double y, double v, double beta, double gamma, double q){
    double w_r;
    if(q ==1){
        w_r=log((1-y)*(1-v)*x*v)+(beta+gamma*y*(1-v))*v*(1-x);
    }
    else{
        w_r = (log(x*v*(1-y)*(1-v))+beta*(1-x)*v)+(-log(1-q*(1-exp(-gamma*y*(1-v)*(1-x)*v))));
    }
    return w_r;
}
double get_w_mutant(double x_r, double y_r, double v_r, double x_m, double y_m, double v_m, double beta, double gamma, double q){
    double z_r;
    double z_m;
    double w_m;
    z_r = (log(x_r*v_r)+beta*(1-x_r)*v_r)+(-log(1-q*(1-exp(-gamma*y_r*(1-v_r)*(1-x_r)*v_r))));
    if(q ==1){
        w_m = log((1-y_m)*(1-v_m)*x_m*v_m)+(beta+gamma*y_m*(1-v_m))*v_m*(1-x_m);
    }
    else{
        w_m = log((1-y_m)*(1-v_m)*(((1-q)*(1-exp(-beta*v_m*(1-x_m)))*(x_r*v_r*exp(beta*(1-x_r)*v_r)*1/(1-q*(1-exp(-gamma*y_r*(1-v_r)*(1-x_r)*v_r))))+x_m*v_m)/(1-q*(1-exp(-(beta+gamma*y_m*(1-v_m))*v_m*(1-x_m))))));
    }
    return w_m;
}

double get_z_resident(double x, double y, double v, double beta, double gamma, double q){
    double z_r;
    z_r = (log(x*v)+beta*(1-x)*v)*(-log(1-q*(1-exp(-gamma*y*(1-v)*(1-x)*v))));
    return z_r;
}

//// main program
int main(void){
    int beta_count; //counter for loop
    int gamma_count;//counter for loop
    int beta_max;   //max value of β
    int gamma_max;  //max value of γ
    double delta_beta;  // stepsize for calculation on β
    double delta_gamma; // stepsize for calculation on γ

    double beta;
    double gamma;
    double q_value;
    int max_step;           //number of perturbation
    double mutation_step;   //standard deviation of mutations
    double mut_rate;        //mutation rate of (v,x,y)

    int pxm;    //determine whether mutation happens on x(if pxm=1 -> mutation/ pxm=0 -> no-mutation)
    int pym;    //determine whether mutation happens on y(if pym=1 -> mutation/ pym=0 -> no-mutation)
    int pvm;    //determine whether mutation happens on v(if pvm=1 -> mutation/ pvm=0 -> no-mutation)

    ///variants for calculation
    //initial values for numerical calculation
    double x_i;
    double y_i;
    double v_i;
    //resident values during calculation
    double x_r;
    double y_r;
    double v_r;
    double w_r;
    //mutant values during calculation
    double x_m;
    double y_m;
    double v_m;
    double w_m;

    ///parameters
    beta_max=30;
    gamma_max=30;
    delta_beta=0.1;
    delta_gamma=0.1;

    max_step=1*pow(10,4);
    mutation_step=0.001;
    mut_rate=0.5;

    //// generate random number
    // if you want to set seed value mamually
    //int seed = 0;
    //mt19937 gen(seed);

    // if you want to change seed value randomly
    random_device rd;
    mt19937 gen(rd());

    /// array for record the result, (β,γ,(β,γ,x,y,v))
    // suppose two sets of initial values and compare the results
    // broad teaching ESS: start from inner points and cannot reach terminal ESS(as far as we explored)
    // only teaching ESS: start from terminal points and can reach terminal ESS
    //q=1.0
    //broad teaching ESS
    vector <vector <vector <double> > > end_values1_10(beta_max/delta_beta+1, vector < vector <double> > (gamma_max/delta_gamma+1, vector <double>(5)));
    //only teaching ESS
    vector <vector <vector <double> > > end_values2_10(beta_max/delta_beta+1, vector < vector <double> > (gamma_max/delta_gamma+1, vector <double>(5)));
    //q=0.95
    //broad teaching ESS
    vector <vector <vector <double> > > end_values1_095(beta_max/delta_beta+1, vector < vector <double> > (gamma_max/delta_gamma+1, vector <double>(5)));
    //only teaching ESS
    vector <vector <vector <double> > > end_values2_095(beta_max/delta_beta+1, vector < vector <double> > (gamma_max/delta_gamma+1, vector <double>(5)));
    //q=0.9
    //broad teaching ESS
    vector <vector <vector <double> > > end_values1_09(beta_max/delta_beta+1, vector < vector <double> > (gamma_max/delta_gamma+1, vector <double>(5)));
    //only teaching ESS
    vector <vector <vector <double> > > end_values2_09(beta_max/delta_beta+1, vector < vector <double> > (gamma_max/delta_gamma+1, vector <double>(5)));
    //q=0.8
    //broad teaching ESS
    vector <vector <vector <double> > > end_values1_08(beta_max/delta_beta+1, vector < vector <double> > (gamma_max/delta_gamma+1, vector <double>(5)));
    //only teaching ESS
    vector <vector <vector <double> > > end_values2_08(beta_max/delta_beta+1, vector < vector <double> > (gamma_max/delta_gamma+1, vector <double>(5)));

    //// q=1.0
    q_value=1;
    ///broad teaching ESS
    //explore the parameter space
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        beta=beta_count*delta_beta;
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            gamma=gamma_count*delta_gamma;
            //initial values
            x_i=0.25;
            y_i=0.25;
            v_i=0.25;
            //set initial values to resident values
            x_r=x_i;
            y_r=y_i;
            v_r=v_i;
            ///numerical calculation(generate perturbation)
            for(int i=0; i<max_step; i++){
                bernoulli_distribution dist_pxm(mut_rate);
                bernoulli_distribution dist_pym(mut_rate);
                bernoulli_distribution dist_pvm(mut_rate);
                pxm = dist_pxm(gen);
                pym = dist_pym(gen);
                pvm = dist_pvm(gen);                    
                //slight mutation(normal distribution)
                normal_distribution<> dist_mut_x(x_r, mutation_step);
                normal_distribution<> dist_mut_y(y_r, mutation_step);
                normal_distribution<> dist_mut_v(v_r, mutation_step);
                //if mutated value fall within the domain, adopt mutated value
                //if not, retain the original value(mutation fails)
                if (pxm ==1){
                    x_m=dist_mut_x(gen);
                    if(x_m <= 0 or x_m >= 1){
                        x_m=x_r;
                    }
                }
                else if(pxm ==0){
                    x_m=x_r;
                }
                if (pym ==1){
                    y_m=dist_mut_y(gen);
                    if(y_m <= 0 or y_m >= 1){
                        y_m=y_r;
                    }
                }
                else if(pym ==0){
                    y_m=y_r;
                }
                if (pvm ==1){
                    v_m=dist_mut_v(gen);
                    if(v_m <= 0 or v_m >= 1){
                        v_m=v_r;
                    }
                }
                else if(pvm ==0){
                    v_m=v_r;
                }

                // comapare fitness between resident one and mutant one(judge the success of invasion)
                w_r=get_w_resident(x_r,y_r,v_r,beta,gamma,q_value);
                w_m=get_w_mutant(x_r,y_r,v_r,x_m,y_m,v_m,beta,gamma,q_value);
                // if fitness of mutant>fitness of resident, invasion success and update the resident values
                if(w_r<w_m){
                    x_r=x_m;
                    y_r=y_m;
                    v_r=v_m;
                }
            }

            // record results in arrays
            end_values1_10[beta_count][gamma_count][0]=beta;
            end_values1_10[beta_count][gamma_count][1]=gamma;
            end_values1_10[beta_count][gamma_count][2]=x_r;
            end_values1_10[beta_count][gamma_count][3]=y_r;
            end_values1_10[beta_count][gamma_count][4]=v_r;
            // reset the initial values
            x_i=0;
            y_i=0;
            v_i=0;
        }
    }

    ///only teaching ESS
    //explore the parameter space
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        beta=beta_count*delta_beta;
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            gamma=gamma_count*delta_gamma;
            //initial values
            x_i=0.99;
            y_i=0.05;
            v_i=0.25;
            //set initial values to resident values
            x_r=x_i;
            y_r=y_i;
            v_r=v_i;
            //numerical calculation(generate perturbation)
            for(int i=0; i<max_step; i++){
                bernoulli_distribution dist_pxm(mut_rate);
                bernoulli_distribution dist_pym(mut_rate);
                bernoulli_distribution dist_pvm(mut_rate);
                pxm = dist_pxm(gen);
                pym = dist_pym(gen);
                pvm = dist_pvm(gen);                    
                //slight mutation(normal distribution)
                normal_distribution<> dist_mut_x(x_r, mutation_step);
                normal_distribution<> dist_mut_y(y_r, mutation_step);
                normal_distribution<> dist_mut_v(v_r, mutation_step);
                //if mutated value fall within the domain, adopt mutated value
                //if not, retain the original value(mutation fails)
                if (pxm ==1){
                    x_m=dist_mut_x(gen);
                    if(x_m <= 0 or x_m >= 1){
                        x_m=x_r;
                    }
                }
                else if(pxm ==0){
                    x_m=x_r;
                }
                if (pym ==1){
                    y_m=dist_mut_y(gen);
                    if(y_m <= 0 or y_m >= 1){
                        y_m=y_r;
                    }
                }
                else if(pym ==0){
                    y_m=y_r;
                }     
                if (pvm ==1){
                    v_m=dist_mut_v(gen);
                    if(v_m <= 0 or v_m >= 1){
                        v_m=v_r;
                    }
                }
                else if(pvm ==0){
                    v_m=v_r;
                }

                // comapare fitness between resident one and mutant one(judge the success of invasion)
                w_r=get_w_resident(x_r,y_r,v_r,beta,gamma,q_value);
                w_m=get_w_mutant(x_r,y_r,v_r,x_m,y_m,v_m,beta,gamma,q_value);
                // if fitness of mutant>fitness of resident, invasion success and update the resident values
                if(w_r<w_m){
                    x_r=x_m;
                    y_r=y_m;
                    v_r=v_m;
                }
            }
            // record results in arrays
            end_values2_10[beta_count][gamma_count][0]=beta;
            end_values2_10[beta_count][gamma_count][1]=gamma;
            end_values2_10[beta_count][gamma_count][2]=x_r;
            end_values2_10[beta_count][gamma_count][3]=y_r;
            end_values2_10[beta_count][gamma_count][4]=v_r;
            // reset the initial values
            x_i=0;
            y_i=0;
            v_i=0;
        }
    }

    //// q=0.95
    q_value=0.95;
    ///broad teaching ESS
    //explore the parameter space
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        beta=beta_count*delta_beta;
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            gamma=gamma_count*delta_gamma;
            //initial values
            x_i=0.25;
            y_i=0.25;
            v_i=0.25;
            //set initial values to resident values
            x_r=x_i;
            y_r=y_i;
            v_r=v_i;
            ///numerical calculation(generate perturbation)
            for(int i=0; i<max_step; i++){
                bernoulli_distribution dist_pxm(mut_rate);
                bernoulli_distribution dist_pym(mut_rate);
                bernoulli_distribution dist_pvm(mut_rate);
                pxm = dist_pxm(gen);
                pym = dist_pym(gen);
                pvm = dist_pvm(gen);                    
                //slight mutation(normal distribution)
                normal_distribution<> dist_mut_x(x_r, mutation_step);
                normal_distribution<> dist_mut_y(y_r, mutation_step);
                normal_distribution<> dist_mut_v(v_r, mutation_step);
                //if mutated value fall within the domain, adopt mutated value
                //if not, retain the original value(mutation fails)
                if (pxm ==1){
                    x_m=dist_mut_x(gen);
                    if(x_m <= 0 or x_m >= 1){
                        x_m=x_r;
                    }
                }
                else if(pxm ==0){
                    x_m=x_r;
                }
                if (pym ==1){
                    y_m=dist_mut_y(gen);
                    if(y_m <= 0 or y_m >= 1){
                        y_m=y_r;
                    }
                }
                else if(pym ==0){
                    y_m=y_r;
                }                    
                if (pvm ==1){
                    v_m=dist_mut_v(gen);
                    if(v_m <= 0 or v_m >= 1){
                        v_m=v_r;
                    }
                }
                else if(pvm ==0){
                    v_m=v_r;
                }

                // comapare fitness between resident one and mutant one(judge the success of invasion)
                w_r=get_w_resident(x_r,y_r,v_r,beta,gamma,q_value);
                w_m=get_w_mutant(x_r,y_r,v_r,x_m,y_m,v_m,beta,gamma,q_value);
                // if fitness of mutant>fitness of resident, invasion success and update the resident values
                if(w_r<w_m){
                    x_r=x_m;
                    y_r=y_m;
                    v_r=v_m;
                }
            }
            // record results in arrays
            end_values1_095[beta_count][gamma_count][0]=beta;
            end_values1_095[beta_count][gamma_count][1]=gamma;
            end_values1_095[beta_count][gamma_count][2]=x_r;
            end_values1_095[beta_count][gamma_count][3]=y_r;
            end_values1_095[beta_count][gamma_count][4]=v_r;
            // reset the initial values
            x_i=0;
            y_i=0;
            v_i=0;
        }
    }

    ///only teaching ESS
    //explore the parameter space
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        beta=beta_count*delta_beta;
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            gamma=gamma_count*delta_gamma;
            //initial values
            x_i=0.99;
            y_i=0.05;
            v_i=0.25;
            //set initial values to resident values
            x_r=x_i;
            y_r=y_i;
            v_r=v_i;
            ///numerical calculation(generate perturbation)
            for(int i=0; i<max_step; i++){
                bernoulli_distribution dist_pxm(mut_rate);
                bernoulli_distribution dist_pym(mut_rate);
                bernoulli_distribution dist_pvm(mut_rate);
                pxm = dist_pxm(gen);
                pym = dist_pym(gen);
                pvm = dist_pvm(gen);                    
                //slight mutation(normal distribution)
                normal_distribution<> dist_mut_x(x_r, mutation_step);
                normal_distribution<> dist_mut_y(y_r, mutation_step);
                normal_distribution<> dist_mut_v(v_r, mutation_step);
                //if mutated value fall within the domain, adopt mutated value
                //if not, retain the original value(mutation fails)
                if (pxm ==1){
                    x_m=dist_mut_x(gen);
                    if(x_m <= 0 or x_m >= 1){
                        x_m=x_r;
                    }
                }
                else if(pxm ==0){
                    x_m=x_r;
                }
                if (pym ==1){
                    y_m=dist_mut_y(gen);
                    if(y_m <= 0 or y_m >= 1){
                        y_m=y_r;
                    }
                }
                else if(pym ==0){
                    y_m=y_r;
                }
                if (pvm ==1){
                    v_m=dist_mut_v(gen);
                    if(v_m <= 0 or v_m >= 1){
                        v_m=v_r;
                    }
                }
                else if(pvm ==0){
                    v_m=v_r;
                }

                // comapare fitness between resident one and mutant one(judge the success of invasion)
                w_r=get_w_resident(x_r,y_r,v_r,beta,gamma,q_value);
                w_m=get_w_mutant(x_r,y_r,v_r,x_m,y_m,v_m,beta,gamma,q_value);
                // if fitness of mutant>fitness of resident, invasion success and update the resident values
                if(w_r<w_m){
                    x_r=x_m;
                    y_r=y_m;
                    v_r=v_m;
                }
            }
            // record results in arrays
            end_values2_095[beta_count][gamma_count][0]=beta;
            end_values2_095[beta_count][gamma_count][1]=gamma;
            end_values2_095[beta_count][gamma_count][2]=x_r;
            end_values2_095[beta_count][gamma_count][3]=y_r;
            end_values2_095[beta_count][gamma_count][4]=v_r;
            // reset the initial values
            x_i=0;
            y_i=0;
            v_i=0;
        }
    }

    //// q=0.9
    q_value=0.9;

    ///broad teaching ESS
    ///explore the parameter space
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        beta=beta_count*delta_beta;
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            gamma=gamma_count*delta_gamma;
            //initial values
            x_i=0.25;
            y_i=0.25;
            v_i=0.25;
            //set initial values to resident values
            x_r=x_i;
            y_r=y_i;
            v_r=v_i;
            ///numerical calculation(generate perturbation)
            for(int i=0; i<max_step; i++){
                bernoulli_distribution dist_pxm(mut_rate);
                bernoulli_distribution dist_pym(mut_rate);
                bernoulli_distribution dist_pvm(mut_rate);
                pxm = dist_pxm(gen);
                pym = dist_pym(gen);
                pvm = dist_pvm(gen);                    
                //slight mutation(normal distribution)
                normal_distribution<> dist_mut_x(x_r, mutation_step);
                normal_distribution<> dist_mut_y(y_r, mutation_step);
                normal_distribution<> dist_mut_v(v_r, mutation_step);
                //if mutated value fall within the domain, adopt mutated value
                //if not, retain the original value(mutation fails)
                if (pxm ==1){
                    x_m=dist_mut_x(gen);
                    if(x_m <= 0 or x_m >= 1){
                        x_m=x_r;
                    }
                }
                else if(pxm ==0){
                    x_m=x_r;
                }
                if (pym ==1){
                    y_m=dist_mut_y(gen);
                    if(y_m <= 0 or y_m >= 1){
                        y_m=y_r;
                    }
                }
                else if(pym ==0){
                    y_m=y_r;
                }            
                if (pvm ==1){
                    v_m=dist_mut_v(gen);
                    if(v_m <= 0 or v_m >= 1){
                        v_m=v_r;
                    }
                }
                else if(pvm ==0){
                    v_m=v_r;
                }

                // comapare fitness between resident one and mutant one(judge the success of invasion)
                w_r=get_w_resident(x_r,y_r,v_r,beta,gamma,q_value);
                w_m=get_w_mutant(x_r,y_r,v_r,x_m,y_m,v_m,beta,gamma,q_value);
                // if fitness of mutant>fitness of resident, invasion success and update the resident values
                if(w_r<w_m){
                    x_r=x_m;
                    y_r=y_m;
                    v_r=v_m;
                }
            }

            // record results in arrays
            end_values1_09[beta_count][gamma_count][0]=beta;
            end_values1_09[beta_count][gamma_count][1]=gamma;
            end_values1_09[beta_count][gamma_count][2]=x_r;
            end_values1_09[beta_count][gamma_count][3]=y_r;
            end_values1_09[beta_count][gamma_count][4]=v_r;
            // reset the initial values
            x_i=0;
            y_i=0;
            v_i=0;
        }
    }

    ///only teaching ESS
    //explore the parameter space
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        beta=beta_count*delta_beta;
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            gamma=gamma_count*delta_gamma;
            //initial values
            x_i=0.99;
            y_i=0.05;
            v_i=0.25;
            //set initial values to resident values
            x_r=x_i;
            y_r=y_i;
            v_r=v_i;
            ///numerical calculation(generate perturbation)
            for(int i=0; i<max_step; i++){
                bernoulli_distribution dist_pxm(mut_rate);
                bernoulli_distribution dist_pym(mut_rate);
                bernoulli_distribution dist_pvm(mut_rate);
                pxm = dist_pxm(gen);
                pym = dist_pym(gen);
                pvm = dist_pvm(gen);                    
                //slight mutation(normal distribution)
                normal_distribution<> dist_mut_x(x_r, mutation_step);
                normal_distribution<> dist_mut_y(y_r, mutation_step);
                normal_distribution<> dist_mut_v(v_r, mutation_step);
                //if mutated value fall within the domain, adopt mutated value
                //if not, retain the original value(mutation fails)
                if (pxm ==1){
                    x_m=dist_mut_x(gen);
                    if(x_m <= 0 or x_m >= 1){
                        x_m=x_r;
                    }
                }
                else if(pxm ==0){
                    x_m=x_r;
                }
                if (pym ==1){
                    y_m=dist_mut_y(gen);
                    if(y_m <= 0 or y_m >= 1){
                        y_m=y_r;
                    }
                }
                else if(pym ==0){
                    y_m=y_r;
                }
                if (pvm ==1){
                    v_m=dist_mut_v(gen);
                    if(v_m <= 0 or v_m >= 1){
                        v_m=v_r;
                    }
                }
                else if(pvm ==0){
                    v_m=v_r;
                }

                // comapare fitness between resident one and mutant one(judge the success of invasion)
                w_r=get_w_resident(x_r,y_r,v_r,beta,gamma,q_value);
                w_m=get_w_mutant(x_r,y_r,v_r,x_m,y_m,v_m,beta,gamma,q_value);
                // if fitness of mutant>fitness of resident, invasion success and update the resident values
                if(w_r<w_m){
                    x_r=x_m;
                    y_r=y_m;
                    v_r=v_m;
                }
            }

            // record results in arrays
            end_values2_09[beta_count][gamma_count][0]=beta;
            end_values2_09[beta_count][gamma_count][1]=gamma;
            end_values2_09[beta_count][gamma_count][2]=x_r;
            end_values2_09[beta_count][gamma_count][3]=y_r;
            end_values2_09[beta_count][gamma_count][4]=v_r;
            // reset the initial values
            x_i=0;
            y_i=0;
            v_i=0;
        }
    }

    //// q=0.8
    q_value=0.8;

    ///broad teaching ESS
    //explore the parameter space
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        beta=beta_count*delta_beta;
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            gamma=gamma_count*delta_gamma;
            //initial values
            x_i=0.25;
            y_i=0.25;
            v_i=0.25;
            //set initial values to resident values
            x_r=x_i;
            y_r=y_i;
            v_r=v_i;
            ///numerical calculation(generate perturbation)
            for(int i=0; i<max_step; i++){
                bernoulli_distribution dist_pxm(mut_rate);
                bernoulli_distribution dist_pym(mut_rate);
                bernoulli_distribution dist_pvm(mut_rate);
                pxm = dist_pxm(gen);
                pym = dist_pym(gen);
                pvm = dist_pvm(gen);                    
                //slight mutation(normal distribution)
                normal_distribution<> dist_mut_x(x_r, mutation_step);
                normal_distribution<> dist_mut_y(y_r, mutation_step);
                normal_distribution<> dist_mut_v(v_r, mutation_step);
                //if mutated value fall within the domain, adopt mutated value
                //if not, retain the original value(mutation fails)
                if (pxm ==1){
                    x_m=dist_mut_x(gen);
                    if(x_m <= 0 or x_m >= 1){
                        x_m=x_r;
                    }
                }
                else if(pxm ==0){
                    x_m=x_r;
                }
                if (pym ==1){
                    y_m=dist_mut_y(gen);
                    if(y_m <= 0 or y_m >= 1){
                        y_m=y_r;
                    }
                }
                else if(pym ==0){
                    y_m=y_r;
                }
                if (pvm ==1){
                    v_m=dist_mut_v(gen);
                    if(v_m <= 0 or v_m >= 1){
                        v_m=v_r;
                    }
                }
                else if(pvm ==0){
                    v_m=v_r;
                }

                // comapare fitness between resident one and mutant one(judge the success of invasion)
                w_r=get_w_resident(x_r,y_r,v_r,beta,gamma,q_value);
                w_m=get_w_mutant(x_r,y_r,v_r,x_m,y_m,v_m,beta,gamma,q_value);
                // if fitness of mutant>fitness of resident, invasion success and update the resident values
                if(w_r<w_m){
                    x_r=x_m;
                    y_r=y_m;
                    v_r=v_m;
                }
            }
            // record results in arrays
            end_values1_08[beta_count][gamma_count][0]=beta;
            end_values1_08[beta_count][gamma_count][1]=gamma;
            end_values1_08[beta_count][gamma_count][2]=x_r;
            end_values1_08[beta_count][gamma_count][3]=y_r;
            end_values1_08[beta_count][gamma_count][4]=v_r;
            // reset the initial values
            x_i=0;
            y_i=0;
            v_i=0;
        }
    }

    ///only teaching ESS
    //explore the parameter space
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        beta=beta_count*delta_beta;
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            gamma=gamma_count*delta_gamma;
            //initial values
            x_i=0.99;
            y_i=0.05;
            v_i=0.25;
            //set initial values to resident values
            x_r=x_i;
            y_r=y_i;
            v_r=v_i;
            ///numerical calculation(generate perturbation)
            for(int i=0; i<max_step; i++){
                bernoulli_distribution dist_pxm(mut_rate);
                bernoulli_distribution dist_pym(mut_rate);
                bernoulli_distribution dist_pvm(mut_rate);
                pxm = dist_pxm(gen);
                pym = dist_pym(gen);
                pvm = dist_pvm(gen);                    
                //slight mutation(normal distribution)
                normal_distribution<> dist_mut_x(x_r, mutation_step);
                normal_distribution<> dist_mut_y(y_r, mutation_step);
                normal_distribution<> dist_mut_v(v_r, mutation_step);
                //if mutated value fall within the domain, adopt mutated value
                //if not, retain the original value(mutation fails)
                if (pxm ==1){
                    x_m=dist_mut_x(gen);
                    if(x_m <= 0 or x_m >= 1){
                        x_m=x_r;
                    }
                }
                else if(pxm ==0){
                    x_m=x_r;
                }
                if (pym ==1){
                    y_m=dist_mut_y(gen);
                    if(y_m <= 0 or y_m >= 1){
                        y_m=y_r;
                    }
                }
                else if(pym ==0){
                    y_m=y_r;
                }
                if (pvm ==1){
                    v_m=dist_mut_v(gen);
                    if(v_m <= 0 or v_m >= 1){
                        v_m=v_r;
                    }
                }
                else if(pvm ==0){
                    v_m=v_r;
                }

                // comapare fitness between resident one and mutant one(judge the success of invasion)
                w_r=get_w_resident(x_r,y_r,v_r,beta,gamma,q_value);
                w_m=get_w_mutant(x_r,y_r,v_r,x_m,y_m,v_m,beta,gamma,q_value);
                // if fitness of mutant>fitness of resident, invasion success and update the resident values
                if(w_r<w_m){
                    x_r=x_m;
                    y_r=y_m;
                    v_r=v_m;
                }
            }
            // record results in arrays
            end_values2_08[beta_count][gamma_count][0]=beta;
            end_values2_08[beta_count][gamma_count][1]=gamma;
            end_values2_08[beta_count][gamma_count][2]=x_r;
            end_values2_08[beta_count][gamma_count][3]=y_r;
            end_values2_08[beta_count][gamma_count][4]=v_r;
            // reset the initial values
            x_i=0;
            y_i=0;
            v_i=0;
        }
    }

    //// export and save results
    ofstream fout10("data_q10.csv");
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            double val_te  = end_values1_10[beta_count][gamma_count][3];
            double val_te2 = end_values2_10[beta_count][gamma_count][3];
            double val_il2 = end_values2_10[beta_count][gamma_count][2];
            double val_il1 = end_values1_10[beta_count][gamma_count][2];
            int map =(val_te > 0.0001 && val_te2 >= 0.00001) ? 2 :
                    (val_te > 0.0001 && val_te2 < 0.00001 && val_il2 >=0.999) ? 3 :
                    (val_te > 0.0001 && val_te2 < 0.00001 && val_il2 < 1) ? 4 :
                    (val_te < 0.0001 && val_il1 > 0.9999) ? 1: 0;
            fout10 << beta_count << "," << gamma_count << "," << 
            end_values1_10[beta_count][gamma_count][2] << "," << end_values1_10[beta_count][gamma_count][3] << ","  << end_values1_10[beta_count][gamma_count][4] << "," <<
            end_values2_10[beta_count][gamma_count][2] << "," << end_values2_10[beta_count][gamma_count][3] << ","  << end_values2_10[beta_count][gamma_count][4] << "," <<
            map << "\n";
        }
    }
    fout10.close();

    ofstream fout095("data_q095.csv");
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            double val_te  = end_values1_095[beta_count][gamma_count][3];
            double val_te2 = end_values2_095[beta_count][gamma_count][3];
            double val_il2 = end_values2_095[beta_count][gamma_count][2];
            double val_il1 = end_values1_095[beta_count][gamma_count][2];
            int map =(val_te > 0.0001 && val_te2 >= 0.00001) ? 2 :
                    (val_te > 0.0001 && val_te2 < 0.00001 && val_il2 >=0.999) ? 3 :
                    (val_te > 0.0001 && val_te2 < 0.00001 && val_il2 < 1) ? 4 :
                    (val_te < 0.0001 && val_il1 > 0.9999) ? 1: 0;
            fout095 << beta_count << "," << gamma_count << "," << 
            end_values1_095[beta_count][gamma_count][2] << "," << end_values1_095[beta_count][gamma_count][3] << ","  << end_values1_095[beta_count][gamma_count][4] << "," <<
            end_values2_095[beta_count][gamma_count][2] << "," << end_values2_095[beta_count][gamma_count][3] << ","  << end_values2_095[beta_count][gamma_count][4] << "," <<
            map << "\n";
        }
    }
    fout095.close();

    ofstream fout09("data_q09.csv");
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            double val_te  = end_values1_09[beta_count][gamma_count][3];
            double val_te2 = end_values2_09[beta_count][gamma_count][3];
            double val_il2 = end_values2_09[beta_count][gamma_count][2];
            double val_il1 = end_values1_09[beta_count][gamma_count][2];
            int map =(val_te > 0.0001 && val_te2 >= 0.00001) ? 2 :
                    (val_te > 0.0001 && val_te2 < 0.00001 && val_il2 >=0.999) ? 3 :
                    (val_te > 0.0001 && val_te2 < 0.00001 && val_il2 < 1) ? 4 :
                    (val_te < 0.0001 && val_il1 > 0.9999) ? 1: 0;
            fout09 << beta_count << "," << gamma_count << "," << 
            end_values1_09[beta_count][gamma_count][2] << "," << end_values1_09[beta_count][gamma_count][3] << ","  << end_values1_09[beta_count][gamma_count][4] << "," <<
            end_values2_09[beta_count][gamma_count][2] << "," << end_values2_09[beta_count][gamma_count][3] << ","  << end_values2_09[beta_count][gamma_count][4] << "," <<
            map << "\n";
        }
    }
    fout09.close();

    ofstream fout08("data_q08.csv");
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            double val_te  = end_values1_08[beta_count][gamma_count][3];
            double val_te2 = end_values2_08[beta_count][gamma_count][3];
            double val_il2 = end_values2_08[beta_count][gamma_count][2];
            double val_il1 = end_values1_08[beta_count][gamma_count][2];
            int map =(val_te > 0.0001 && val_te2 >= 0.00001) ? 2 :
                    (val_te > 0.0001 && val_te2 < 0.00001 && val_il2 >=0.999) ? 3 :
                    (val_te > 0.0001 && val_te2 < 0.00001 && val_il2 < 1) ? 4 :
                    (val_te < 0.0001 && val_il1 > 0.9999) ? 1: 0;
            fout08 << beta_count << "," << gamma_count << "," << 
            end_values1_08[beta_count][gamma_count][2] << "," << end_values1_08[beta_count][gamma_count][3] << ","  << end_values1_08[beta_count][gamma_count][4] << "," <<
            end_values2_08[beta_count][gamma_count][2] << "," << end_values2_08[beta_count][gamma_count][3] << ","  << end_values2_08[beta_count][gamma_count][4] << "," <<
            map << "\n";
        }
    }
    fout08.close();

    //// visualize using gunuplot
    FILE *gnuplot = popen("gnuplot -persistent", "w");
    
    /// set output style(png style)
    //fprintf(gnuplot, "set terminal pngcairo size 1200,1200 enhanced font 'Arial,16'\n");
    //fprintf(gnuplot, "set output 'Figure6.png'\n");
    /// set output style(eps style)
    fprintf(gnuplot, "set terminal epscairo size 12,12 enhanced font 'Arial,24'\n");
    fprintf(gnuplot, "set output 'Figure6.eps'\n");

    fprintf(gnuplot, "set multiplot layout 2,2\n");

    ///q=1
    fprintf(gnuplot, "set size 0.45,0.45\n");
    fprintf(gnuplot, "set origin 0.05,0.5\n");
    fprintf(gnuplot, "set label '(a) {/:Italic q} = 1' at screen 0.03, 0.92 \n");
    fprintf(gnuplot, "set xlabel \'Effect of teaching, {/:Italic {/Symbol g}}\' offset 0, 0.5\n");
    fprintf(gnuplot, "set ylabel \'Efficiency of observation learning, {/:Italic {/Symbol b}}\' offset -2,0\n");
    fprintf(gnuplot, "set tics out \n");
    fprintf(gnuplot, "set tics scale 0.5 \n");
    fprintf(gnuplot, "set xtics offset 0, 0.5\n");

    // setting for heatmap views
    fprintf(gnuplot, "set view map\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "set palette defined (0 '#FFFFFF', 1 '#808080', 2 '#90EE90', 3 '#AFEEEE', 4 '#008B8B')\n");            // color palette（0:white, 1:whitesmoke, 2:lightgreen, 3:lightcyan, 4:darkcyan）
    fprintf(gnuplot, "unset colorbox\n");  
    
    //scaling ticks
    fprintf(gnuplot, "set ytics (");
    for (beta_count=0; beta_count<=beta_max/delta_beta; beta_count+=5/delta_beta){
        fprintf(gnuplot, "\"%.0f\" %d", beta_count*delta_beta, beta_count);
        if (beta_count != beta_max/delta_beta+1) fprintf(gnuplot, ", ");
    }
    fprintf(gnuplot, ")\n");
    fprintf(gnuplot, "set xtics (");
    for (gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count+=5/delta_beta){
        fprintf(gnuplot, "\"%.0f\" %d", gamma_count*delta_gamma, gamma_count);
        if (gamma_count != gamma_max/delta_gamma+1) fprintf(gnuplot, ", ");
    }
    fprintf(gnuplot, ")\n");

    /// send data in data block format
    // 1)declare data block
    fprintf(gnuplot, "$DATAQ10 << EOD\n");
    // 2)export result
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        beta=beta_count*delta_beta;
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            gamma=gamma_count*delta_gamma;
            // teaching evolves (initial value 1)
            double val_te=end_values1_10[beta_count][gamma_count][3];
            // teaching evolves (initial value 2)
            double val_te2=end_values2_10[beta_count][gamma_count][3];
            // social learning & teaching do not evolve (initial value 2)
            double val_il2=end_values2_10[beta_count][gamma_count][2];
            // social learning & teaching do not evolve (initial value 1)
            double val_il1=end_values1_10[beta_count][gamma_count][2];
            //classify region
            int map =(val_te > 0.0001 && val_te2 >= 0.00001) ? 2 : (val_te > 0.0001 && val_te2 < 0.00001 && val_il2 >=0.999) ? 3 : (val_te > 0.0001 && val_te2 < 0.00001 && val_il2 < 1) ? 4 : (val_te < 0.0001 && val_il1 > 0.9999) ? 1: 0;  
            fprintf(gnuplot, "%d ", map);
        }
        fprintf(gnuplot, "\n");
    }
    // 3)end data block
    fprintf(gnuplot, "EOD\n");
    // 4)plot data block using splot
    fprintf(gnuplot, "splot $DATAQ10 matrix with image\n");
    fflush(gnuplot);
    
    ///q=0.95
    fprintf(gnuplot, "set size 0.45,0.45\n");
    fprintf(gnuplot, "set origin 0.52,0.5\n");
    fprintf(gnuplot, "set label '(b) {/:Italic q} = 0.95' at screen 0.5, 0.92\n");
    fprintf(gnuplot, "set xlabel \'Effect of teaching, {/:Italic {/Symbol g}}\' \n");
    fprintf(gnuplot, "set ylabel \'Efficiency of observation learning, {/:Italic {/Symbol b}}\' offset -2,0 \n");
    fprintf(gnuplot, "set tics out \n");
    fprintf(gnuplot, "set tics scale 0.5 \n");

    // setting for heatmap views
    fprintf(gnuplot, "set view map\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "set palette defined (0 '#FFFFFF', 1 '#808080', 2 '#90EE90', 3 '#AFEEEE', 4 '#008B8B')\n");            // color palette（0:white, 1:whitesmoke, 2:lightgreen, 3:lightcyan, 4:darkcyan）
    fprintf(gnuplot, "unset colorbox\n");                           // カラーバーオフ    
    
    //scaling ticks
    fprintf(gnuplot, "set ytics (");
    for (beta_count=0; beta_count<=beta_max/delta_beta; beta_count+=5/delta_beta){
        fprintf(gnuplot, "\"%.0f\" %d", beta_count*delta_beta, beta_count);
        if (beta_count != beta_max/delta_beta+1) fprintf(gnuplot, ", ");
    }
    fprintf(gnuplot, ")\n");
    fprintf(gnuplot, "set xtics (");
    for (gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count+=5/delta_beta){
        fprintf(gnuplot, "\"%.0f\" %d", gamma_count*delta_gamma, gamma_count);
        if (gamma_count != gamma_max/delta_gamma+1) fprintf(gnuplot, ", ");
    }
    fprintf(gnuplot, ")\n");

    /// send data in data block format
    // 1)declare data block
    fprintf(gnuplot, "$DATAQ095 << EOD\n");
    // 2)export result
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        beta=beta_count*delta_beta;
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            gamma=gamma_count*delta_gamma;
            // teaching evolves (initial value 1)
            double val_te=end_values1_095[beta_count][gamma_count][3];
            // teaching evolves (initial value 2)
            double val_te2=end_values2_095[beta_count][gamma_count][3];
            // social learning & teaching do not evolve (initial value 2)
            double val_il2=end_values2_095[beta_count][gamma_count][2];
            // social learning & teaching do not evolve (initial value 1)
            double val_il1=end_values1_095[beta_count][gamma_count][2];
            // classify region
            int map =(val_te > 0.0001 && val_te2 >= 0.00001) ? 2 : (val_te > 0.0001 && val_te2 < 0.00001 && val_il2 >=0.999) ? 3 : (val_te > 0.0001 && val_te2 < 0.00001 && val_il2 < 1) ? 4 : (val_te < 0.0001 && val_il1 > 0.9999) ? 1: 0;  
            fprintf(gnuplot, "%d ", map);
        }
        fprintf(gnuplot, "\n");
    }
    // 3)end data block
    fprintf(gnuplot, "EOD\n");
    // 4)plot data block using splot
    fprintf(gnuplot, "splot $DATAQ095 matrix with image\n");
    fflush(gnuplot);
    
    ////q=0.9
    fprintf(gnuplot, "set size 0.45,0.45\n");
    fprintf(gnuplot, "set origin 0.05,0.05\n");
    fprintf(gnuplot, "set tics front\n");
    fprintf(gnuplot, "set label '(c) {/:Italic q} = 0.9' at screen 0.03, 0.47 front\n");
    fprintf(gnuplot, "set xlabel \'Effect of teaching, {/:Italic {/Symbol g}}\' front\n");
    fprintf(gnuplot, "set ylabel \'Efficiency of observation learning, {/:Italic {/Symbol b}}\' offset -2,0 front\n");
    fprintf(gnuplot, "set tics out \n");
    fprintf(gnuplot, "set tics scale 0.5 \n");

    // setting for heatmap views
    fprintf(gnuplot, "set view map\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "set palette defined (0 '#FFFFFF', 1 '#808080', 2 '#90EE90', 3 '#AFEEEE', 4 '#008B8B')\n");            // color palette（0:white, 1:whitesmoke, 2:lightgreen, 3:lightcyan, 4:darkcyan）
    fprintf(gnuplot, "unset colorbox\n");   
    
    //scaling ticks
    fprintf(gnuplot, "set ytics (");
    for (beta_count=0; beta_count<=beta_max/delta_beta; beta_count+=5/delta_beta){
        fprintf(gnuplot, "\"%.0f\" %d", beta_count*delta_beta, beta_count);
        if (beta_count != beta_max/delta_beta+1) fprintf(gnuplot, ", ");
    }
    fprintf(gnuplot, ")\n");
    fprintf(gnuplot, "set xtics (");
    for (gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count+=5/delta_beta){
        fprintf(gnuplot, "\"%.0f\" %d", gamma_count*delta_gamma, gamma_count);
        if (gamma_count != gamma_max/delta_gamma+1) fprintf(gnuplot, ", ");
    }
    fprintf(gnuplot, ")\n");

    /// send data in data block format
    // 1)declare data block
    fprintf(gnuplot, "$DATAQ090 << EOD\n");
    // 2)export result
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        beta=beta_count*delta_beta;
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            gamma=gamma_count*delta_gamma;
            // teaching evolves (initial value 1)
            double val_te=end_values1_09[beta_count][gamma_count][3];
            // teaching evolves (initial value 2)
            double val_te2=end_values2_09[beta_count][gamma_count][3];
            // social learning & teaching do not evolve (initial value 2)
            double val_il2=end_values2_09[beta_count][gamma_count][2];
            // social learning & teaching do not evolve (initial value 1)
            double val_il1=end_values1_09[beta_count][gamma_count][2];
            //classify region
            int map =(val_te > 0.0001 && val_te2 >= 0.00001) ? 2 : (val_te > 0.0001 && val_te2 < 0.00001 && val_il2 >=0.999) ? 3 : (val_te > 0.0001 && val_te2 < 0.00001 && val_il2 < 1) ? 4 : (val_te < 0.0001 && val_il1 > 0.9999) ? 1: 0;  
            fprintf(gnuplot, "%d ", map);
        }
        fprintf(gnuplot, "\n");
    }
    // 3)end data block
    fprintf(gnuplot, "EOD\n");
    // 4)plot data block using splot
    fprintf(gnuplot, "splot $DATAQ090 matrix with image\n");
    fflush(gnuplot);

    ///q=0.8
    fprintf(gnuplot, "set size 0.45,0.45\n");
    fprintf(gnuplot, "set origin 0.52,0.05\n");
    fprintf(gnuplot, "set label '(d) {/:Italic q} = 0.8' at screen 0.5, 0.47\n");
    fprintf(gnuplot, "set xlabel \'Effect of teaching, {/:Italic {/Symbol g}}\' \n");
    fprintf(gnuplot, "set ylabel \'Efficiency of observation learning, {/:Italic {/Symbol b}}\' offset -2,0 \n");
    fprintf(gnuplot, "set tics out \n");
    fprintf(gnuplot, "set tics scale 0.5 \n");

    // setting for heatmap views
    fprintf(gnuplot, "set view map\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "set palette defined (0 '#FFFFFF', 1 '#808080', 2 '#90EE90', 3 '#AFEEEE', 4 '#008B8B')\n");            // color palette（0:white, 1:whitesmoke, 2:lightgreen, 3:lightcyan, 4:darkcyan）
    fprintf(gnuplot, "unset colorbox\n");   
    
    //scaling ticks
    fprintf(gnuplot, "set ytics (");
    for (beta_count=0; beta_count<=beta_max/delta_beta; beta_count+=5/delta_beta){
        fprintf(gnuplot, "\"%.0f\" %d", beta_count*delta_beta, beta_count);
        if (beta_count != beta_max/delta_beta+1) fprintf(gnuplot, ", ");
    }
    fprintf(gnuplot, ")\n");
    fprintf(gnuplot, "set xtics (");
    for (gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count+=5/delta_beta){
        fprintf(gnuplot, "\"%.0f\" %d", gamma_count*delta_gamma, gamma_count);
        if (gamma_count != gamma_max/delta_gamma+1) fprintf(gnuplot, ", ");
    }
    fprintf(gnuplot, ")\n");

    /// send data in data block format
    // 1)declare data block
    fprintf(gnuplot, "$DATAQ080 << EOD\n");
    // 2)export result
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        beta=beta_count*delta_beta;
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            gamma=gamma_count*delta_gamma;
            // teaching evolves (initial value 1)
            double val_te=end_values1_08[beta_count][gamma_count][3];
            // teaching evolves (initial value 2)
            double val_te2=end_values2_08[beta_count][gamma_count][3];
            // social learning & teaching do not evolve (initial value 2)
            double val_il2=end_values2_08[beta_count][gamma_count][2];
            // social learning & teaching do not evolve (initial value 1)
            double val_il1=end_values1_08[beta_count][gamma_count][2];
            //classify region
            int map =(val_te > 0.0001 && val_te2 >= 0.00001) ? 2 : (val_te > 0.0001 && val_te2 < 0.00001 && val_il2 >=0.999) ? 3 : (val_te > 0.0001 && val_te2 < 0.00001 && val_il2 < 1) ? 4 : (val_te < 0.0001 && val_il1 > 0.9999) ? 1: 0;  
            fprintf(gnuplot, "%d ", map);
        }
        fprintf(gnuplot, "\n");
    }
    // 3)end data block
    fprintf(gnuplot, "EOD\n");
    // 4)plot data block using splot
    fprintf(gnuplot, "splot $DATAQ080 matrix with image\n");
    fflush(gnuplot);

    fprintf(gnuplot, "unset multiplot\n");
    fprintf(gnuplot, "unset output\n");
}