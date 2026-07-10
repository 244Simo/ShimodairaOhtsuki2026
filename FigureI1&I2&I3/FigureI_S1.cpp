# include <stdio.h>
# include <iostream>
# include <random>
# include <cmath>
# include <fstream>
# include <vector>
using namespace std;

////Function prototype
double get_z();
double get_w();
////Function
double get_z(double yf, double vf, double zf, double xp, double yp, double vp, double beta, double gamma){
    double z;
    z = zf*(1-exp(-(beta+gamma*yf*(1-vf))*vp*(1-xp)))+xp*vp;
    return z;
}
double get_w(double yf, double vf, double zf, double xp, double yp, double vp, double beta, double gamma, double epsilon){
    double z;
    double w;
    z = zf*(1-exp(-(beta+gamma*yf*(1-vf))*vp*(1-xp)))+xp*vp;
    w = 1+epsilon*(z*(1-yp)*(1-vp));
    return w;
}

//// main program
int main(void){
    double beta;
    double gamma;
    double q_value;
    double epsilon;

    double mut_rate;        //mutation rate of (v,x,y)
    double mutation_step;   //standard deviation of mutations
    int pxm;    //determine whether mutation happens on x(if pxm=1 -> mutation/ pxm=0 -> no-mutation)
    int pym;    //determine whether mutation happens on y(if pym=1 -> mutation/ pym=0 -> no-mutation)
    int pvm;    //determine whether mutation happens on v(if pvm=1 -> mutation/ pvm=0 -> no-mutation)
    int vtp;    //determine whether learner chooses own parent as a role model(if vtp=1 -> parent/ vtp=0 -> non-parent)
    int role_model_index;   // index of a role model

    int N;      //number of agent(individual)
    int g_max;  //steps on simulation

    ///variants for calculation
    //(x,y,v)-values of an individual[index] i in the focal generation
    double x_i;
    double y_i;
    double v_i;
    //(x,y,v)-values of an individual[index] i in the previous generation
    double xf_i;
    double yf_i;
    double vf_i;
    double zf_i;
    //(x,y,v)-values of an individual[index] i in the focal generation after mutation
    double xm_i;
    double ym_i;
    double vm_i;

    /// parameters
    // For figure I1
    beta=1;
    gamma=20;
    q_value=0.9;
    epsilon=0.1;
    N=1000;
    g_max=1*pow(10,5);
    mut_rate=0.1;
    mutation_step=0.01;

    //// generate random number
    // if you want to set seed value mamually
    //int seed = 0;
    //mt19937 gen(seed);

    // if you want to change seed value randomly
    random_device rd;
    mt19937 gen(rd());

    /// array
    // focal generation array(index of an individual, traits of each individual)
    vector< vector<double> > individuals(N, vector<double>(7));
    // parental(previous) generation array(index of an individual, traits of each individual)
    vector< vector<double> > prev_individuals(N, vector<double>(7));
    // all generation array(nuber of generation, index of an individual, traits of each individual)
    vector< vector <vector <double> > > generation_list(g_max,vector<vector<double> > (N,vector<double>(7)));

    ///arraies related to generation update(Wright-Fisher like process)
    //fitness array(fitness of individuals in the parental generation)
    vector<double> ind_fit(N);
    //array for recording indexes of parents in focal generation(parents are chosen based on weighted fitness)
    vector<double> parent_list(N);

    // array for random oblique transmission when parents does not chosen as role models
    vector< vector<double> > random_rolemodels(N, vector<double>(7));

    ///set initial values for each individual
    //Initial values are given based on a uniform distribution
    uniform_real_distribution<double> dist1(0.0,1.0);
    //put values into the array
    //0:x, 1:y, 2:v, 3:z, 4:w, 5:rate of vertical transmission, 6:parent number(previous generation)
    for(int i=0; i<N; i++){
        // put (x,y,v) values
        for(int j=0; j<3; j++){
            prev_individuals[i][j]=dist1(gen);
        }
        // obtain (x,y,v) values for calculation
        x_i=prev_individuals[i][0];
        y_i=prev_individuals[i][1];
        v_i=prev_individuals[i][2];
        // Since there are no role models in the first(initial) generation,
        // we assume that the individuals in the first generation have parents who have the same strategy
        // and the individuals always choose own parents as a role model (vertical transmission)
        prev_individuals[i][3]=x_i*v_i*exp(-(beta+gamma*y_i*(1-v_i))*(v_i*(1-x_i)));
        prev_individuals[i][4]=(1-v_i)*(1-y_i)*prev_individuals[i][3];
        prev_individuals[i][5]=1;
        prev_individuals[i][6]=i;
        // reset variants
        x_i=0,y_i=0,v_i=0;
    }
    // put information of initial generation to all generation array
    for(int i=0; i<N; i++){
        for(int j=0; j<7; j++){
            generation_list[0][i][j]=prev_individuals[i][j];
        }
    }

    /// update generations
    for(int g=1; g<g_max; g++){
        // For checking the progress during the simulation 
        /*
        if(g % 1000 ==0){
            cout << "generation is" << g << endl;
        }
        */
        
        /// determine individuals in the previous generation who reproduce
        // update array of each individual
        for(int i=0; i<N; i++){
            ind_fit[i]=prev_individuals[i][4];
        }
        // return index weighted by fitness(Wright-Fisher like update)
        discrete_distribution<> WFP(ind_fit.begin(), ind_fit.end());
        for(int i=0; i<N; i++){
            parent_list[i]=WFP(gen);
        }
        /// inheritance, learning, and obtain fitness
        // inherit strategy(x,y,v) from parents(previous generation) to own offspring(focal generation) -> mutation in offspring
        // calculate z-value (cultural level) and w-value(fitness)
        for(int i=0; i<N; i++){
            //1. record the information that who is parent of an individuals[i]
            int parent_index=parent_list[i];
            individuals[i][6]=parent_index;
            //2. inherit strategy (x,y,v) from the parent(prev_individuals) to the individuals[i]
            xf_i=prev_individuals[parent_index][0];
            yf_i=prev_individuals[parent_index][1];
            vf_i=prev_individuals[parent_index][2];
            individuals[i][0]=xf_i;
            individuals[i][1]=yf_i;
            individuals[i][2]=vf_i;
            //3. obtaing x,y,v value of the individuals[i] before mutation
            x_i=individuals[i][0];
            y_i=individuals[i][1];
            v_i=individuals[i][2];
            //4. mutation
            //4-1. determine whether mutation happens or not based on bernoulli distribution
            bernoulli_distribution dist_pxm(mut_rate);
            bernoulli_distribution dist_pym(mut_rate);
            bernoulli_distribution dist_pvm(mut_rate);
            pxm = dist_pxm(gen);
            pym = dist_pym(gen);
            pvm = dist_pvm(gen);
            //4-2. if pxm=1 or pym=1 or pvm=1, mutation happens
            //if mutated value fall within the domain, adopt mutated value
            //if not, retain the original value(mutation fails
            //x mutation
            if(pxm ==1){
                // normal distribution (mean: x_i, standard deviation: mutation_step)
                normal_distribution<> dist_mut_x(x_i, mutation_step);
                xm_i = dist_mut_x(gen);
                // if mutated value(xm_i) fall within the domain, we update the strategy(x_i)
                if(0<=xm_i && xm_i<=1){
                    individuals[i][0]=xm_i;
                    x_i=xm_i;
                }
            }
            //y mutation
            if(pym ==1){
                // // normal distribution (mean: y_i, standard deviation: mutation_step)
                normal_distribution<> dist_mut_y(y_i, mutation_step);
                ym_i = dist_mut_y(gen);
                // if mutated value(ym_i) fall within the domain, we update the strategy(y_i)
                if(0<=ym_i && ym_i<=1){
                    individuals[i][1]=ym_i;
                    y_i=ym_i;
                }
            }
            //v mutation
            if(pvm ==1){
                // normal distribution (mean: v_i, standard deviation: mutation_step)
                normal_distribution<> dist_mut_v(v_i, mutation_step);
                vm_i = dist_mut_v(gen);
                // if mutated value(vm_i) fall within the domain, we update the strategy(v_i)
                if(0<=vm_i && vm_i<=1){
                    individuals[i][2]=vm_i;
                    v_i=vm_i;
                }
            }
            //5. determine whether the individuals[i] chooses own parent as a role model based on bernoulli distribution
            bernoulli_distribution dist_vtrate(q_value);
            vtp=dist_vtrate(gen);
            individuals[i][5]=vtp;  //vtp=1 -> choose a parent, vtp=0 -> choose a random individual except own parent 
            //5-1. Choose own parent as a role model
            if(vtp==1){
                //obtain the z-value of the parent(role model)(x,y,v values have already updated in the process 2)
                zf_i=prev_individuals[parent_index][3];
                //calculate z-value and w-value of the individuals[i](in the focal generation)
                individuals[i][3]=get_z(yf_i,vf_i,zf_i,x_i,y_i,v_i,beta,gamma);
                individuals[i][4]=get_w(yf_i,vf_i,zf_i,x_i,y_i,v_i,beta,gamma,epsilon);
            }
            //5-2. Choose non-parent as a role model
            else if(vtp==0){
                // prepare array of individuals in parental generation exept own parent
                random_rolemodels=prev_individuals;
                random_rolemodels.erase(random_rolemodels.begin()+individuals[i][6]);
                //randomly choose a role model based on uniform distibution
                uniform_int_distribution<int> dist_role_model(0,N-2);
                role_model_index=dist_role_model(gen);
                //obtain the y,v,z- values of the role model
                //(note: we cannot neglect y&v-values in scenario 1 because role models always perform teaching)
                zf_i=random_rolemodels[role_model_index][3];
                yf_i=random_rolemodels[role_model_index][1];
                vf_i=random_rolemodels[role_model_index][2];
                //calculate z-value and w-value of the individuals[i](in the focal generation)
                individuals[i][3]=get_z(yf_i,vf_i,zf_i,x_i,y_i,v_i,beta,gamma);
                individuals[i][4]=get_w(yf_i,vf_i,zf_i,x_i,y_i,v_i,beta,gamma, epsilon);     
            }
        }   //all individuals in the focal generation has updated

        /// record values(information) in array
        //update previous generation array (information of focal generation array in one generation will be used as an information of previous generation in next generation)
        prev_individuals=individuals;
        // update all generation array
        for(int i=0; i<N; i++){
            for(int j=0; j<7; j++){
                generation_list[g][i][j]=individuals[i][j];
            }
        }
        // reset the focal generation array(this array is used to record results during one generation)
        for(int i=0; i<N; i++){
            for(int j=0; j<7; j++){
                individuals[i][j]=0;
            }
        }
    }

    ///export the result in csv file
    // generation_list[g][i][...]
    ofstream fout("timelapse_data.csv");
    int x_min = 0;
    int x_max = g_max-1;
    for(int s = x_min; s < x_max; s++){
        for(int i = 0; i < N; i++){
            double sl = (1-generation_list[s][i][0]) * generation_list[s][i][2];
            double il = generation_list[s][i][0] * generation_list[s][i][2];
            double ex = (1-generation_list[s][i][1]) * (1-generation_list[s][i][2]);
            double te = generation_list[s][i][1] * (1-generation_list[s][i][2]);
            double z  = generation_list[s][i][3];
            double w  = generation_list[s][i][4];
            fout << s << ","
                << sl << ","
                << il << ","
                << ex << ","
                << te << ","
                << z  << ","
                << w  << "\n";
        }
    }
    fout.close();
    // end program
    return 0;
}




