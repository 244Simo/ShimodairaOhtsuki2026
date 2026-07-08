# include <stdio.h>
# include <iostream>
# include <random>
# include <cmath>
# include <fstream>
# include <vector>
using namespace std;

////Function prototype
//For scenario1(General teaching)
double get_w_resident1();
double get_w_mutant1();
double get_z_resident1();
//For scenario2(Exclusive teaching)
double get_w_resident2();
double get_w_mutant2();
double get_z_resident2();
////Function
//For scenario1(General teaching)
double get_w_resident1(double x, double y, double v, double beta, double gamma, double q){
    double w_r;
    w_r = 1+0.1*(x*v*(1-y)*(1-v)*exp((beta+gamma*y*(1-v))*(1-x)*v));
    return w_r;
}
double get_w_mutant1(double x_r, double y_r, double v_r, double x_m, double y_m, double v_m, double beta, double gamma, double q){
    double z_r;
    double z_m;
    double w_m;
    z_r = x_r*v_r*exp((beta+gamma*y_r*(1-v_r))*(1-x_r)*v_r);
    z_m = ((1-q)*(1-exp(-(beta+gamma*y_r*(1-v_r))*v_m*(1-x_m)))*z_r+x_m*v_m)/(1-q*(1-exp(-(beta+gamma*y_m*(1-v_m))*v_m*(1-x_m))));
    w_m = 1+0.1*((1-y_m)*(1-v_m)*z_m);
    return w_m;
}
double get_z_resident1(double x, double y, double v, double beta, double gamma, double q){
    double z_r;
    z_r = x*v*exp((beta+gamma*y*(1-v))*(1-x)*v);
    return z_r;
}
//For scenario2(Exclusive teaching)
double get_w_resident2(double x, double y, double v, double beta, double gamma, double q){
    double w_r;
    w_r = 1+0.1*(x*v*(1-y)*(1-v)*exp(beta*(1-x)*v)*1/(1-q*(1-exp(-gamma*y*(1-v)*(1-x)*v))));
    return w_r;
}
double get_w_mutant2(double x_r, double y_r, double v_r, double x_m, double y_m, double v_m, double beta, double gamma, double q){
    double z_r;
    double z_m;
    double w_m;
    z_r = x_r*v_r*exp(beta*(1-x_r)*v_r)*1/(1-q*(1-exp(-gamma*y_r*(1-v_r)*(1-x_r)*v_r)));
    z_m = ((1-q)*(1-exp(-beta*v_m*(1-x_m)))*z_r+x_m*v_m)/(1-q*(1-exp(-(beta+gamma*y_m*(1-v_m))*v_m*(1-x_m))));
    w_m = 1+0.1*((1-y_m)*(1-v_m)*z_m);
    return w_m;
}
double get_z_resident2(double x, double y, double v, double beta, double gamma, double q){
    double z_r;
    z_r = x*v*exp(beta*(1-x)*v)*1/(1-q*(1-exp(-gamma*y*(1-v)*(1-x)*v)));
    return z_r;
}

//// main program
int main(void){
    int beta_count; //counter for loop
    int gamma_count;//counter for loop
    int beta_max;   //max value of β
    int gamma_max;  //max value of γ
    double delta_beta;
    double delta_gamma;

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

    /// parameters
    beta_max=10;
    gamma_max=30;
    delta_beta=0.2;
    delta_gamma=0.2;
    q_value=0.8;
    max_step=pow(10,4);
    mutation_step=0.001;
    mut_rate=0.5;

    //// generate random number
    // if you want to set seed value mamually
    // int seed = 0;
    // mt19937 gen(seed);

    // if you want to change seed value randomly
    random_device rd;
    mt19937 gen(rd());

    // array for record the result, (β,γ,(β,γ,x,y,v))
    // For Scenario 1
    vector <vector <vector <double> > > end_values1(beta_max/delta_beta+1, vector < vector <double> > (gamma_max/delta_gamma+1, vector <double>(5)));
    // For Scenario 2
    vector <vector <vector <double> > > end_values2(beta_max/delta_beta+1, vector < vector <double> > (gamma_max/delta_gamma+1, vector <double>(5)));

    ///explore the parameter space
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        beta=beta_count*delta_beta;
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            gamma=gamma_count*delta_gamma;                
            // Scenario 1
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
                w_r=get_w_resident1(x_r,y_r,v_r,beta,gamma,q_value);
                w_m=get_w_mutant1(x_r,y_r,v_r,x_m,y_m,v_m,beta,gamma,q_value);
                // if fitness of mutant>fitness of resident, invasion success and update the resident values
                if(w_r<w_m){
                    x_r=x_m;
                    y_r=y_m;
                    v_r=v_m;
                }
            }
            // record results in arrays
            end_values1[beta_count][gamma_count][0]=beta;
            end_values1[beta_count][gamma_count][1]=gamma;
            end_values1[beta_count][gamma_count][2]=x_r;
            end_values1[beta_count][gamma_count][3]=y_r;
            end_values1[beta_count][gamma_count][4]=v_r;
            // reset the initial values
            x_i=0;
            y_i=0;
            v_i=0;

            /// Scenario2
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
                w_r=get_w_resident2(x_r,y_r,v_r,beta,gamma,q_value);
                w_m=get_w_mutant2(x_r,y_r,v_r,x_m,y_m,v_m,beta,gamma,q_value);
                // if fitness of mutant>fitness of resident, invasion success and update the resident values
                if(w_r<w_m){
                    x_r=x_m;
                    y_r=y_m;
                    v_r=v_m;
                }
            }

            // record results in arrays
            end_values2[beta_count][gamma_count][0]=beta;
            end_values2[beta_count][gamma_count][1]=gamma;
            end_values2[beta_count][gamma_count][2]=x_r;
            end_values2[beta_count][gamma_count][3]=y_r;
            end_values2[beta_count][gamma_count][4]=v_r;
            // reset the initial values
            x_i=0;
            y_i=0;
            v_i=0;
        }
    }

    ////export the result
    /// Scenario 1
    ofstream fout1("result_scenario1.csv");
    fout1 << "beta,gamma,x_value,y_value,v_value\n";
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            fout1 << end_values1[beta_count][gamma_count][0] << "," 
                 << end_values1[beta_count][gamma_count][1] << ","
                 << end_values1[beta_count][gamma_count][2] << ","
                 << end_values1[beta_count][gamma_count][3] << ","
                 << end_values1[beta_count][gamma_count][4] << "\n";
        }
    }
    fout1.close();

    /// Scenario 2
    ofstream fout2("result_scenario2.csv");
    fout2 << "beta,gamma,x_value,y_value,v_value\n";
    for(beta_count=0; beta_count<=beta_max/delta_beta; beta_count++){
        for(gamma_count=0; gamma_count<=gamma_max/delta_gamma; gamma_count++){
            fout2 << end_values2[beta_count][gamma_count][0] << "," 
                 << end_values2[beta_count][gamma_count][1] << ","
                 << end_values2[beta_count][gamma_count][2] << ","
                 << end_values2[beta_count][gamma_count][3] << ","
                 << end_values2[beta_count][gamma_count][4] << "\n";
        }
    }
    fout2.close();

    //// visualize using gunuplot
    FILE *gnuplot = popen("gnuplot -persistent", "w");
    
    /// set output style(png style)
    //fprintf(gnuplot, "set terminal pngcairo size 1200,1800 enhanced font 'Arial,12'\n");
    //fprintf(gnuplot, "set output 'FigureH1.png'\n");
    /// set output style(eps style)
    fprintf(gnuplot, "set terminal epscairo size 10,12 enhanced font 'Arial,20'\n");
    fprintf(gnuplot, "set output 'FigureH1.eps'\n");

    fprintf(gnuplot, "set multiplot layout 6,2\n");
    fprintf(gnuplot, "set rmargin 6\n");

    //// Scenario1
    fprintf(gnuplot, "set label '(a) Scenario 1' at screen 0, 0.985 front font 'Arial,26'\n");

    /// Social learning
    fprintf(gnuplot, "set size 0.5,0.175\n");
    fprintf(gnuplot, "set origin 0,0.795\n");
    fprintf(gnuplot, "set label '(a-1) {/:Italic{/Symbol t}}_{SL}' at screen 0, 0.96 front font 'Arial,24'\n");
    fprintf(gnuplot, "set ylabel '{/:Italic β}' offset 2,0\n");
    fprintf(gnuplot, "set xlabel '{/:Italic γ}' offset 0,0.7\n");
    fprintf(gnuplot, "set xtics offset 0, 0.2\n");
    fprintf(gnuplot, "set ytics offset +2, 0\n");
    fprintf(gnuplot, "set tics out \n");
    fprintf(gnuplot, "set ytics scale 0.5 \n");
    
    // setting for heatmap views
    fprintf(gnuplot, "set view map\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "set palette viridis\n");
    fprintf(gnuplot, "set cbrange [0.2:0.4]\n");
    
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
    fprintf(gnuplot, "$DATAS1 << EOD\n");
    // 2)export result
    for (beta_count = 0; beta_count <= beta_max/delta_beta; beta_count++) {
        beta = beta_count * delta_beta;
        for (gamma_count = 0; gamma_count <= gamma_max/delta_gamma; gamma_count++) {
            gamma = gamma_count * delta_gamma;
            double val = (1 - end_values1[beta_count][gamma_count][2])
                       * end_values1[beta_count][gamma_count][4];
            fprintf(gnuplot, "%f ", val);
        }
        fprintf(gnuplot, "\n");
    }
    // 3)end data block
    fprintf(gnuplot, "EOD\n");
    // 4)plot data block using splot
    fprintf(gnuplot, "splot $DATAS1 matrix with image\n");
    fflush(gnuplot);


    /// Individual leaening
    fprintf(gnuplot, "set label '(a-2) {/:Italic{/Symbol t}}_{IL}' at screen 0, 0.80 front font 'Arial,24'\n");
    fprintf(gnuplot, "set size 0.5,0.175\n");
    fprintf(gnuplot, "set origin 0,0.635\n");
    fprintf(gnuplot, "set view map\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "set palette viridis\n");
    fprintf(gnuplot, "set cbrange [0:0.5]\n");
    
    /// send data in data block format
    // 1)declare data block
    fprintf(gnuplot, "$DATAI1 << EOD\n");
    // 2)export result
    for (beta_count = 0; beta_count <= beta_max/delta_beta; beta_count++) {
        beta = beta_count * delta_beta;
        for (gamma_count = 0; gamma_count <= gamma_max/delta_gamma; gamma_count++) {
            gamma = gamma_count * delta_gamma;
            double val = end_values1[beta_count][gamma_count][2]
                       * end_values1[beta_count][gamma_count][4];
            fprintf(gnuplot, "%f ", val);
        }
        fprintf(gnuplot, "\n");
    }
    // 3)end data block
    fprintf(gnuplot, "EOD\n");
    // 4)plot data block using splot
    fprintf(gnuplot, "splot $DATAI1 matrix with image\n");
    fflush(gnuplot);

    /// Exploitation
    fprintf(gnuplot, "set label '(a-3) {/:Italic{/Symbol t}}_{EX}' at screen 0, 0.64 front font 'Arial,24'\n");
    fprintf(gnuplot, "set size 0.5,0.175\n");
    fprintf(gnuplot, "set origin 0,0.475\n");
    fprintf(gnuplot, "set view map\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "set palette viridis\n");
    fprintf(gnuplot, "set cbrange [0.3:0.6]\n");

    /// send data in data block format
    // 1)declare data block
    fprintf(gnuplot, "$DATAE1 << EOD\n");
    // 2)export result
    for (beta_count = 0; beta_count <= beta_max/delta_beta; beta_count++) {
        beta = beta_count * delta_beta;
        for (gamma_count = 0; gamma_count <= gamma_max/delta_gamma; gamma_count++) {
            gamma = gamma_count * delta_gamma;
            double val = (1- end_values1[beta_count][gamma_count][4])
                       * (1- end_values1[beta_count][gamma_count][3]);
            fprintf(gnuplot, "%f ", val);
        }
        fprintf(gnuplot, "\n");
    }
    // 3)end data block
    fprintf(gnuplot, "EOD\n");
    // 4)plot data block using splot
    fprintf(gnuplot, "splot $DATAE1 matrix with image\n");
    fflush(gnuplot);

    /// Teaching
    fprintf(gnuplot, "set label '(a-4) {/:Italic{/Symbol t}}_{TE}' at screen 0, 0.48 front font 'Arial,24'\n");
    fprintf(gnuplot, "set size 0.5,0.175\n");
    fprintf(gnuplot, "set origin 0,0.315\n");
    fprintf(gnuplot, "set view map\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "set palette viridis\n");
    fprintf(gnuplot, "set cbrange [0:0.3]\n");

    /// send data in data block format
    // 1)declare data block
    fprintf(gnuplot, "$DATAT1 << EOD\n");
    // 2)export result
    for (beta_count = 0; beta_count <= beta_max/delta_beta; beta_count++) {
        beta = beta_count * delta_beta;
        for (gamma_count = 0; gamma_count <= gamma_max/delta_gamma; gamma_count++) {
            gamma = gamma_count * delta_gamma;
            double val = (1- end_values1[beta_count][gamma_count][4])
                       * end_values1[beta_count][gamma_count][3];
            fprintf(gnuplot, "%f ", val);
        }
        fprintf(gnuplot, "\n");
    }
    // 3)end data block
    fprintf(gnuplot, "EOD\n");
    // 4)plot data block using splot
    fprintf(gnuplot, "splot $DATAT1 matrix with image\n");
    fflush(gnuplot);

    /// z-value(cultural level)
    fprintf(gnuplot, "set label '(a-5) {/:Italic z}' at screen 0, 0.32 front font 'Arial,24'\n");
    fprintf(gnuplot, "set size 0.5,0.175\n");
    fprintf(gnuplot, "set origin 0,0.16\n");
    fprintf(gnuplot, "set view map\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "load 'plasma.pal'\n");  
    fprintf(gnuplot, "set pm3d map\n");
    fprintf(gnuplot, "set cbrange [0.4:1.6]\n");

    /// send data in data block format
    // 1)declare data block
    fprintf(gnuplot, "$DATAZ1 << EOD\n");
    // 2)export result
    for (beta_count = 0; beta_count <= beta_max/delta_beta; beta_count++) {
        beta = beta_count * delta_beta;
        for (gamma_count = 0; gamma_count <= gamma_max/delta_gamma; gamma_count++) {
            gamma = gamma_count * delta_gamma;
            double z_value;
            z_value=get_z_resident1(end_values1[beta_count][gamma_count][2],end_values1[beta_count][gamma_count][3],end_values1[beta_count][gamma_count][4],beta,gamma,q_value);
            fprintf(gnuplot, "%f ", z_value);
        }
        fprintf(gnuplot, "\n");
    }
    // 3)end data block
    fprintf(gnuplot, "EOD\n");
    // 4)plot data block using splot
    fprintf(gnuplot, "splot $DATAZ1 matrix with image\n");
    fflush(gnuplot);
    
    // w(fitness)
    fprintf(gnuplot, "set label '(a-6) {/:Italic w}' at screen 0, 0.165 front font 'Arial,24'\n");
    fprintf(gnuplot, "set size 0.5,0.175\n");
    fprintf(gnuplot, "set origin 0,0\n");
    fprintf(gnuplot, "set view map\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "load 'plasma.pal'\n");
    fprintf(gnuplot, "set pm3d map\n");
    fprintf(gnuplot, "set cbrange [1.02:1.04]\n");

    /// send data in data block format
    // 1)declare data block
    fprintf(gnuplot, "$DATAW1 << EOD\n");
    // 2)export result
    for (beta_count = 0; beta_count <= beta_max/delta_beta; beta_count++) {
        beta = beta_count * delta_beta;
        for (gamma_count = 0; gamma_count <= gamma_max/delta_gamma; gamma_count++) {
            gamma = gamma_count * delta_gamma;
            double w_value;
            w_value=get_w_resident1(end_values1[beta_count][gamma_count][2],end_values1[beta_count][gamma_count][3],end_values1[beta_count][gamma_count][4],beta,gamma,q_value);
            fprintf(gnuplot, "%f ", w_value);
        }
        fprintf(gnuplot, "\n");
    }
    // 3)end data block
    fprintf(gnuplot, "EOD\n");
    // 4)plot data block using splot
    fprintf(gnuplot, "splot $DATAW1 matrix with image\n");


    //// Scenario2
    /// Social learning
    fprintf(gnuplot, "set label '(b) Scenario 2' at screen 0.5, 0.985 front font 'Arial,26'\n");
    fprintf(gnuplot, "set label '(b-1) {/:Italic{/Symbol t}}_{SL}' at screen 0.5, 0.96 front font 'Arial,24'\n");
    fprintf(gnuplot, "set size 0.5,0.175\n");
    fprintf(gnuplot, "set origin 0.52,0.795\n");
    
    // setting for heatmap views
    fprintf(gnuplot, "set view map\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "set palette viridis\n");
    fprintf(gnuplot, "set cbrange [0.2:0.4]\n");
    
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
    fprintf(gnuplot, "$DATAS2 << EOD\n");
    // 2)export result
    for (beta_count = 0; beta_count <= beta_max/delta_beta; beta_count++) {
        beta = beta_count * delta_beta;
        for (gamma_count = 0; gamma_count <= gamma_max/delta_gamma; gamma_count++) {
            gamma = gamma_count * delta_gamma;
            double val = (1 - end_values2[beta_count][gamma_count][2])
                       * end_values2[beta_count][gamma_count][4];

            fprintf(gnuplot, "%f ", val);
        }
        fprintf(gnuplot, "\n");
    }
    // 3)end data block
    fprintf(gnuplot, "EOD\n");
    // 4)plot data block using splot
    fprintf(gnuplot, "splot $DATAS2 matrix with image\n");
    fflush(gnuplot);

    /// Individual leaening
    fprintf(gnuplot, "set label '(b-2) {/:Italic{/Symbol t}}_{IL}' at screen 0.5, 0.8 front font 'Arial,24'\n");
    fprintf(gnuplot, "set size 0.5,0.175\n");
    fprintf(gnuplot, "set origin 0.52,0.635\n");
    fprintf(gnuplot, "set view map\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "set palette viridis\n");
    fprintf(gnuplot, "set cbrange [0:0.5]\n");

    /// send data in data block format
    // 1)declare data block
    fprintf(gnuplot, "$DATAI2 << EOD\n");
    // 2)export result
    for (beta_count = 0; beta_count <= beta_max/delta_beta; beta_count++) {
        beta = beta_count * delta_beta;
        for (gamma_count = 0; gamma_count <= gamma_max/delta_gamma; gamma_count++) {
            gamma = gamma_count * delta_gamma;
            double val = end_values2[beta_count][gamma_count][2]
                       * end_values2[beta_count][gamma_count][4];
            fprintf(gnuplot, "%f ", val);
        }
        fprintf(gnuplot, "\n");
    }
    // 3)end data block
    fprintf(gnuplot, "EOD\n");
    // 4)plot data block using splot
    fprintf(gnuplot, "splot $DATAI2 matrix with image\n");
    fflush(gnuplot);

    /// Exploitation
    fprintf(gnuplot, "set label '(b-3) {/:Italic{/Symbol t}}_{EX}' at screen 0.5, 0.64 front font 'Arial,24'\n");
    fprintf(gnuplot, "set size 0.5,0.175\n");
    fprintf(gnuplot, "set origin 0.52,0.475\n");
    fprintf(gnuplot, "set view map\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "set palette viridis\n");
    fprintf(gnuplot, "set cbrange [0.3:0.6]\n");
    
    /// send data in data block format
    // 1)declare data block
    fprintf(gnuplot, "$DATAE2 << EOD\n");
    // 2)export result
    for (beta_count = 0; beta_count <= beta_max/delta_beta; beta_count++) {
        beta = beta_count * delta_beta;
        for (gamma_count = 0; gamma_count <= gamma_max/delta_gamma; gamma_count++) {
            gamma = gamma_count * delta_gamma;
            double val = (1- end_values2[beta_count][gamma_count][4])
                       * (1- end_values2[beta_count][gamma_count][3]);

            fprintf(gnuplot, "%f ", val);
        }
        fprintf(gnuplot, "\n");
    }
    // 3)end data block
    fprintf(gnuplot, "EOD\n");
    // 4)plot data block using splot
    fprintf(gnuplot, "splot $DATAE2 matrix with image\n");
    fflush(gnuplot);
    
    /// Teaching
    fprintf(gnuplot, "set label '(b-4) {/:Italic{/Symbol t}}_{TE}' at screen 0.5, 0.48 front font 'Arial,24'\n");
    fprintf(gnuplot, "set size 0.5,0.175\n");
    fprintf(gnuplot, "set origin 0.52,0.315\n");
    fprintf(gnuplot, "set view map\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "set palette viridis\n");
    fprintf(gnuplot, "set cbrange [0:0.3]\n");

    /// send data in data block format
    // 1)declare data block
    fprintf(gnuplot, "$DATAT2 << EOD\n");
    // 2)export result
    for (beta_count = 0; beta_count <= beta_max/delta_beta; beta_count++) {
        beta = beta_count * delta_beta;
        for (gamma_count = 0; gamma_count <= gamma_max/delta_gamma; gamma_count++) {
            gamma = gamma_count * delta_gamma;
            double val = (1- end_values2[beta_count][gamma_count][4])
                       * end_values2[beta_count][gamma_count][3];
            fprintf(gnuplot, "%f ", val);
        }
        fprintf(gnuplot, "\n");
    }
    // 3)end data block
    fprintf(gnuplot, "EOD\n");
    // 4)plot data block using splot
    fprintf(gnuplot, "splot $DATAT2 matrix with image\n");
    fflush(gnuplot);

    /// z-value(cultural level)
    fprintf(gnuplot, "set label '(b-5) {/:Italic z}' at screen 0.5, 0.32 front font 'Arial,24'\n");
    fprintf(gnuplot, "set size 0.5,0.175\n");
    fprintf(gnuplot, "set origin 0.52,0.16\n");
    fprintf(gnuplot, "set view map\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "load 'plasma.pal'\n");  
    fprintf(gnuplot, "set pm3d map\n");
    fprintf(gnuplot, "set cbrange [0.4:1.6]\n");

    /// send data in data block format
    // 1)declare data block
    fprintf(gnuplot, "$DATAZ2 << EOD\n");
    // 2)export result
    for (beta_count = 0; beta_count <= beta_max/delta_beta; beta_count++) {
        beta = beta_count * delta_beta;
        for (gamma_count = 0; gamma_count <= gamma_max/delta_gamma; gamma_count++) {
            gamma = gamma_count * delta_gamma;
            double z_value;
            z_value=get_z_resident2(end_values2[beta_count][gamma_count][2],end_values2[beta_count][gamma_count][3],end_values2[beta_count][gamma_count][4],beta,gamma,q_value);
            fprintf(gnuplot, "%f ", z_value);
        }
        fprintf(gnuplot, "\n");
    }
    // 3)end data block
    fprintf(gnuplot, "EOD\n");
    // 4)plot data block using splot
    fprintf(gnuplot, "splot $DATAZ2 matrix with image\n");
    fflush(gnuplot);
    
    /// w(fitness)
    fprintf(gnuplot, "set label '(b-6) {/:Italic w}' at screen 0.5, 0.165 front font 'Arial,24'\n");
    fprintf(gnuplot, "set size 0.5,0.175\n");
    fprintf(gnuplot, "set origin 0.52,0\n");
    fprintf(gnuplot, "set view map\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "load 'plasma.pal'\n");
    fprintf(gnuplot, "set pm3d map\n");
    fprintf(gnuplot, "set cbrange [1.02:1.04]\n");

    /// send data in data block format
    // 1)declare data block
    fprintf(gnuplot, "$DATAW2 << EOD\n");
    // 2)export result
    for (beta_count = 0; beta_count <= beta_max/delta_beta; beta_count++) {
        beta = beta_count * delta_beta;
        for (gamma_count = 0; gamma_count <= gamma_max/delta_gamma; gamma_count++) {
            gamma = gamma_count * delta_gamma;
            double w_value;
            w_value=get_w_resident2(end_values2[beta_count][gamma_count][2],end_values2[beta_count][gamma_count][3],end_values2[beta_count][gamma_count][4],beta,gamma,q_value);
            fprintf(gnuplot, "%f ", w_value);
        }
        fprintf(gnuplot, "\n");
    }
    // 3)end data block
    fprintf(gnuplot, "EOD\n");
    // 4)plot data block using splot
    fprintf(gnuplot, "splot $DATAW2 matrix with image\n");

    fprintf(gnuplot, "unset multiplot\n");
    fprintf(gnuplot, "unset output\n");
    fflush(gnuplot);
}