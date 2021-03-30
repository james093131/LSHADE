#include<stdio.h>
#include<fstream>
#include<iostream>
#include<sstream>
#include<stdlib.h>
#include <string.h>
#include<time.h>
#include<vector>
#include <algorithm>
#include<math.h>
#include<float.h>
#include <sys/stat.h>
#include <random>
#include <chrono>



using namespace std;

typedef vector<char> c1d;
typedef vector<c1d> c2d;
typedef vector<c2d> c3d;
typedef vector<c3d> c4d;
typedef vector<int> i1d;
typedef vector<i1d> i2d;
typedef vector<i2d> i3d;
typedef vector<i3d> i4d;
typedef vector<double>d1d;
typedef vector<d1d> d2d;
typedef vector<d2d> d3d;
typedef vector<d3d> d4d;


class LSHADE{
    public :
        d1d Run_Result;
        d1d Run_Iteration_Result;
        double max;
        double min;
    public:
        void RUN(int RUN,int ITER,int POP,int DIM,int A,int H,int pbest,const char *F)
        {

            Run_Result.resize(RUN);
            Run_Iteration_Result.resize(ITER,0);
            int r = 0;
            srand( time(NULL) );
            double START = clock();

            while(r<RUN)
            {
                srand( time(NULL) );
                INI(ITER,POP,DIM,A,H);

                Particle_INI(DIM,F);

                int iteration = 0;
                while(iteration < ITER)
                {
                    RANK();
                    Mutation_Selection(pbest,A,H,iteration,F);

                    iteration ++;
                }
                Run_Result[r] = Current_Best;
                r++;
            }
            double END = clock();
            OUT( RUN, ITER,DIM,START,END,F);
        }
    private:
        d2d Particle;
        d1d Objective_Value;
        d2d Archieve;
        d2d H_Table;
        d1d Objective_Rank_INDEX;

        int Archieve_Coef;
        int H_Table_Coef;
        // int NOW_POP;
        double Current_Best;
    
    private:
    void INI(int ITER,int POP,int DIM,int A,int H)
    {
        Particle.clear();
        Particle.swap((Particle));

        Objective_Value.clear();
        Objective_Value.swap(Objective_Value);

        Archieve.clear();
        Archieve.swap(Archieve);
        
        Archieve_Coef = 0;
        H_Table_Coef = 0;

        H_Table.clear();
        H_Table.swap(H_Table);

        Objective_Rank_INDEX.clear();
        Objective_Rank_INDEX.swap(Objective_Rank_INDEX);

        Particle.resize(POP,d1d(DIM));

        Objective_Value.resize(POP);

        Archieve.resize(A,d1d(DIM));

        H_Table.resize(H,d1d(2,0.5));

        Objective_Rank_INDEX.resize(POP);
        for(int i=0;i<POP;i++)
        {
            Objective_Rank_INDEX[i] = i;
        }

        Current_Best  = DBL_MAX;
        // NOW_POP = POP;
    }

   
    void Particle_INI(int DIM,const char*F)
    {
         if(F == std::string("A"))
            {
                for(int i=0;i<Particle.size();i++)
                {
                    ACKLEY(DIM,i);
                }
            }
        else if (F == std::string("R"))
        {
            for(int i=0;i<Particle.size();i++)
            {
                RASTRIGIN(DIM,i);
            }
        }
        else if(F ==std::string("RO"))
        {
            for(int i=0;i<Particle.size();i++)
            {
                ROSENBROCK(DIM,i);  
            }
        }
        else if(F ==std::string("S"))
        {
            for(int i=0;i<Particle.size();i++)
            {
                SPHERE(DIM,i);  
            }
        }
        else if(F ==std::string("M"))
        {
            for(int i=0;i<Particle.size();i++)
            {
                Michalewicz(DIM,i);
            }
        }
        
    }
    void pairsort(d1d a, d1d &b) 
    { 
        int n =a.size();
        pair<double, double> pairt[n]; 
    
        // Storing the respective array 
        // elements in pairs. 
        for (int i = 0; i < n; i++)  
        { 
            pairt[i].first = a[i]; 
            pairt[i].second = b[i]; 
        } 
    
        // Sorting the pair array. 
        sort(pairt, pairt + n); 
        
        // Modifying original arrays 
        for (int i = 0; i < n; i++)  
        { 
            b[i] = pairt[i].second; 
        } 
    } 
    void RANK()
    {
      pairsort(Objective_Value,Objective_Rank_INDEX); 
    }
    double Normal_Distribution(double mean)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine g1 (seed);            
        std::normal_distribution<double> distribution (mean,0.1);
        return distribution(g1);
    }
    double Cauchy_Distribution(double mean)
    {
        std::default_random_engine g2;
        std::cauchy_distribution<double> distribution(mean,0.1);
        return distribution(g2);
    }
    void Mutation_Selection(int pbest,int A,int H,int iter,const char*F)
    {
        d2d S_Table;
        d2d V(Particle.size(),d1d(Particle[0].size()));
        d1d CR_Table(Particle.size());
        d1d F_Table(Particle.size());

        for(int i=0;i<Particle.size();i++)
        {
            int r = rand() % ((H_Table.size()-1) - 0 + 1) + 0;

            double CR = Normal_Distribution(H_Table[r][0]);
            double F = Cauchy_Distribution(H_Table[r][1]);
            CR_Table[i]  = CR;
            F_Table[i] = F;

            for(int j=0;j<Particle[i].size();j++)
            {
                double a = ((float(rand()) / float(RAND_MAX)) * (1.0 - 0.0)) + 0.0;
                if (a < CR)
                {
                    V[i][j] = Mutation( pbest, CR, F,i,j);
                }
                  
                else 
                {
                    V[i][j] = Particle[i][j];
                }
                

                if(V[i][j] < min)
                    V[i][j] =min;
                else if(V[i][j] > max)
                    V[i][j]  = max;
            }
            
        }
        Evaluation(V,CR_Table,F_Table,A,H,iter,F);

    }

    double Mutation(int pbest,double CR,double F,int Pop,int Dim)
    {
        int X_P_best = Current_To_Pbest(pbest);

        int r1 = rand() % ( (Objective_Value.size() -1) - 0 + 1) + 0;
        while( r1 == Pop)
        {
            r1 = rand() % ( (Objective_Value.size() -1) - 0 + 1) + 0;
        }
        
        
        double  r2_bar ;
        int r2 = rand() % ( (Objective_Value.size() -1 + Archieve_Coef) - 0 + 1) + 0;
       
        while (r2 ==r1 || r2 == Pop)
        {
            r2 = rand() % ( (Objective_Value.size() -1 + Archieve_Coef) - 0 + 1) + 0;
        }
        
        if(r2 < Objective_Value.size())
        {
            r2_bar =  Particle[r2][Dim];
        }
        else{
            r2_bar = Archieve[r2-Objective_Value.size()][Dim];
        }
    
        
        
        double V = Particle[Pop][Dim] + F * (Particle[X_P_best][Dim] - Particle[Pop][Dim]) + F *(Particle[r1][Dim] - r2_bar); 
        return V;
    }
    int  Current_To_Pbest(int pbest)
    {
        int r = rand() % ( (pbest-1) - 0 + 1) + 0;
        return Objective_Rank_INDEX[r] ;
    }

    void Update_Htable(int H,d2d S_table)//next cccccode 
    {
        double S_CR_SUM = 0;
        double S_F_SUM = 0;
        double S_FIT_SUM = 0;
        for(int i=0;i<S_table.size();i++)
        {
            S_CR_SUM += S_table[i][0];
            S_F_SUM += S_table[i][1];
            S_FIT_SUM += S_table[i][2];
        }
        double s1 = 0;
        double s2 = 0;
        double f1 = 0;
        double f2 = 0;

        for(int i=0;i<S_table.size();i++)
        {
            s1 += (S_table[i][2]/S_FIT_SUM)*pow(S_table[i][0],2);
            s2 += (S_table[i][2]/S_FIT_SUM)* S_table[i][0];

            f1 += (S_table[i][2]/S_FIT_SUM)*pow(S_table[i][1],2);
            f2 += (S_table[i][2]/S_FIT_SUM)* S_table[i][1];
        }

        if(H_Table_Coef == H)
            H_Table_Coef = 0;

        H_Table[H_Table_Coef][0] = s1/s2;

        H_Table[H_Table_Coef][1] = f1/f2; 

        H_Table_Coef++;

    }

    void Evaluation(d2d V,d1d CR,d1d F,int A,int H,int iter,const char *K)
    {
        d2d S_Table;
        for(int i=0;i<V.size();i++)
        {
            double V_Objective_Value = Function_Evaluate(V[i].size(),V[i],K);
            if( V_Objective_Value < Objective_Value[i])
            {
                d1d X;
                X.push_back(CR[i]);
                X.push_back(F[i]);
                X.push_back(Objective_Value[i] - V_Objective_Value);
                S_Table.push_back(X);
                
            
                Particle[i].assign(V[i].begin(), V[i].end());
                
                Objective_Value[i] = V_Objective_Value;

                if(Current_Best > Objective_Value[i])
                {
                    Current_Best = Objective_Value[i];
                }

                if(Archieve_Coef == A)///這邊要改成random pop
                {
                    int r = rand() % ( (Archieve.size() -1) - 0 + 1) + 0;
                    Archieve[r].assign(V[i].begin(), V[i].end());

                }
                else{
                    Archieve[Archieve_Coef].assign(V[i].begin(), V[i].end());
                    Archieve_Coef ++;
                }
               
            }

        }
    
        Update_Htable(H,S_Table);
        Run_Iteration_Result[iter] += Current_Best;
        cout<<"# "<<iter<<' '<<Current_Best<<endl;
    }

    void ACKLEY(int DIM,int index) //random initial in ACKLEY Function and using RADVIZ calculate 2 dimension coordinates
    {
        max = 40.0;
        min = -40.0;
   
        for(int i=0;i<DIM;i++)
        {
            double a = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
            Particle[index][i] = a;

        }
        double F = ACKLEY_OBJECTIVE_VALUE(DIM,Particle[index]);
        Objective_Value[index] = F;
    }
  
    double ACKLEY_OBJECTIVE_VALUE(int DIM,d1d arr) //Calculate the objective value at ACKLEY Function 
    {
        double sum1= 0;
        double sum2 = 0;
        for(int i=0;i<DIM;i++)
        {

            sum1 += pow(arr[i],2);
            sum2 += cos(2*M_PI*arr[i]);

        }
        double F = -20*(exp((-0.2)*sqrt(sum1/DIM)))-exp(sum2/DIM)+20+exp(1);
        return F;
    }
    double RASTRIGIN_OBJECTIVE_VALUE(int DIM,d1d arr)
        {
            double sum1= 0;
            double sum2 = 0;
            for(int i=0;i<DIM;i++)
            {

                sum1 += pow(arr[i],2);
                sum2 += cos(2*M_PI*arr[i]);

            }
            double F =  sum1 - 10*sum2 +10*DIM;
            return F;
        }
          
        
    void RASTRIGIN(int DIM,int index) //random initial in RASTRIGIN Function and using RADVIZ calculate 2 dimension coordinates
    {
         max = 5.12;
         min = -5.12;

        for(int i=0;i<DIM;i++)
        {
            double a = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
            Particle[index][i] = a;

        }

        double F = RASTRIGIN_OBJECTIVE_VALUE(DIM,Particle[index]);
        Objective_Value[index] = F;
    }
    double ROSENBROCK_OBJECTIVE_VALUE(int DIM,d1d arr)
        {
            double sum1 = 0;
            double sum2 = 0;
            for(int i=1;i<DIM;i++)
            {
                sum1 += pow (arr[i] - pow(arr[i-1],2),2) ;
                sum2 += pow(arr[i]-1 ,2);
            // cout<<"S "<<sum1<<' '<<sum2<<endl;

            }
            double F =  100*sum1 +sum2;
            return F;
        }
        void ROSENBROCK(int DIM,int index)
        {
             max = 10;
             min = -5;
      
            for(int i=0;i<DIM;i++)
            {
                double a = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
                Particle[index][i]= a;
            }

            double F = ROSENBROCK_OBJECTIVE_VALUE(DIM,Particle[index]);
            Objective_Value[index] = F;
        }

        double SPHERE_OBJECTIVE_VALUE(int DIM,d1d arr)
        {
            double sum1= 0;
            for(int i=0;i<DIM;i++)
            {

                sum1 += pow(arr[i],2);


            }
            double F =  sum1;
            return F;
        }
          
        
        void SPHERE(int DIM,int index) //random initial in RASTRIGIN Function and using RADVIZ calculate 2 dimension coordinates
        {
             max = M_PI;
             min = 0;

            for(int i=0;i<DIM;i++)
            {
                double a = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
                Particle[index][i] = a;

            }

            double F = SPHERE_OBJECTIVE_VALUE(DIM,Particle[index]);
            Objective_Value[index] = F;
        }
        double Michalewicz_OBJECTIVE_VALUE(int DIM,d1d arr)
        {
            double sum = 0;
            for(int i=0;i<DIM;i++)
            {
                double cal1 = 0;
                double cal2 = 0;
                cal1 = sin(2*M_PI*arr[i]);

                double X = pow(arr[i],2);
                cal2 = sin(i*X/M_PI);
                cal2 = pow(cal2,20);

                sum += cal1*cal2;
            }
            double F =  -sum;
            return F;
        }
          
        
        void Michalewicz(int DIM,int index) //random initial in RASTRIGIN Function and using RADVIZ calculate 2 dimension coordinates
        {
             max = 5.12;
             min = -5.12;

            for(int i=0;i<DIM;i++)
            {
                double a = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
                Particle[index][i] = a;

            }

            double F = Michalewicz_OBJECTIVE_VALUE(DIM,Particle[index]);
            Objective_Value[index] = F;
        }

    void INI_FUNCTION(int DIM,int index,const char *F)
        {
            if(F == std::string("A"))
            {
                ACKLEY(DIM,index);
            }
            else if (F == std::string("R"))
                RASTRIGIN(DIM,index);
            else if(F ==std::string("RO"))
                ROSENBROCK(DIM,index);  
            else if(F ==std::string("S"))
                SPHERE(DIM,index);  
            else if(F ==std::string("M"))
                Michalewicz(DIM,index);

        }
        double FUNCTION(int DIM,d1d arr,const char *F)
        {
            double R = 0.0;
            if(F == std::string("A"))
            {
                R = ACKLEY_OBJECTIVE_VALUE(DIM,arr);
            }
            else if (F == std::string("R"))
                R = RASTRIGIN_OBJECTIVE_VALUE(DIM,arr);
           
            else if (F == std::string("RO"))
                R = ROSENBROCK_OBJECTIVE_VALUE(DIM,arr);
            else if (F == std::string("S"))
                R = SPHERE_OBJECTIVE_VALUE(DIM,arr);
            else if(F ==std::string("M"))
                R = Michalewicz_OBJECTIVE_VALUE(DIM,arr);
            return R;
        }    
    double Function_Evaluate(int DIM,d1d arr,const char *F)
    {
        if(F == std::string("A"))
        {
            return ACKLEY_OBJECTIVE_VALUE(DIM,arr);
        }
        else if (F == std::string("R"))
            return RASTRIGIN_OBJECTIVE_VALUE(DIM,arr);
           
        else if (F == std::string("RO"))
            return ROSENBROCK_OBJECTIVE_VALUE(DIM,arr);
        else if (F == std::string("S"))
            return SPHERE_OBJECTIVE_VALUE(DIM,arr);
        else if(F ==std::string("M"))
            return  Michalewicz_OBJECTIVE_VALUE(DIM,arr);
    }
    void OUT(int run ,int iteration,int dim,double START,double END,const char *F)
    {
        double BEST = Run_Result[0];
        double AVG = 0;
        for(int i=0;i<Run_Result.size();i++)
        {   
            AVG += Run_Result[i];
            if (Run_Result[i] < BEST)
                BEST = Run_Result[i];
        }
        AVG = AVG /run;
        for(int i=0;i<iteration;i++)
        {
            cout<<i+1<<' '<<Run_Iteration_Result[i]/run <<endl;
        }
        string FUN;
        if(F == std::string("A"))
             FUN = "Ackley";
        else if(F == std::string("R"))
            FUN = "Rastrigin";
        else if(F == std::string("RO"))
            FUN = "Rosenbrock";
        else if(F == std::string("S"))
            FUN = "Sphere";  
        else if(F == std::string("M"))
            FUN = "Michalewicz";  

        cout<<"# Testing Function : "<<FUN<<endl;
        cout<<"# Run : "<<run<<endl;
        cout<<"# Iteration :"<<iteration<<endl;
        cout<<"# DIM  "<<dim<<endl;
        cout<<"# Best Objective Value "<<BEST<<endl;
        cout<<"# Average Objective Value "<<AVG<<endl;
        cout<<"# Execution Time :"<<(END - START) / CLOCKS_PER_SEC<<"(s)"<<endl;
    }

      // void Linear_Reduction(int ITER,int iteration,int NOW_POP)
    // {
    //     int MAX_NFE = ITER;
    //     int NFE = iteration;
    //     int nmin = 5;
    //     int ninit = NOW_POP;

    //     int new_pop = (nmin - ninit)/MAX_NFE*NFE + ninit;

    //     if (new_pop < NOW_POP)
    //     {
    //         int DELETE = NOW_POP - new_pop;
    //         for(int i=0;i<DELETE;i++)
    //         {

    //         }

    //     }
    // }
};