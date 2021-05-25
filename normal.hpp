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
        d1d Run_Evaluation_Result;
        d2d Record_point;
        d1d Record_Objective_Value;
        double max;
        double min;
    public:
        void RUN(int RUN,int ITER,int POP,int DIM,int A,int H,int pbest,const char *F,int OUTPUT_NODE_QUANTITY,const char *MODE,int EVALUATION)
        {

            Run_Result.resize(RUN);
            Run_Evaluation_Result.resize(EVALUATION/500,0);
            int r = 0;
            srand( time(NULL) );
            double START = clock();

            while(r<RUN)
            {
                srand( time(NULL) );
                INI(ITER,POP,DIM,A,H,OUTPUT_NODE_QUANTITY);

                Particle_INI(DIM,F);

                int iteration = 0;
                while(EVA < EVALUATION)
                {
                    RANK();
                    Mutation_Selection(pbest,A,H,iteration,F);
                    RANK();

                    if(MODE ==std::string("L"))//LSHADE
                    {
                        Particle = Linear_Reduction(ITER, iteration,Objective_Value.size());
                    }
                    else if(MODE ==std::string("R")){//Record
                        Record_Point(OUTPUT_NODE_QUANTITY/ITER);

                        if (iteration % 1000 == 0)
                            OUTPUT_RECORD_NODE(F,DIM,iteration);
                    }
                

                    iteration ++;
                }
                Run_Result[r] = Current_Best;
                r++;
            }
            double END = clock();
            OUTPUT_RECORD_NODE(F,DIM,ITER);
            OUT( RUN, ITER,DIM,POP,A,H,pbest,EVALUATION,START,END,F);
        }
    private:
        d2d Particle;
        d1d Objective_Value;
        d2d Archieve;
        d2d H_Table;
        d1d Objective_Rank_INDEX;

        int Archieve_Coef;
        int H_Table_Coef;
        int EVA ;
        // int NOW_POP;
        double Current_Best;
        int CURRENT_RECORD_NODE; //現在存了多少個點了
    
    private:
    void INI(int ITER,int POP,int DIM,int A,int H,int OUTPUT_NODE_QUANTITY)
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

        Record_point.clear();
        Record_point.swap(Record_point);

        Record_Objective_Value.clear();
        Record_Objective_Value.swap(Record_Objective_Value);

        Particle.resize(POP,d1d(DIM));

        Objective_Value.resize(POP);

        Archieve.resize(A,d1d(DIM));

        H_Table.resize(H,d1d(2,0.5));

        Record_point.resize(OUTPUT_NODE_QUANTITY,d1d(DIM));

        Record_Objective_Value.resize(OUTPUT_NODE_QUANTITY);

        Objective_Rank_INDEX.resize(POP);
        for(int i=0;i<POP;i++)
        {
            Objective_Rank_INDEX[i] = i;
        }

        Current_Best  = DBL_MAX;

        CURRENT_RECORD_NODE = 0;
        
        EVA = 0;

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
        else if(F ==std::string("SP"))
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
        else if (F ==std::string("B"))
        {
            for(int i=0;i<Particle.size();i++)
            {
                Bent_Cigar(DIM,i);
            }
        }
        else if (F==std::string("S"))
        {
            for(int i=0;i<Particle.size();i++)
            {
                Schaffer_F7(DIM,i);
            }
        }
        else if(F ==std::string("Z"))
        {
            for(int i=0;i<Particle.size();i++)
            {
                Zakharov(DIM,i);
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
            pairt[i].second = i; 
        } 
       
        // Sorting the pair array. 
        sort(pairt, pairt + n); 
       
        // Modifying original arrays 
        for (int i = 0; i < n; i++)  
        { 
            b[i] = pairt[i].second; 
        } 
    } 
    void RANK()//sort descending
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
           
            EVA++;
            if(EVA % 500 == 0)
            {
                Run_Evaluation_Result[(EVA-500)/500] += Current_Best;
                cout<<"# "<<EVA<<' '<<Current_Best<<endl;
            }
        }
        
        Update_Htable(H,S_Table);
       

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

         double Bent_Cigar_OBJECTIVE_VALUE(int DIM,d1d arr)
        {
            double sum1= 0;
            double sum2 = 0;
            sum1 = pow(arr[0],2);
            for(int i=1;i<DIM;i++)
            {

                sum2 += pow(arr[i],2);

            }
            double F =  sum1 +1000000*sum2;
            return F;
        }
        void Bent_Cigar(int DIM,int index)
        {
            max = 100;
            min = -100;
      
            for(int i=0;i<DIM;i++)
            {
                double a = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
                Particle[index][i] = a;
            }

            double F = Bent_Cigar_OBJECTIVE_VALUE(DIM,Particle[index]);
            Objective_Value[index] = F;
        }
        double Schaffer_F7_OBJECTIVE_VALUE(int DIM,d1d arr)
        {
            double F = 0;
            for(int i=0;i<DIM-1;i++)
            {
                double si = sqrt( pow(arr[i],2)+pow(arr[i+1],2) );
                double F1 = sqrt(si)* (sin(50*pow(si,0.2))+1);
                F += pow(F1/(DIM-1) , 2);

            }
            
            return F;
        }
        void Schaffer_F7(int DIM,int index)
        {
            max = 100;
            min = -100;
      
            for(int i=0;i<DIM;i++)
            {
                double a = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
                Particle[index][i] = a;
            }

            double F = Schaffer_F7_OBJECTIVE_VALUE(DIM,Particle[index]);
            Objective_Value[index] = F;
        }
        double Zakharov_OBJECTIVE_VALUE(int DIM,d1d arr)
        {
            double sum1 = 0;
            double sum2 = 0;
            for(int i=0;i<DIM;i++)
            {
                sum1 += pow(arr[i],2);
                sum2 += 0.5*(i+1)*arr[i];

            }
            double F =  sum1 +pow(sum2,2)+pow(sum2,4);
            return F;
        }
        void Zakharov(int DIM,int index)
        {
            max = 10;
            min = -5;
      
            for(int i=0;i<DIM;i++)
            {
                double a = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
                Particle[index][i] = a;
            }

            double F = Zakharov_OBJECTIVE_VALUE(DIM,Particle[index]);
            Objective_Value[index] = F;
        }
        double Griewank_OBJECTIVE_VALUE(int DIM,d1d arr)
        {
            double sum1 = 0;
            double sum2 = 0;
            for(int i=0;i<DIM;i++)
            {
                sum1 += pow(arr[i],2);
                sum2 *= cos( (arr[i]/sqrt(i)) );

            }
            double F =  sum1/4000 - sum2 +1;
            return F;
        }
        void Griewank(int DIM,int index)
        {
            max = 600;
            min = -600;
      
            for(int i=0;i<DIM;i++)
            {
                double a = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
                Particle[index][i] = a;
            }

            double F = Griewank_OBJECTIVE_VALUE(DIM,Particle[index]);
            Objective_Value[index] = F;
        }
        // double SchafferN2_OBJECTIVE_VALUE(int DIM,d1d arr)
        // {
        //     double sum1 = 0;
        //     double sum2 = 0;
        //     for(int i=0;i<DIM;i++)
        //     {
        //         sum1 += pow(arr[i],2);
        //         sum2 += pow(arr[i],2);

        //     }
        //     double F = 0.5+ sin(pow(sum1,2));
        //     return F;
        // }
        double Schwefel_OBJECTIVE_VALUE(int DIM,d1d arr)
        {
            double sum1 = 0;
            double sum2 = 0;
            for(int i=0;i<DIM;i++)
            {
                sum1 += arr[i]*sin( sqrt(arr[i]) ) ;

            }
            double F =  418.9829*DIM - sum1;
            return F;
        }
        void Schwefel(int DIM,int index)
        {
            max = 500;
            min = -500;
      
            for(int i=0;i<DIM;i++)
            {
                double a = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
                Particle[index][i] = a;
            }

            double F = Schwefel_OBJECTIVE_VALUE(DIM,Particle[index]);
            Objective_Value[index] = F;
        }
          //  double BOHACHEVSKY_OBJECTIVE_VALUE(int DIM,d1d arr)
        // {
        //只有2D函式還沒寫好
        //     double sum1 = 0;
        //     double sum2 = 0;
        //     for(int i=0;i<DIM;i++)
        //     {
        //         sum1 += arr[i]*sin( sqrt(arr[i]) ) ;

        //     }
        //     double F =  418.9829*DIM - sum1;
        //     return F;
        // }
        double SUM_SQUARES_OBJECTIVE_VALUE(int DIM,d1d arr)
        {
            double sum1 = 0;
            for(int i=0;i<DIM;i++)
            {
                sum1 += i*pow(arr[i],2);

            }
            double F = sum1;
            return F;
        }
        void SUM_SQUARES(int DIM,int index)
        {
            max = 10;
            min = -10;
      
            for(int i=0;i<DIM;i++)
            {
                double a = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
                Particle[index][i] = a;
            }

            double F = SUM_SQUARES_OBJECTIVE_VALUE(DIM,Particle[index]);
            Objective_Value[index] = F;
        }
        // double Booth_OBJECTIVE_VALUE(int DIM,d1d arr)
        // {
        //只有2D函式還沒寫好
        //     double sum1 = 0;
        //     for(int i=0;i<DIM;i++)
        //     {
        //         sum1 += i*pow(arr[i],2);

        //     }
        //     double F = sum1;
        //     return F;
        // }
        double POWELL_OBJECTIVE_VALUE(int DIM,d1d arr)
        {
            double sum1 = 0;
            for(int i=1;i<DIM/4;i++)
            {
                double temp = 0.0;
                temp += pow( (arr[i*4-3]+10*arr[i*4-2]) ,2);
                temp += 5*pow( (arr[i*4-1] - arr[i*4]) ,2);
                double t3 = (arr[i*4-2] - 2*arr[i*4-1]);
                temp += pow(t3,4);
                double t4 = (arr[i*4-3] + 10*arr[i*4]);
                temp += 10*pow(t4,4);
                sum1 += temp;

            }
            double F = sum1;
            return F;
        }
        void POWELL(int DIM,int index)
        {
            max = 5;
            min = -4;
      
            for(int i=0;i<DIM;i++)
            {
                double a = ((float(rand()) / float(RAND_MAX)) * (max - min)) + min;
                Particle[index][i] = a;
            }

            double F = POWELL_OBJECTIVE_VALUE(DIM,Particle[index]);
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
            else if(F ==std::string("SP"))
                SPHERE(DIM,index);  
            else if(F ==std::string("M"))
                Michalewicz(DIM,index);
            else if(F ==std::string("B"))
                Bent_Cigar(DIM,index);
            else if(F==std::string("S"))
                Schaffer_F7(DIM,index);
            else if(F ==std::string("Z"))
                Zakharov(DIM,index);
            else if(F ==std::string("G"))
                Griewank(DIM,index);
            else if(F ==std::string("SC"))
                Schwefel(DIM,index);
            else if(F ==std::string("SS"))
                SUM_SQUARES(DIM,index);
            else if(F ==std::string("P"))
                POWELL(DIM,index);
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
            else if (F == std::string("SP"))
                R = SPHERE_OBJECTIVE_VALUE(DIM,arr);
            else if(F ==std::string("M"))
                R = Michalewicz_OBJECTIVE_VALUE(DIM,arr);
            else if(F ==std::string("B"))
                R = Bent_Cigar_OBJECTIVE_VALUE(DIM,arr);
             else if(F==std::string("S"))
                R = Schaffer_F7_OBJECTIVE_VALUE(DIM,arr);
            else if(F ==std::string("Z"))
                R = Zakharov_OBJECTIVE_VALUE(DIM,arr);
            else if(F ==std::string("G"))
                R = Griewank_OBJECTIVE_VALUE(DIM,arr);
            else if(F ==std::string("SC"))
                R = Schwefel_OBJECTIVE_VALUE(DIM,arr);
            else if(F ==std::string("SS"))
                R = SUM_SQUARES_OBJECTIVE_VALUE(DIM,arr);
            else if(F ==std::string("P"))
                R = POWELL_OBJECTIVE_VALUE(DIM,arr);
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
        else if (F == std::string("SP"))
            return SPHERE_OBJECTIVE_VALUE(DIM,arr);
        else if(F ==std::string("M"))
            return  Michalewicz_OBJECTIVE_VALUE(DIM,arr);
        else if(F ==std::string("B"))
            return Bent_Cigar_OBJECTIVE_VALUE(DIM,arr);
        else if(F==std::string("S"))
            return Schaffer_F7_OBJECTIVE_VALUE(DIM,arr);
        else if(F ==std::string("Z"))
            return Zakharov_OBJECTIVE_VALUE(DIM,arr);
        else if(F ==std::string("G"))
            return Griewank_OBJECTIVE_VALUE(DIM,arr);
        else if(F ==std::string("SC"))
            return Schwefel_OBJECTIVE_VALUE(DIM,arr);
        else if(F ==std::string("SS"))
            return SUM_SQUARES_OBJECTIVE_VALUE(DIM,arr);
        else if(F ==std::string("P"))
            return  POWELL_OBJECTIVE_VALUE(DIM,arr);

    }
    
    void OUT(int run ,int iteration,int dim,int POP,int A,int H,int pbest,int evaluation,double START,double END,const char *F)
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
        for(int i=0;i<Run_Evaluation_Result.size();i++)
        {
            cout<<i*500+500<<' '<<Run_Evaluation_Result[i]/run <<endl;
        }
        string FUN;        
        if(F == std::string("A"))
             FUN = "Ackley";
        else if(F == std::string("R"))
            FUN = "Rastrigin";
        else if(F == std::string("RO"))
            FUN = "Rosenbrock";
        else if(F == std::string("SP"))
            FUN = "Sphere";  
        else if(F == std::string("M"))
            FUN = "Michalewicz";  
        else if(F ==std::string("B"))
            FUN = "Bent Cigar";
        else if(F==std::string("S"))
            FUN = "Schaffer_F7";
        else if(F ==std::string("Z"))
            FUN = "Zakharov";
        else if(F ==std::string("G"))
            FUN = "Griewank";
        else if(F ==std::string("SC"))
            FUN = "Schwefel";
        else if(F ==std::string("SS"))
            FUN = "SUM SQUARES";
        else if(F ==std::string("P"))
            FUN = "Powell";
        cout<<"# Testing Function : "<<FUN<<endl;
        cout<<"# Run : "<<run<<endl;
        cout<<"# Iteration :"<<iteration<<endl;
        cout<<"# Evaluation : "<<evaluation<<endl;
        cout<<"# DIM : "<<dim<<endl;
        cout<<"# Population Size : "<<POP<<endl;
        cout<<"# Achieve Size : "<<A<<endl;
        cout<<"# Historical Size : "<<H<<endl;
        cout<<"# Pbest : "<<pbest<<endl;
        cout<<"# Best Objective Value "<<endl<<BEST<<endl;
        cout<<"# Average Objective Value "<<endl<<AVG<<endl;
        cout<<"# Execution Time :"<<endl<<(END - START) / CLOCKS_PER_SEC<<"(s)"<<endl;
    }
    void Record_Point(int EACH_Iteration_Record)
    {
        int len = Particle.size();
        i1d index(len,0);
        for(int i=0;i<len;i++)
        {
            index[i] = i;
        }   

        i1d CHOOSE(EACH_Iteration_Record,0);
        int k = 0;
        while(k != CHOOSE.size())
        {
            int x = rand() % ( (index.size()-1) - 0 + 1) + 0;
            if(index[x]!=-1)
            {
            CHOOSE[k] = index[x];
            index[x] = -1;
            k++;
            }
        }
        for(int i=0;i<CHOOSE.size();i++)
        {
            Record_point[CURRENT_RECORD_NODE].assign(Particle[CHOOSE[i]].begin(),Particle[CHOOSE[i]].end());
            Record_Objective_Value[CURRENT_RECORD_NODE] = Objective_Value[CHOOSE[i]];
            CURRENT_RECORD_NODE++;
        }   
    }
     void OUTPUT_RECORD_NODE(const char *F,int DIM,int ITER)
        {
            fstream file;
           
            string Z = "RECORD/";
            string A = Z+F+to_string(DIM)+"_"+to_string(ITER)+".txt";
            file.open(A,ios::out);
            
            for(int i=0;i<CURRENT_RECORD_NODE;i++)
            {
                for(int j=0;j<Record_point[i].size();j++)
                {
                    file << Record_point[i][j]<<' ';
                }

                file<<Record_Objective_Value[i]<<endl;
            }
        
           
        }

    
    d2d Linear_Reduction(int ITER,int iteration,int NOW_POP)
    {
        int MAX_NFE = ITER;
        int NFE = iteration;
        int nmin = 5;
        int ninit = NOW_POP;

        int new_pop = (nmin - ninit)*NFE/MAX_NFE + ninit;
        if (new_pop < NOW_POP)
        {
                // cout<<"T"<<endl;
                // for(int i=0;i<Objective_Value.size();i++)
                // {
                //     // cout<<i<<' '<<<endl;
                //     cout<<i<<' '<<Objective_Rank_INDEX[i]<<' '<<Objective_Value[Objective_Rank_INDEX[i]]<<endl;
                // }
                
            
            d2d P;
            int DELETE = NOW_POP - new_pop;
            d1d DELETE_INDEX(DELETE,0);
            for(int i=0;i<DELETE;i++)
            {
                DELETE_INDEX[i] = Objective_Rank_INDEX[Objective_Rank_INDEX.size()-1-i];
            }
            sort(DELETE_INDEX.begin(),DELETE_INDEX.end(),greater<int>());
            for(int i=0;i<DELETE;i++)
            {
                
                Objective_Value.erase(Objective_Value.begin()+DELETE_INDEX[i]);
                Objective_Rank_INDEX.erase(Objective_Rank_INDEX.begin()+DELETE_INDEX[i]);
            }
            for(int i=0;i<Particle.size();i++)
            {
                if(i != DELETE_INDEX[i])
                {
                    d1d A;
                    A.assign(Particle[i].begin(),Particle[i].end());
                    P.push_back(A);
                }

            }
            // cout<<endl<<iteration<<' '<<Objective_Value.size()<<' '<<Objective_Rank_INDEX.size()<<' '<<P.size()<<endl;   
            return P;         
        }
        return Particle;
    }
};