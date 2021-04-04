#include "normal.hpp"

int main(int argc,const char *argv[])
{
    int run = atoi(argv[1]);//just run
    int iteration = atoi(argv[2]);//just iteration
    int pop = atoi(argv[3]);
    int DIM = atoi(argv[4]);
    int A = atoi(argv[5]);
    int H = atoi(argv[6]);
    int pbest = atoi(argv[7]);
    const char *F = argv[8];
    int OUTPUT_NODE_QUANTITY = atoi(argv[9]);
    if( argc > 1 )
    {
        LSHADE lshade;
        lshade.RUN( run,iteration,pop, DIM, A, H, pbest,F,OUTPUT_NODE_QUANTITY);
    }
      
    else 
    {
        LSHADE lshade;
        const char *K ="A";
        lshade.RUN(1,2000,100,50,100,100,4,K,10000);

    }
   
    return 0;
}