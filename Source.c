// define necessary libraries

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>



// function for printing out network 1 and 2
void printNetworks(int np, int nip, int *W1, int *W2, int o1, int o2)
{
	
    printf("First Network \n");
    for (int i = 0; i < np; i++)
    {
        // Prints network 1 weights
        printf("[");
        for (int j = 0; j < nip; j++)
        {
            printf("%d", W1[i*nip+j]);
            if (j < nip-1)
                printf(", ");
        }
        printf("]\n");
    }
    printf("\nSecond Network \n");
    for (int i = 0; i < np; i++)
    {
        // Print network 2 weights
        printf("[");
        for (int j = 0; j < nip; j++)
        {
            printf("%d", W2[i*nip+j]);
            if (j < nip-1)
                printf(", ");
        }
        printf("]\n");
    }

    // network out
    printf("\nOutput 1: %d", o1);
    printf("\nOutput 2: %d\n", o2);
}



// function for generate random -1 +1 into N length array
void generateRandomInput( int *x,int N)
{
    int i=0;
    for (i = 0; i < N; i++)
    {
        x[i] = 1;
        if (rand() % 2 == 0)
            x[i] = -1;
    }
}

// change network weight according to learning algorithm
void updateWeights(int np, int bnd, int nip, int nout, int *pout, int *W, int *inp)
{
    for (int i = 0; i < np; i++) // loop for each perceptron
    {
        if (pout[i] == nout)
        {
            for (int j = 0; j < nip; j++) // loop for each inputs into concerned perceptron
            {   // calculate desired weight
                int w = (W[i*nip+j] - (pout[i] * inp[(i*nip)+j]));
                // if the calculated weight is out of border than set it as border
                if (w < -bnd)
                {
                    w = -bnd;
                }
                else if (w > bnd)
                {
                    w = bnd;
                }
                // set calculated weight
                W[i*nip+j] = w;
            }
        }
    }
}

// Method generates outputs for each perceptron in the network and then overall network output
int calculateOutput(int np, int nip, int *pout, int *W, int *inp)
{
    int nout = 1;
    // Outer loop handles num of perceptrons
    for (int i = 0; i < np; i++) // loop for each perceptron
    {
        // set initial output of perceptron
        pout[i] = 0;
        // calculate sum of weighted inputs
        for (int j = 0; j < nip; j++)  // loop for each inputs into concerned perceptron
        {
            // Getting perceptron output
            pout[i] = pout[i] + (W[(i*nip)+j] * inp[(i*nip)+j]);
        }
        // saturate calculated output as -1 or 1 according to threshold 0.

        if (pout[i] <= 0)
            pout[i] = -1;
        else
            pout[i] = 1;

        // calculate cumulative multiplication of every perceptron output in the network
        nout = nout * pout[i];
    }
    return nout;
}

// function for generate in range of [-b +b] random weight for n length array
void setRandomWeight(int *W,int n,int b)
{
    int i=0;
    for (i=0;i<n;i++)
        *(W+i)=(rand() % (2 * b + 1) - b);
}
// take necessary inputs from user
void takeInputs(int *nInput, int *nPerceptron,int *boundWeight, int *nAttacker,int *nInputPerc)
{
	printf("Number of inputs : ");
	fflush(stdout);
	scanf("%d", nInput);

    printf("Total numbers of perceptrons: ");
    fflush(stdout);
    scanf("%d", nPerceptron);

    if (*nInput % *nPerceptron !=0)
	*nInput= *nInput + *nPerceptron - (*nInput % *nPerceptron);	   
    *nInputPerc= *nInput /  *nPerceptron;   

    printf("Note, due to the input per perceptron should be integer, adjusted number of input is : %d\n",*nInput);

    printf("Positive boundaries of the weights: ");
    fflush(stdout);
    scanf("%d", boundWeight);

    printf("Number of attackers: ");
    fflush(stdout);
    scanf("%d", nAttacker);
}

int main(void)
{

	// initialize  comm size and each process id for MPI
	int comm_size;
	int proc_id;
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);


	// define necessary params these are number of inputs, positive boundary of weights
	// number of perceptron, number of inputs per perseptron and number of attackers and number of attacker per process
	int nInput, boundWeight, nPerceptron,nInputPerc, nAttacker,nAttackerProcess,sMax=50,iterMax=1000000;
	// time variable to measure the time costs
	clock_t begin, end;
	double time_taken;

	// seed number of random generator depends on time to take diiferent randm number for each time
	// but also depends on process d to take differen random number for each process
	srand (time(NULL)+proc_id);

	// take inputs from user if it is 0th process
	if (proc_id == 0)
	{
	    takeInputs(&nInput,&nPerceptron,&boundWeight, &nAttacker,&nInputPerc);
	}



	// Right now just 0th process know the given params.
	//So we need to broadcast the given params to the other processes
	MPI_Bcast(&nInput, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nPerceptron, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&boundWeight, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nAttacker, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nInputPerc, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// every processes need to calculate their attacker size
	// if given number of attacker is divisible by number of process no problem
	// but unless we need to spread the remains to the first n process
	nAttackerProcess=nAttacker/comm_size;
	if (proc_id < nAttacker % comm_size)
	    nAttackerProcess++;


	// create first and second layer weight matrix for attacked network, it has nPerceptron rows and nInputPerc column
	int *W1 = malloc(nPerceptron*nInputPerc*sizeof(int));
	int *W2 =  malloc(nPerceptron*nInputPerc*sizeof(int));
	// set initial weight by random in range of [-boundWeight +boundWeight]
	setRandomWeight(W1,nPerceptron*nInputPerc,boundWeight);
	setRandomWeight(W2,nPerceptron*nInputPerc,boundWeight);

	// create  attacker network weight matrix, it has 3 dimension nAttackerProcess,nPerceptron and nInputPerc column
	int *aW =  malloc(nAttackerProcess*nPerceptron*nInputPerc*sizeof(int));
	setRandomWeight(aW,nAttackerProcess*nPerceptron*nInputPerc,boundWeight);

	// Synchronisation iter, total iter and input array
	int iter = 0;
	int S = 0;
	// how many times attacker synchorised
	int* attackS =  calloc(nAttackerProcess, sizeof(int));
	int* inputsArray = malloc(nInput*sizeof(int));

	// network 1 and 2 outputs
	int net1out = 0;
	int net2out = 0;

	// perceptron outputs for network 1 and network 2
	int* pOutputs1 = malloc(nPerceptron*sizeof(int));
	int* pOutputs2 = malloc(nPerceptron*sizeof(int));



	// attackers' perceptron outputs and network outputs
	// there is number of attacker network output
	int *OutputsAttackers = malloc(nAttackerProcess*sizeof(int));
	// there is number of attacker times number of perceptron,  perceptron output
	int **pOutputsAttackers = malloc(nAttackerProcess*sizeof(int*));
	for (int i = 0; i < nAttackerProcess; i++)
		pOutputsAttackers[i] = malloc(nPerceptron*sizeof(int));

	// Store the starting clock time
	if (proc_id == 0)
	    begin = clock();


	while ((S < sMax) && (iter < iterMax))
	{
	    if (proc_id == 0)
	    {
	        // Generate inputs as -1 or +1
	        generateRandomInput( inputsArray,nInput);
	        // calculate network output
	        net1out = calculateOutput(nPerceptron, nInputPerc, pOutputs1, W1, inputsArray);
	        net2out = calculateOutput(nPerceptron, nInputPerc, pOutputs2, W2, inputsArray);
	    }
	    // Print initial network weights and outputs
	    if (proc_id == 0 && iter == 0)
	    {
	        printf("Initial state of networks\n");
	        printNetworks(nPerceptron, nInputPerc, W1, W2, net1out, net2out);
	    }

	    // Broadcasting the inputsArray to all other processes
	    MPI_Bcast(&(inputsArray[0]), nInput, MPI_INT, 0, MPI_COMM_WORLD);
	    MPI_Bcast(&net1out, 1, MPI_INT, 0, MPI_COMM_WORLD);
	    MPI_Bcast(&net2out, 1, MPI_INT, 0, MPI_COMM_WORLD);

	    // calculate Output for each attacker network and also its perceptrons output
	    for (int i = 0; i < nAttackerProcess; i++)
	    {
	        OutputsAttackers[i] = calculateOutput(nPerceptron, nInputPerc, pOutputsAttackers[i], &aW[i*nPerceptron*nInputPerc], inputsArray);
	    }

	    if (proc_id == 0)
	    {
	        // Check network outputs and apply learning rule if neccessary
	        if (net1out == net2out)
	        {
	            // update  1 and 2nd network weights
	            updateWeights(nPerceptron, boundWeight, nInputPerc, net1out, pOutputs1, W1, inputsArray);
	            updateWeights(nPerceptron, boundWeight, nInputPerc, net2out, pOutputs2, W2, inputsArray);
	            // count how many times both networks output is the same consecutively
	            S++;
	        }
	        else
	            S=0;
	    }

	    // ATTACKER: Check outputs for each attacker and apply learning rule if equal
	    for (int i = 0; i < nAttackerProcess; i++) // loop for each attacker
	    {
	        // if attackers output is the same with network 1
	        if (OutputsAttackers[i] == net1out)
	        {   // then both net1 and net 2 is the same apply weight cahnge
	            if (net1out == net2out)
	            {
	                updateWeights(nPerceptron, boundWeight, nInputPerc, OutputsAttackers[i], pOutputsAttackers[i], &aW[i*nPerceptron*nInputPerc], inputsArray);
	            }
	            // count how many times attacker output and net1 is the same consecutively
	            attackS[i] = attackS[i] + 1;
	        }
	        else
	        {   // attacker is not succeed so go to learning section
	            if (net1out == net2out)
	            {
	                int w = 0;
	                int wMin = 0;
	                int selectedPercp = 0;

	                for (int j = 0; j < nPerceptron; j++)
	                {
	                    for (int k = 0; k < nInputPerc; k++)
	                    {
	                        // Getting absolute minimum perceptron weight * input
	                        w = abs(aW[i*nPerceptron*nInputPerc+ j*nInputPerc+ k] * inputsArray[(j*nInputPerc)+k]);
	                        if ((w < wMin) || ((j == 0) && (k == 0)))
	                        {
	                            wMin = w;
	                            selectedPercp = j;
	                        }

	                    }
	                }

	                // Flipping the sign and assuming output of A
	                pOutputsAttackers[i][selectedPercp] = pOutputsAttackers[i][selectedPercp] * -1;
	                OutputsAttackers[i] = net1out;

	                // Applying learning rule
	                updateWeights(nPerceptron, boundWeight, nInputPerc, OutputsAttackers[i], pOutputsAttackers[i], &aW[i*nPerceptron*nInputPerc], inputsArray);
	            }

	            attackS[i] = 0;
	        }
	    }

	    iter++;

	    // Broadcast 'S' so other processes know when to break out of loop
	    MPI_Bcast(&S, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}

	if (proc_id == 0)
	{
	    // Print final state of networks and see how they are same
	    printf("\nFinal State of Network 1 and 2\n\n");
	    printNetworks(nPerceptron, nInputPerc, W1, W2, net1out, net2out);

	    // if network 1 and network 2 are synchronised
	    if (S == sMax)
	    {
	        printf("\n Network1 and Network2 are SYNCHRONISED !!\n");
	        printf("Synchronised iter: %d\n", S);
	        printf("Total iter: %d", iter);

	        // Take ending clock time and print time taken
	        end = clock();
	        time_taken = ((double)(end - begin) / CLOCKS_PER_SEC)*1000;
	        printf("\nTime Taken: %.3f ms\n", time_taken);
	    }
	    else
	    {
	        printf("\n Network1 and Network2 are NOT SYNCHRONISED !!!!\n");
	        printf("Total iter: %d\n", iter);
	    }

	}

	// wait other process
	MPI_Barrier(MPI_COMM_WORLD);

	// in every process atackred index starts 0 . but in general index every process should hve their own offset value
	// 0.process' ofset is 0, 1st process's ofset is number of process 0, 2nd process's ofset is sum of number of process 0 and 1. so on so forth.

	int offset = 0;
	for (int i=0;i<=proc_id;i++)
	{
	    int tmp=0;
	    for (int j=0;j<i;j++)
	    {
	        tmp=nAttacker/comm_size;
	        if (j < nAttacker % comm_size)
	            tmp++;
	    }
	    offset+=tmp;
	}


	// Print the id number of synchrorised attack network with network1.
	printf("---  List of syncronised attackers --------\n");
	for (int i = 0; i < nAttackerProcess; i++) 
	{

	    if (attackS[i] >= sMax)
	    {
	        printf("- Attacker %d Synchronised with Network 1 for '%d' iter \n", i+offset, attackS[i]);
	        //localSynced++;
	    }
	}

	// Releasing memory
	if (proc_id == 0)
	{
	    free(W1);
	    free(W2);
	}

	free(inputsArray);
	free(pOutputs1);
	free(pOutputs2);
	free(aW);
	free(attackS);


	MPI_Finalize();
	return 0;
}
