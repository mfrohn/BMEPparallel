#include <mutex>
#include <queue>
#include <array>
#include "/usr/local/opt/libomp/include/omp.h"
#include <math.h>
#include <vector>
#include <chrono>
#include <thread>
#include <random>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include "xprb_cpp.h"
#include "xprs.h"

/* NAMESPACES */
using namespace std;   
using std::scientific;
using std::fixed; 
using namespace ::dashoptimization;

/* DEFINITIONS & CONSTANTS */
#define CURRENTVERSION 	   "7.14 - Golden Master"
#define CURRENTVERSIONDATE "July 10 2022"

#define MIN(X, Y) (((X) <= (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) >= (Y)) ? (X) : (Y))

#define MAX_INT    0x7FFFFFFF			//int infinity 
#define MAX_DOUBLE 0x1.FFFFFFFFFFFFFP+1023	//double infinity 
#define INF 	   MAX_DOUBLE			//infinity - it is set to the double one  
#define MAX_TIME   3600				//maximum running time to be considered 1h 
//#define MAX_TIME   -1				//no max running time
#define PRECISION  16                           //output precision
#define MAX_TOLL   0.0000000001                 //10^-10

//These constants are just used for formatting console output text
const char separator    = ' ';  
const int nameWidth     = 30;
const int numWidth      = 25;
const int numWidth2     = 12;
double TOLL=1;
double EPS=0.5;

/* Global variables */
mutex m; 					// Global locker for mutual exclusive access to data

struct EDGE{
	int i; 
	int j;
	long double weight;
	double prio;	// insertion priority
};

struct OPTIMUM_STRUCT{
	double tree_len;					//Length of the best so far tree 
	vector <struct EDGE> tree; 				//Best so far tree encoded as a list of edges
	long long OverallNumber_of_Primal_Solutions;	 	//...
	vector <vector <struct EDGE> > primal_bound_list;
	int core;
};

OPTIMUM_STRUCT Optimum;  		  			//Global Optimum

struct NODE{
	vector <struct EDGE> partial_tree;  
	double lower_bound;
	int start_edge;
	int end_edge;
	int taxon;
	bool empty;
	int prio;	// queue priority
};

struct CompareNodes{
	bool operator()(const NODE& u, const NODE& v){
		return u.prio < v.prio;
	}
};

priority_queue<NODE,vector<NODE>,CompareNodes> Queue; 

struct STATS{
	bool      outoftime;		/* Check if the algorithm is running out of time   	*/
	long long nodecounter;		/* B&B Node Counter					*/	
	long long int_sol_counter;	/* Integer solution counter    				*/
	long long LPCounter;      	/* LP calls executed 				 	*/
	long long DualCounter;      	/* Dual calls executed 					*/
	double    PardiRootLB;		/* Pardi's alpha root lower bound 			*/
	double    PardiGammaRootLB;	/* Pardi's gamma root lower bound 			*/
	double    KraftRootLB;		/* Formulation 2 root lower bound			*/
	double    ManifoldRootLB;	/* Formulation 3 root lower bound			*/
	double    LPRootLB;		/* Formulation 4 root lower bound 			*/
	double    RootLB;		/* Root LB: LP or Pardi depending on bound type 	*/
	double    RootUB;		/* Root upper bound 				 	*/
	double    GAP;		  	/* Root Gap     				   	*/
	double 	  start_time;       	/* start time (sec) of the global search		*/
	double    end_time;         	/* end time (sec) of the global search			*/
	bool	  core_runnnig;     	/* Specifies if the core is idle or running    		*/
};

int 	  boundtype=2;			/*  1: Pardi's one; 2: LP + Lagrangian; */
double  **dist=NULL;			/*  Distance Matrix Pointer 		*/ 	
STATS    *StatsArray=NULL; 	 	/*  Statistics array, one entry per core. By default, core 0 stores RootLB and UB */

double   *precomputedpow;       	/*  precomputed powers 2^-k  		*/
double   *precomputedpowRes;     	/*  precomputed powers 2^-k/res_factor  */
double    res_factor = 1.2446;
bool 	  obj_rescaling=false;
double 	 *beta=NULL;	   	    	/*  Pardi's beta_k vector		       */	
double  **beta_all=NULL; 		/*  Pardi's beta vector (for gamma bound only) */
double   *bgamma=NULL;
int 	  verbose=1;            	/*  Control how much output to write in the console */
long 	  OverallCoreNumber;		/*  Overall number of cores available in the machine */

// entropy lower bound information
double    opt_alpha=-1;			
double    opt_entropybound=-MAX_DOUBLE;
double    K=0.0;

int	  readingPattern=1;
double    bestPrimal;
bool      readtree=false;
bool 	  output_for_experiment=true;
int       torder=0,eorder=0,qrule=3,n_p=3,s_p=2;
string    namefile;
string    namefile2;
double    user_max_time = 3599;

//bool Out_of_Time(){if(user_max_time == MAX_TIME) return false; else return (omp_get_wtime()-StatsArray[0].start_time > user_max_time) ? true : false;}
bool Out_of_Time(){return (omp_get_wtime()-StatsArray[0].start_time > user_max_time) ? true : false;}
bool ThereAreIdleCores(int cores){lock_guard<mutex> lg(m); for(register unsigned int i=0; i<cores; i++) if(StatsArray[i].core_runnnig==false) return true; return false;}
bool AllIdleCores(int cores){lock_guard<mutex> lg(m); for(register unsigned int i=0; i<cores; i++) if(StatsArray[i].core_runnnig==true) return false; return true;}

/* My Classes */
#include "InputHandler.cpp"
#include "OutputHandler.cpp"
#include "Solver.cpp"
InputHandler *IH;

void Precompute_Pardi_Beta_Array(register unsigned int n){
	cout << "Precompute data for Pardi's gamma bound...";
  for(register unsigned int k=4; k<=n; k++){ 
  	beta[k]=INF;
    int counter=0;
    for(register unsigned int i=1; i<=k-2; i++)
   	  for(register unsigned int j=i+1; j<=k-1; j++){
	 			if (beta[k] > dist[i][k]+dist[j][k]-dist[i][j]) 
	 				beta[k] = dist[i][k]+dist[j][k]-dist[i][j];
	 			counter++;
  	  	beta_all[k][counter]=dist[i][k]+dist[j][k]-dist[i][j];
	 		}
	 	// sort beta_all[k] non-decreasing
  	for(int i=1; i<=counter-1; i++){
			int min_idx = i;
			for(int j=i+1; j<=counter; j++) if(beta_all[k][min_idx]>beta_all[k][j]) min_idx=j;
			double min = beta_all[k][min_idx];
			while(min_idx > i){beta_all[k][min_idx]=beta_all[k][min_idx-1]; min_idx--;}
			beta_all[k][i]=min;
		} 	
		bgamma[k]=0.0;
		for(register unsigned int i=1; i<=k-3; i++) 
			bgamma[k]+=beta_all[k][i]*precomputedpow[n-1-i];
		bgamma[k]+=beta_all[k][k-2]*precomputedpow[n-1-(k-3)];
  }
  cout<<"done!"<<endl;	
}

void GlobalSearch(int n, int core){
	Solver *CoreSolver=new Solver(core,n,boundtype);
	while(true){
		if(StatsArray[core].outoftime==false){
			NODE node=CoreSolver->RetrieveNode_and_Delete_from_Queue();
			if(node.empty==false){
				StatsArray[core].core_runnnig=true;
				CoreSolver->SetNode_and_Run(&node);
				StatsArray[core].core_runnnig=false;
			}	
			if(Queue.empty()==true && AllIdleCores(OverallCoreNumber)==true) break;
		}
		else{
			break;	
		}	
	}
	delete CoreSolver;
}

void Initializer_and_Terminator(char * const datafile, int taxa, bool DSrescaling, int Max_Cores_To_USE){
	int n;
	//Only Master allocates memory
	#pragma omp master
	{
		n=taxa; 
		IH->ReadDistanceMatrix(datafile,namefile,taxa,dist,n,readingPattern);
		if(readtree) IH->ReadPrimalSolution(datafile,namefile2,n);
		precomputedpow=new double[n+1]; 
		precomputedpowRes=new double[n+1]; 
		for(register unsigned int i=0; i<=n; i++){precomputedpow[i]=pow(2.0,i); precomputedpowRes[i]=pow(2.0,(double)i/res_factor);}
		beta=new double[n+1]; 
		bgamma=new double[n+1];
		beta_all=new double*[n+1]; for(register unsigned int i=1; i<=n; i++) beta_all[i]=new double[n*n+1];
		IH->MachineEpsilon();
		TOLL=MAX_TOLL;
		cout<<"Setting tolerance to "<<TOLL<<endl;	
		IH->StabilityChecker(n,dist);
		if(DSrescaling)IH->ActivateEntropyAnalysis=true;
		if(IH->ActivateEntropyAnalysis)IH->EntropyAnalysis(n,dist,DSrescaling);	
	}

	OverallCoreNumber=MIN(omp_get_num_procs(),Max_Cores_To_USE);
	#pragma omp master
	{	
		cout<<"Available computing cores: \x1b[92m"<<OverallCoreNumber<<"\x1b[0m"<<endl;
		StatsArray=new STATS[OverallCoreNumber]; 
		if(!readtree){
			Optimum.tree_len=MAX_DOUBLE;
			Optimum.tree.clear();}
		Optimum.OverallNumber_of_Primal_Solutions=0;
	}
	#pragma omp parallel for 
		for(register unsigned int i=0; i<OverallCoreNumber; i++){
			StatsArray[i].outoftime=false;
			StatsArray[i].nodecounter=0;
			StatsArray[i].int_sol_counter=0;
			StatsArray[i].LPCounter=0;
			StatsArray[i].DualCounter=0;
			StatsArray[i].RootLB=-MAX_DOUBLE;
			StatsArray[i].RootUB=MAX_DOUBLE;
			StatsArray[i].start_time=0;
			StatsArray[i].end_time=0;
			StatsArray[i].core_runnnig=false;
		} 
	

	stringstream ss; ss<<Optimum.tree_len; string sss=ss.str(); int nDigits = MAX(sss.size(),25); 
 	std::cout << std::fixed << std::setprecision(PRECISION)<<endl;
 	Solver *Core0Solver;  /* temp solver  */
	#pragma omp master 
	{
		StatsArray[0].start_time = omp_get_wtime();
		Core0Solver=new Solver(0,n,boundtype);
		Core0Solver->SetPrimal_and_Taxa_Order();
		bestPrimal=Optimum.tree_len;
		Precompute_Pardi_Beta_Array(n);	
		Core0Solver->SetLower_Bounds(boundtype);
		if(boundtype==1){
			cout<<"\nAs requested, the lower bound used during the global search will be Pardi's one."<<endl;
			StatsArray[0].RootLB=StatsArray[0].PardiGammaRootLB;
		} 		
		else{
			//cout<<"\nAs requested, the lower bound used during the global search will be based on linear programming."<<endl; 
			if(opt_entropybound>StatsArray[0].LPRootLB && opt_entropybound>StatsArray[0].PardiGammaRootLB){
				cout << "The entropy bound is the best lower bound at the root node." << endl;
				StatsArray[0].RootLB=opt_entropybound;
			}
			if(StatsArray[0].LPRootLB>=opt_entropybound && StatsArray[0].LPRootLB>=StatsArray[0].PardiGammaRootLB){
				cout << "The LP bound is the best lower bound at the root node." << endl;
				StatsArray[0].RootLB=StatsArray[0].LPRootLB;
			}
			if(StatsArray[0].PardiGammaRootLB>=opt_entropybound && StatsArray[0].PardiGammaRootLB>=StatsArray[0].LPRootLB){
				cout << "The Pardi Gamma bound is the best lower bound at the root node." << endl;
				StatsArray[0].RootLB=StatsArray[0].PardiGammaRootLB;
			}
		}
		cout<<endl;
		cout<<"Size of the starting Branch-&-Bound node-queue: "<<Queue.size()<<endl;
		StatsArray[0].GAP = abs(Optimum.tree_len-StatsArray[0].RootLB)/Optimum.tree_len*100;
 		cout<<"Primal upper bound: \t\t\t\t"<< left << setw(numWidth+PRECISION) << setfill(separator) << std::setprecision(PRECISION)<<Optimum.tree_len<<endl; 
		cout<<"Pardi's alpha bound at root node: \t\t"<< left << setw(numWidth+PRECISION) << setfill(separator) << std::setprecision(PRECISION)<<StatsArray[0].PardiRootLB<<endl; 
		cout<<"Pardi's gamma bound at root node: \t\t"<< left << setw(numWidth+PRECISION) << setfill(separator) << std::setprecision(PRECISION)<<StatsArray[0].PardiGammaRootLB<<endl; 
		if(boundtype!=1){
			cout<<"Manifold lower bound at root node: \t\t"  << left << setw(numWidth+PRECISION) << setfill(separator) << std::setprecision(PRECISION)<<StatsArray[0].ManifoldRootLB<<endl; 
			cout<<"Kraft lower bound at root node:  \t\t"  << left << setw(numWidth+PRECISION) << setfill(separator) << std::setprecision(PRECISION)<<StatsArray[0].KraftRootLB<<endl; 
			cout<<"LP lower bound at root node:  \t\t\t"  << left << setw(numWidth+PRECISION) << setfill(separator) << std::setprecision(PRECISION)<<StatsArray[0].LPRootLB<<endl; 
		}
		if(opt_entropybound>=0)cout<<"Entropy bound at root node:  \t\t"  << left << setw(numWidth+PRECISION) << setfill(separator) << std::setprecision(PRECISION)<<opt_entropybound<<" (alpha="<<opt_alpha<<")"<<endl; 
		//cout<<"+++ "<< StatsArray[0].PardiGammaRootLB << " & " << StatsArray[0].LPRootLB << " & " << entropy_bound <<endl;
		cout<<"Percentage gap at the root node: \t\t" << left << setfill(separator) << std::setprecision(PRECISION)<<std::fixed<<"\x1b[93m"<<StatsArray[0].GAP<<"\x1b[0m"<<endl;
		delete Core0Solver;
	}
	//TOLL=0.000000001;
	if(StatsArray[0].GAP > TOLL){
		#pragma omp master 
		{
			cout<<"\nStarting parallel global search on ";
			auto start = std::chrono::system_clock::now();
			std::time_t start_date = std::chrono::system_clock::to_time_t(start);
			cout<<std::ctime(&start_date); 
			if(user_max_time == MAX_TIME) cout<<"Setting no time limit for the global search."<<endl;
			else cout<<"Maximum runtime limit set to "<<std::setprecision(2)<<user_max_time<<" seconds."<<std::setprecision(PRECISION)<<endl<<endl;
			cout << "- ";	
		 	cout << left << setw(numWidth2)<< setfill(separator) << "# "; 
			cout << left << setw(numWidth+PRECISION) << setfill(separator) << "Primal Bound"; 
			cout << left << setw(numWidth+PRECISION) << setfill(separator) << "Gap(%)"; 
			cout << left << setw(numWidth) << setfill(separator) << "Time (sec)";
			cout << left << setw(numWidth) << setfill(separator) << "Queue Size";
			cout << left << setw(numWidth) << setfill(separator) << "NNI";
			cout <<endl;
			cout<<"* "
				<<setw(numWidth2)<< setfill(separator)<<std::fixed<<Optimum.OverallNumber_of_Primal_Solutions
				<<setw(numWidth+PRECISION) << setfill(separator)<<std::setprecision(PRECISION)<<Optimum.tree_len<<std::setprecision(PRECISION)
				<<setw(numWidth+PRECISION) << setfill(separator)<<std::fixed<< abs(Optimum.tree_len-StatsArray[0].RootLB)/Optimum.tree_len*100
				<<setprecision(8)
				<<setw(numWidth) << setfill(separator)<<std::fixed<<omp_get_wtime() - StatsArray[0].start_time
				<<setw(numWidth) << setfill(separator)<<Queue.size()
				<<setw(numWidth) << setfill(separator)<<"no"<<endl;

		}	
		#pragma omp barrier

		#pragma omp parallel shared(Optimum, StatsArray, dist, precomputedpow, beta, verbose, OverallCoreNumber) num_threads(OverallCoreNumber)
		{
			int core = omp_get_thread_num(); 
			GlobalSearch(n,core); //Running parallel Search;
		}
	}
	else{
		#pragma omp master 
		{
			cout<<"The solution is optimal."<<endl;
		}	
	}

	//Let's wait for every core to finish
	#pragma omp barrier
	#pragma omp master
	{ 
		StatsArray[0].end_time = omp_get_wtime();
		cout<<"\nGlobal search completed on ";
		auto start = std::chrono::system_clock::now();
		std::time_t start_date = std::chrono::system_clock::to_time_t(start);
		cout<<std::ctime(&start_date); 
		if ((StatsArray[0].end_time-StatsArray[0].start_time) >= user_max_time) cout<<"Attention: maximum runtime reached. The optimality of the primal bound is not guaranteed."<<endl;

		cout<<"\n+Global search statistics"<<endl;
		if ((StatsArray[0].end_time-StatsArray[0].start_time) < user_max_time){
			cout << left << setw(nameWidth) << setfill(separator) <<"Optimum: ";  
		}
		else cout << left << setw(nameWidth) << setfill(separator) <<"+Best primal bound value: ";	
		cout << left << setw(5)<< setfill(separator) <<"\x1b[92m"<<std::setprecision(PRECISION)<<Optimum.tree_len<<"\x1b[0m"<<endl;
		//cout<<"+++"<<endl;
		cout << "found on core " << Optimum.core << endl;
		cout << left << setw(nameWidth) << setfill(separator) <<"+Nodes: " << endl;
		long long OverallNodes =0; 
		for(register unsigned int i=0; i<OverallCoreNumber; i++) cout << "Core " << i << ": " << StatsArray[i].nodecounter << endl;
		for(register unsigned int i=0; i<OverallCoreNumber; i++) OverallNodes+= StatsArray[i].nodecounter; 
		cout << left << setw(numWidth+PRECISION)  << setfill(separator) <<std::setprecision(PRECISION)<<OverallNodes<<endl;
		cout << left << setw(nameWidth) << setfill(separator) <<"+Solution time (sec): "; 
		cout <<setprecision(8);
		cout << left << setw(numWidth+PRECISION)  << setfill(separator) << StatsArray[0].end_time-StatsArray[0].start_time; 
		if((StatsArray[0].end_time-StatsArray[0].start_time) > 3600) cout<<"("<<(StatsArray[0].end_time-StatsArray[0].start_time)/3600<<" hours)"<<endl;
		else cout<<endl;	
		Optimum.primal_bound_list.push_back(Optimum.tree);

		if(output_for_experiment){
			stringstream ss,ss1,ss2; 
			ss<<round(Optimum.tree_len); 
			string opt = ss.str(); 
			int dynamic_precision_opt = 15-opt.length(); 
			ss1<<round(StatsArray[0].RootLB); 
			opt = ss1.str();
			int dynamic_precision_lp = 15-opt.length();
			ss2<<OverallNodes; 
			opt = ss2.str();
			int dynamic_precision_nodes = 15-opt.length();
			cout<<"+++  "<<std::setprecision(PRECISION)<< fixed
						 <<namefile            <<right<< setw(30-namefile.length())    			   <<left<<" "
						 <<n                   <<left << setw(4)                       			   <<left<<" "
						 //<<torder+1              <<left << setw(4)                       			   <<left<<" "
						 //<<eorder+1              <<left << setw(4)                       			   <<left<<" "
						 //<<qrule+1               <<left << setw(4)                       			   <<left<<" "
						 //<<n_p                 <<left << setw(4)                       			   <<left<<" "
						 //<<s_p                 <<left << setw(4)                       			   <<left<<" "
						 <<Optimum.tree_len    <<left << setw(dynamic_precision_opt)   			   <<left<<" "
						 <<bestPrimal					 <<left << setw(dynamic_precision_lp)    			   <<left<<" "
						 <<StatsArray[0].RootLB<<left << setw(dynamic_precision_lp)    			   <<left<<" "
						 <<abs(Optimum.tree_len-StatsArray[0].RootLB)/Optimum.tree_len*100 << setw(15)<<left<<" "
						 <<OverallNodes        <<left << setw(dynamic_precision_nodes) 			   <<left<<" "
						 <<Optimum.OverallNumber_of_Primal_Solutions-1<<left<<setw(dynamic_precision_nodes)<<left<<" "
						 <<StatsArray[0].end_time-StatsArray[0].start_time
						 <<endl;
		}
		else{
			OutputHandler *OHandler = new OutputHandler(n);
			//OHandler->ComputeD3Hierarchy(Optimum.tree,namefile);
			//OHandler->OutputPrimalBoundList(Optimum.primal_bound_list,namefile);
			//OHandler->OutputSpectralAnalysis(Optimum.primal_bound_list,dist);
			delete OHandler;
		}	

		//Only Master frees memory	
 		delete[] beta;
 		delete[] bgamma;
 		for(register unsigned int i=1;i<=n;i++) delete [] beta_all[i]; delete [] beta_all;
		delete[] precomputedpow;
		delete[] precomputedpowRes;
		delete[] StatsArray;
		if(dist!=NULL){for(register unsigned int i=1;i<=n;i++) delete [] dist[i]; delete [] dist;}		
		delete IH;		
		cout<<"\nQuitting."<<endl<<endl;	
	}	
}

void InputParser(int argc, char * const argv[]){
	int taxa_to_consider=MAX_INT;
	bool DSrescaling=false;
	int cores_to_use=MAX_INT;
	IH->InputSettler(argc, argv, taxa_to_consider,DSrescaling,readtree,cores_to_use,boundtype,obj_rescaling,torder,eorder,qrule,n_p,s_p);
	Initializer_and_Terminator(argv[1],taxa_to_consider,DSrescaling,cores_to_use);
}


int main (int argc, char * const argv[]){
		// std::cout << "\x1B[2J\x1B[H"; //Clear screen
		cout<<"BMEP - Parallel solver"<<endl;
		cout<<"Version: "<<CURRENTVERSION<<" | release date: "<<CURRENTVERSIONDATE<<"."<<endl<<endl;
		cout<<"Authors: Daniele Catanzaro (1) | Martin Frohn (2) | Raffaele Pesenti (3) | Laurence A. Wolsey (1)"<<endl;
		cout<<"(1) Center for Operations Research and Econometrics (CORE), Université Catholique de Louvain (UCL);"<<endl;
		cout<<"(2) Department of Mathematics and Computer Science, Eindhoven University of Technology (TU/e);"<<endl;
		cout<<"(3) Department of Management, Università Ca Foscari, Venezia."<<endl<<endl;
		cout<<xbgetversion()<<endl;
		// cout<<"MAX Double: "<<MAX_DOUBLE<<" > "<<pow(2.0,100)<<"?";
		// if(MAX_DOUBLE > pow(2.0,100)) cout<<"yes"<<endl; else cout<<"false"<<endl; exit(0);

		IH=new InputHandler();
		if(argc<=1) IH->HowTo();
		else InputParser(argc,argv);
		return 0;	
}
/*Stat rosa pristina nomine 
  Nomina nuda tenemus */
