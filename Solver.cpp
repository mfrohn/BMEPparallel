class Solver{
private:
	/* Engine variables */
	int Core;
	int n;
	vector <struct EDGE> tree;   
	vector <struct EDGE> *old_tree; 
	int **Tau; 
	double **Wtau;
	double **M; 
	int *mem;	
	double len;
	int* alpha_sigma;
	int *pi;
	/* **************** */
	XPRBprob problem;							/*          XPRESS Problem			 */
	XPRBvar ***path;							/*            XPRESS vars			 */
	XPRBctr *Kraft,**Unicity,Third;				/*      XPRESS Core Constraints		 */
	// Manifold version
	XPRBprob problem2;							/*          XPRESS Problem			 */
	XPRBvar ***path2;							/*            XPRESS vars			 */
	XPRBctr **Unicity2,Third2;
	// int ***DistTaxon;							/* Topological distances when taxon is inserted */
	double **DualValuesAtLeaf;                  /* This Matrix store at row i the dual values (lambda and mu) relative to taxon i */
	void LoadProblem();
	void LoadProblem2();
	void Explore(int, int, int);
	void insertPartition(int, double, int);
	void DFS(int, double);
	void Add(int, int);
	void Restore(int);
	void RandomizeTaxaOrder();
	static bool compare_edge_prio(EDGE, EDGE);
	int LBEdgeOrder(int);
	void RandomizeEdgeOrder(int);
	void RestoreEdgeOrder(int);
	void ComputeTopologicalMatrix(int);
	void ComputeWeightedOptimalTopologicalMatrix();
	void PhylogeneticDiversityTaxaOrder();
	double LagrangianDualLowerBound(int);
	double FastLinearProgrammingLowerBound(int);
	double EntropyBound(int, int, double* &);
	void EntropyAnalysis(int);
	double KraftLowerBound();
	double PardiLowerBound(int);
	double PardiGammaLowerBound(int);
	//double LinLog(int);
	bool Pruned(int, int, double);
	bool Push(int, int, int);
	void UpdateOptimum();
	void UpdateOptimum2(bool);
	void SetStartingTree();
	void ShortcutTree_and_reorder_Taxa();
	void PrintTree();
	// void PrintM();
	void NNISwap(int, int, int, int);
	void NNI(int,bool);

public:
	NODE RetrieveNode_and_Delete_from_Queue();
	void SetNode_and_Run(NODE *node);
	void SetPrimal_and_Taxa_Order();
	void SetLower_Bounds(int);
	// variables for PrimalDualLB
	int n_rows_con_eq, n_cols_con_eq;
	double **constraints_eq;
	double *dual_vars;
	double *old_dual_vars;
	double ***primal_vars;
	double ***dual_gaps;
	double *min_d;
	void PrimalDualLB();
	void NJtree();
	void SWAA();
	void SWAA2();
	// void SetStartingQueue(int taxon,int edge_cardinality);
	Solver(){}
	Solver(int CoreNum,int taxa,int boundtype){
		Core=CoreNum; 
		n=taxa;
		old_tree=new vector <struct EDGE>[n+1];
		Tau=new int*[n+1];  for(register unsigned int i=0; i<=n;  i++)   Tau[i]=new int[n+1];
	    Wtau=new double*[n+1];  for(register unsigned int i=0; i<=n;  i++)   Wtau[i]=new double[n+1];
	    M=new double*[n-1];  for(register unsigned int i=0; i<=n-2;  i++) M[i]=new double[7];  
	    mem=new int[n-1];
	    if(boundtype!=1){
	   		path=new XPRBvar**[n+1]; for(register unsigned int i=1; i<=n-1; i++) {path[i]=new  XPRBvar*[n+1]; for(register unsigned int j=i+1; j<=n; j++) path[i][j]=new  XPRBvar[n];}
	    	path2=new XPRBvar**[n+1]; for(register unsigned int i=1; i<=n-1; i++) {path2[i]=new  XPRBvar*[n+1]; for(register unsigned int j=i+1; j<=n; j++) path2[i][j]=new  XPRBvar[n];}
	    	Kraft=new XPRBctr[n+1]; 
	    	Unicity=new XPRBctr*[n+1]; for(register unsigned int i=1; i<=n; i++) Unicity[i]=new XPRBctr[n+1];
	    	Unicity2=new XPRBctr*[n+1]; for(register unsigned int i=1; i<=n; i++) Unicity2[i]=new XPRBctr[n+1];
	    	DualValuesAtLeaf=new double*[n+1]; 
	    	for(register unsigned int i=1; i<=n; i++) DualValuesAtLeaf[i]=new double[n+1];
	    	n_rows_con_eq=n*(n-1)/2+n+1; n_cols_con_eq=(n-2)*n*(n-1)/2;
			constraints_eq=new double*[n_rows_con_eq+1]; for(register unsigned int i=1; i<=n_rows_con_eq; i++) {constraints_eq[i]=new double[n_cols_con_eq+1]; for(register unsigned int j=1; j<=n_cols_con_eq; j++) constraints_eq[i][j]=0;}
			dual_vars=new double[n_rows_con_eq+1]; min_d=new double[n+1];
			old_dual_vars=new double[n_rows_con_eq+1]; min_d=new double[n+1];
			dual_gaps=new double**[n+1]; for(register unsigned int i=1; i<=n; i++) {dual_gaps[i]=new double*[n+1]; for(register unsigned int j=i+1; j<=n; j++) dual_gaps[i][j]=new double[n-1];} 
			primal_vars=new double**[n+1]; for(register unsigned int i=1; i<=n; i++) {primal_vars[i]=new double*[n+1]; for(register unsigned int j=i+1; j<=n; j++) primal_vars[i][j]=new double[n-1];} 
			LoadProblem(); 
			LoadProblem2();
		}
		SetStartingTree(); 
		ComputeTopologicalMatrix(3); 
	}		
	~Solver(){
	  tree.clear();
	  for(register unsigned int i=0; i<=n; i++) old_tree[i].clear();
	  delete [] old_tree;
	  for(register unsigned int i=0; i<=n; i++)    delete [] Tau[i]; delete [] Tau; 
	  for(register unsigned int i=0; i<=n; i++)    delete [] Wtau[i]; delete [] Wtau; 
	  for(register unsigned int i=0; i<=n-2; i++)  delete [] M[i]; delete [] M; 
	  delete [] mem; 
	  if(boundtype!=1){
	  	for(register unsigned int i=1; i<=n-1; i++) {for(register unsigned int j=i+1; j<=n; j++) delete [] path[i][j];  delete [] path[i];}  delete [] path;
	  	for(register unsigned int i=1; i<=n-1; i++) delete [] path2[i];  delete [] path2;
	  	delete[] Kraft; 
	  	for(register unsigned int i=1; i<=n; i++) delete[] Unicity[i]; delete[] Unicity;
	  	for(register unsigned int i=1; i<=n; i++) delete[] Unicity2[i]; delete[] Unicity2;
	  	for(register unsigned int i=1; i<=n; i++) delete[] DualValuesAtLeaf[i]; delete[] DualValuesAtLeaf; 
	  	for(register unsigned int i=1; i<=n_rows_con_eq; i++) delete[] constraints_eq[i]; delete[] constraints_eq;
	  	delete[] dual_vars; delete[] old_dual_vars; delete[] min_d;
	  	for(register unsigned int i=1; i<=n; i++) {for(register unsigned int j=i+1; j<=n; j++) delete [] dual_gaps[i][j];  delete [] dual_gaps[i];}  delete [] dual_gaps;
	  	for(register unsigned int i=1; i<=n; i++) {for(register unsigned int j=i+1; j<=n; j++) delete [] primal_vars[i][j];  delete [] primal_vars[i];}  delete [] primal_vars;
	  }
	  // for(register unsigned int i=0; i<=n; i++){for(register unsigned int j=0; j<=n; j++) delete [] DistTaxon[i][j]; delete [] DistTaxon[i];} delete [] DistTaxon; 
	 }
};

void Solver::PrintTree(){
	cout<<"Core: "<<Core<<" printing tree..."<<endl;
	for(register unsigned int i=0; i<(int)tree.size(); i++) cout<<tree[i].i<<"\t"<<tree[i].j<<endl;
}


void Solver::SetStartingTree(){
	struct EDGE x; x.prio=-1;
	x.i=1; x.j=n+1; tree.push_back(x);
	x.i=2; x.j=n+1; tree.push_back(x);
	x.i=3; x.j=n+1; tree.push_back(x);
	if(eorder==1){LBEdgeOrder(3);}	
}

void Solver::Add(int taxon, int e){
	struct EDGE edge1, edge2, edge3; 	
	edge1.i=taxon;      edge1.j=n+taxon-2;	edge1.prio=-1;
	edge2.i=tree[e].i;  edge2.j=n+taxon-2;	edge2.prio=-1;
	edge3.i=tree[e].j;  edge3.j=n+taxon-2;	edge3.prio=-1;

	tree.erase(tree.begin()+e);
	tree.push_back(edge1);
	tree.push_back(edge2);
	tree.push_back(edge3);
}

void Solver::Restore(int e){
	struct EDGE edge;
	long tsize=tree.size();
	edge.i=tree[tsize-2].i;
	edge.j=tree[tsize-1].i;

	tree.pop_back(); 
	tree.pop_back(); 
	tree.pop_back();
	tree.insert(tree.begin()+e,edge);
}

void Solver::ComputeTopologicalMatrix(int taxon){
	//M[i][1],M[i][2],M[i][3] store the nodes adjacent to the internal node i+n
	int tsize=tree.size();
	for(int i=0; i<=n-2;i++) mem[i]=1;
	for(int e=0; e<tsize; e++){
		int temp1=tree[e].i; int temp2=tree[e].j;	
		if (temp1 > n) {M[temp1-n][mem[temp1-n]]=temp2; mem[temp1-n]++;}
		if (temp2 > n) {M[temp2-n][mem[temp2-n]]=temp1; mem[temp2-n]++;}
	}
	register int CurrentNode,VisitedNode; 
	len=0;
	Tau[0][0] = 0; // Modified 2/6/2009
	for(int j=0; j<tsize; j++){
		if(tree[j].i<=taxon){
			int i = tree[j].i;
			int father=tree[j].j; 
			
			Tau[i][i]=0; // Modified 2/6/2009 - rimodificato il 29/9/2020
			Tau[i][0]=0; // modificato il 29/9/2020

			M[father-n][4]=1; M[father-n][5]=i; M[father-n][6]=0; 
			CurrentNode=father; 
			while(true){
				if (M[CurrentNode-n][6]<3){
					M[CurrentNode-n][6]++; VisitedNode=M[CurrentNode-n][(int)M[CurrentNode-n][6]];
					if (VisitedNode != M[CurrentNode-n][5]){
						if (VisitedNode > n){
							M[VisitedNode-n][4]=M[CurrentNode-n][4]+1;
							M[VisitedNode-n][5]=CurrentNode;
							M[VisitedNode-n][6]=0;
							CurrentNode=VisitedNode;
						}
						else{    
							 // Modified 2/6/2009 - rimodificato il 29/9/2020
							 	Tau[i][VisitedNode]=M[CurrentNode-n][4]+1;	
							 	if(Tau[i][VisitedNode]>Tau[i][0]){
							 		Tau[i][0] = Tau[i][VisitedNode];        
							 		if (Tau[i][0] > Tau[0][0]) Tau[0][0] = Tau[i][0];
							 	} 
							 // End Modified 2/6/2009 - 29/9/2020
							 if(obj_rescaling) len+=dist[i][VisitedNode]*precomputedpow[n]/precomputedpowRes[Tau[i][VisitedNode]];
							 else len+=dist[i][VisitedNode]*precomputedpow[n-Tau[i][VisitedNode]]; 
							  	
						}
					}
				}
				else{
					if(CurrentNode == father) break;
					else CurrentNode=M[CurrentNode-n][5];
				}				
			}
      	}
    }  	
}

void Solver::ComputeWeightedOptimalTopologicalMatrix(){
	cout<<"Calculate the distance between taxa given by the phylogeny...";
	//M[i][1],M[i][2],M[i][3] store the nodes adjacent to the internal node i+n
	double** edgeweights = new double*[2*n-1]; for(register unsigned int i=0; i<=2*n-2;i++){edgeweights[i]=new double[4]; for(register unsigned int j=1; j<=3; j++) edgeweights[i][j]=-1;}
	for(int i=0; i<=n-2;i++) mem[i]=1;
	int tsize=Optimum.tree.size();
	for(int e=0; e<tsize; e++){
		int temp1=Optimum.tree[e].i; int temp2=Optimum.tree[e].j;	
		if (temp1 > n) {M[temp1-n][mem[temp1-n]]=temp2; 
			edgeweights[temp1][mem[temp1-n]]=Optimum.tree[e].weight;
			mem[temp1-n]++;}
		if (temp2 > n) {M[temp2-n][mem[temp2-n]]=temp1; 
			edgeweights[temp2][mem[temp2-n]]=Optimum.tree[e].weight;
			mem[temp2-n]++;}
		if (temp1 <= n) edgeweights[temp1][1]=Optimum.tree[e].weight;
		if (temp2 <= n) edgeweights[temp2][1]=Optimum.tree[e].weight;
	}
	//for(register unsigned int i=1; i<=2*n-2;  i++)
	//	cout << edgeweights[i][1] << " " << edgeweights[i][2] << " " << edgeweights[i][3] << endl;
	register int CurrentNode,VisitedNode;
	for(int e=0; e<tsize; e++){
		bool process=false; int i, father;
		if(Optimum.tree[e].i<=n){
			i = Optimum.tree[e].i;
			father = Optimum.tree[e].j;
			process=true;}
		else if(Optimum.tree[e].j<=n){
			i = Optimum.tree[e].j;
			father = Optimum.tree[e].i;
			process=true;}
		if(process){ 
			Wtau[i][i] = 0;
			/* M[internal-n][4] = path-length from taxon i to internal
			   M[internal-n][5] = predecessor of internal
			   M[internal-n][6] = neighbours of internal visited */
			M[father-n][4]=edgeweights[i][1]; M[father-n][5]=i; M[father-n][6]=0; 
			CurrentNode=father; 
			while(true){
				if (M[CurrentNode-n][6]<3){
					M[CurrentNode-n][6]++; VisitedNode=M[CurrentNode-n][(int)M[CurrentNode-n][6]];
					if (VisitedNode != M[CurrentNode-n][5]){
						if (VisitedNode > n){ // if neighbour of CurrentNode is internal
							//cout << "set path length from " << i << " to " << VisitedNode << " = " << M[CurrentNode-n][4] << " + " << edgeweights[CurrentNode][(int)M[CurrentNode-n][6]] << endl;
							M[VisitedNode-n][4]=M[CurrentNode-n][4]+edgeweights[CurrentNode][(int)M[CurrentNode-n][6]];
							M[VisitedNode-n][5]=CurrentNode;
							M[VisitedNode-n][6]=0;
							CurrentNode=VisitedNode;
						} else // neighbour of CurrentNode is external
							//cout << "set path length from " << i << " to " << VisitedNode << " = " << M[CurrentNode-n][4] << " + " << edgeweights[VisitedNode][1] << endl;
							Wtau[i][VisitedNode]=M[CurrentNode-n][4]+edgeweights[VisitedNode][1];	
					}
				} else{ // continue recursion
					if(CurrentNode == father) break;
					else CurrentNode=M[CurrentNode-n][5];
				}				
			}
      	}
    }
    //for(int i=1; i<=n;i++){for(int j=1; j<=n;j++) cout << Wtau[i][j] << " ";cout << endl;}
    for(register unsigned int i=0; i<=2*n-2; i++) delete [] edgeweights[i]; delete [] edgeweights;
    cout<<"done!"<<endl;  	
}

void Solver::PhylogeneticDiversityTaxaOrder(){
	/*Assumption: ComputeWeightedOptimalTopologicalMatrix() was already applied, i.e.,
	 			  Wtau contains the estimated evolutionary distance between pairs of taxa
	              (of the current optimal primal solution) */
	cout<<"Calculate the Phylogenetic Diversity Taxa Order...";
	double **new_dist=new double*[n+1]; for(register unsigned int i=0; i<=n;i++) new_dist[i] = new double[n+1];
	pi=new int[n+1]; for(register unsigned int i=1; i<=n;i++) pi[i]=i;
	// find the pair of taxa at maximum estimated evolutionary distance (and put them as 1st and 2nd)
	double max=0.0; int max_i,max_j; for(int i=1;i<=n;i++)for(int j=i+1;j<=n;j++) if(Wtau[i][j]>max){max=Wtau[i][j];max_i=i;max_j=j;}
	int dummy=pi[1]; pi[1]=pi[max_i]; pi[max_i]=dummy;
	dummy=pi[2]; pi[2]=pi[max_j]; pi[max_j]=dummy;
	// insert the remaining taxa into new_tau in the appropriate order
	for(int taxon=3;taxon<=n-1;taxon++){
		double max_score=0.0; int max_score_taxon_idx;
		// calculate the accumlated estimated evolutionary distance of all taxa with the taxa inserted in new_tau
		for(int i=taxon;i<=n;i++){
			double score=0.0; for(int j=1;j<=n;j++) score+=Wtau[pi[i]][pi[j]];
			if(score>max_score){max_score=score; max_score_taxon_idx=i;}
		}
		// insert the taxon with maximum accumlated estimated evolutionary distance into new_tau
		int dummy=pi[max_score_taxon_idx]; pi[max_score_taxon_idx]=pi[taxon]; pi[taxon]=dummy;
	}
	for(register unsigned int i=1; i<=n;i++)
		for(register unsigned int j=1; j<=n;j++)
			new_dist[i][j] = dist[pi[i]][pi[j]];
	for(register unsigned int i=0; i<=n; i++) delete [] dist[i]; delete [] dist;
	dist=new_dist;
	cout<<"done!"<<endl;
	cout <<"phylogenetic diversity order:"<<endl;
}

void Solver::LoadProblem(){
	problem.setMsgLevel(0);
 	XPRSsetintcontrol(problem.getXPRSprob(), XPRS_PRESOLVE,    0);
    XPRSsetintcontrol(problem.getXPRSprob(), XPRS_CUTSTRATEGY, 0); 
	XPRSsetintcontrol(problem.getXPRSprob(), XPRS_HEURSTRATEGY,0); 
	XPRSsetintcontrol(problem.getXPRSprob(), XPRS_THREADS,  1);
	XPRSsetintcontrol(problem.getXPRSprob(), XPRS_LPTHREADS,1);
	XPRSsetintcontrol(problem.getXPRSprob(), XPRS_BARTHREADS,1);
	//Path variables only for i<j
	for(register unsigned int i=1; i<=n-1; i++){
		for(register unsigned int j=i+1; j<=n; j++){
			for(register unsigned int k=2; k<=n-1; k++){
				path[i][j][k]=problem.newVar("path",XPRB_BV);
			}
		}
	}	
	//Obj function
	XPRBlinExp Obj; 
	for(register unsigned int i=1; i<=n-1; i++){
		for(register unsigned int j=i+1; j<=n; j++){
			for(register unsigned int k=2; k<=n-1; k++){
				double coeff; 
				if(obj_rescaling) coeff = 2*dist[i][j]*precomputedpow[n]/precomputedpowRes[k];
				else coeff = 2*dist[i][j]*precomputedpow[n-k];				
				Obj+=path[i][j][k]*coeff;
			}
		}
	} 
	problem.setObj(Obj); 
	//Convexity constraint
	for(register unsigned int i=1; i<=n; i++){
		for(register unsigned int j=i+1; j<=n; j++){
			XPRBlinExp C2;
			for(register unsigned int k=2; k<=n-1; k++)  C2+=path[i][j][k]; 
			Unicity[i][j]=problem.newCtr(C2 == 1);
		}
		for(register unsigned int j=1; j<=i-1; j++){
			XPRBlinExp C2;
			for(register unsigned int k=2; k<=n-1; k++)  C2+=path[j][i][k]; 
			Unicity[j][i]=problem.newCtr(C2 == 1);
		}
	}

	//Kraft equality
	for(register unsigned int i=1; i<=n; i++){
		XPRBlinExp C1; 
		for(register unsigned int j=1; j<=i-1; j++) for(register unsigned int k=2; k<=n-1; k++) C1+=path[j][i][k]*precomputedpow[n-k]; 
		for(register unsigned int j=i+1; j<=n; j++) for(register unsigned int k=2; k<=n-1; k++) C1+=path[i][j][k]*precomputedpow[n-k]; 
		Kraft[i]=problem.newCtr(C1 == precomputedpow[n-1]);
	}

	//Third constraint
	XPRBlinExp TC3;
	for(register unsigned int i=1; i<=n-1; i++){ 
		for(register unsigned int j=i+1; j<=n; j++){ 
			for(register unsigned int k=2; k<=n-1; k++){
			   double coeff = 2*k*precomputedpow[n-k];  
			   TC3+=path[i][j][k]*coeff;
			}
		}
	}		
	Third=problem.newCtr(TC3 == (2*n-3)*precomputedpow[n]);
}


void Solver::LoadProblem2(){
	problem2.setMsgLevel(0);
 	XPRSsetintcontrol(problem2.getXPRSprob(), XPRS_PRESOLVE,    0);
    XPRSsetintcontrol(problem2.getXPRSprob(), XPRS_CUTSTRATEGY, 0); 
	XPRSsetintcontrol(problem2.getXPRSprob(), XPRS_HEURSTRATEGY,0); 
	XPRSsetintcontrol(problem2.getXPRSprob(), XPRS_THREADS,  1);
	XPRSsetintcontrol(problem2.getXPRSprob(), XPRS_LPTHREADS,1);
	XPRSsetintcontrol(problem2.getXPRSprob(), XPRS_BARTHREADS,1);
	//Path variables only for i<j
	for(register unsigned int i=1; i<=n-1; i++){
		for(register unsigned int j=i+1; j<=n; j++){
			for(register unsigned int k=2; k<=n-1; k++){
				path2[i][j][k]=problem2.newVar("path2",XPRB_PL);
			}
		}
	}
	//Obj function
	XPRBlinExp Obj; 
	for(register unsigned int i=1; i<=n-1; i++){
		for(register unsigned int j=i+1; j<=n; j++){
			for(register unsigned int k=2; k<=n-1; k++){
				Obj+=path2[i][j][k]*2*dist[i][j]*precomputedpow[n-k];
			}
		}
	} 
	problem2.setObj(Obj); 
	//Convexity constraint
	for(register unsigned int i=1; i<=n; i++){
		for(register unsigned int j=i+1; j<=n; j++){
			XPRBlinExp C2;
			for(register unsigned int k=2; k<=n-1; k++)  C2+=path2[i][j][k]; 
			Unicity2[i][j]=problem2.newCtr(C2 == 1);
		}
		for(register unsigned int j=1; j<=i-1; j++){
			XPRBlinExp C2;
			for(register unsigned int k=2; k<=n-1; k++)  C2+=path2[j][i][k]; 
			Unicity2[j][i]=problem2.newCtr(C2 == 1);
		}
	}

	//Third constraint
	XPRBlinExp TC3;
	for(register unsigned int i=1; i<=n-1; i++){ 
		for(register unsigned int j=i+1; j<=n; j++){ 
			for(register unsigned int k=2; k<=n-1; k++){
			   double coeff = 2*k*precomputedpow[n-k];  
			   TC3+=path2[i][j][k]*coeff;
			}
		}
	}		
	Third2=problem2.newCtr(TC3 == (2*n-3)*precomputedpow[n]);
}

void Solver::PrimalDualLB(){
	// Kraft constraints
	int start=0;
	for(int h=1; h<=n-1; h++){
		for(int i=1; i<=n; i++)
			for(int j=1+start; j<=(n-2)*(n-h)+start; j++)
				if(i==h) constraints_eq[i][j] = 1;
				else if((i>h)&&(j>=(i-h-1)*(n-2)+1+start)&&(j<=(i-h)*(n-2)+start)) constraints_eq[i][j] = 1;
		start+=(n-2)*(n-h);}
	// phylogenetic manifold constraint
	for(int j=1; j<=n_cols_con_eq; j++) constraints_eq[n+1][j] = 1;
	// convexity constraints
	for(int i=2+n; i<=n*(n-1)/2+n+1; i++)
		for(double j=1; j<=n_cols_con_eq; j++){
			int idx = ceil(j/(n-2));
			if(i-n-1==idx) constraints_eq[i][(int)j] = 1;}
	// set proper coefficients
	for(double j=1.0; j<=n_cols_con_eq; j++){
		int ell = j+1-floor((j-1)/(n-2))*(n-2);
		for(int i=1; i<=n; i++) constraints_eq[i][(int)j] *= pow(2.0,1-ell);
		constraints_eq[n+1][(int)j] *= ell*pow(2.0,-ell)/(2*n-3);}
	// Print constraints
	/*for(int i=1; i<=n_rows_con_eq; i++){
		for(int j=1; j<=n_cols_con_eq; j++) cout << setprecision(2) << constraints_eq[i][j] << " ";
		cout << endl;
	}*/
	// INITIALIZATION
	// dual_vars[1..n] = mu, dual_vars[n+1] = nu, dual_vars[n+2...n_rows_con_eq] = lambda
	cout << "initialize dual variables:" << endl;
	//for(int i=1; i<=n-1; i++)for(int j=i+1; j<=n; j++) dist[i][j] *= pow(10,6);
	// set nu
	dual_vars[n+1] = INF;
	for(int i=1; i<=n-1; i++)for(int j=i+1; j<=n; j++)
		if(dist[i][j]*(2*n-3)/(n-1)<dual_vars[n+1]) dual_vars[n+1] = dist[i][j]*(2*n-3)/(n-1);
	dual_vars[n+1] = 0;
	cout << "nu:" << endl << dual_vars[n+1] << endl;
	// set mu
	cout << "mu:" << endl;
	for(int i_hat=1; i_hat<=n; i_hat++){
		dual_vars[i_hat]=INF; 
		for(int i=1; i<=n-1; i++)
			for(int j=i+1; j<=n; j++)
				for(int ell=2; ell<=n-1; ell++)
					if((i==i_hat)&&(dist[i][j]-ell*dual_vars[n+1]/(2*n-3) < dual_vars[i_hat])) 
						dual_vars[i_hat]=dist[i][j]-ell*dual_vars[n+1]/(2*n-3);
					else if ((j==i_hat)&&(dist[i][j]-ell*dual_vars[n+1]/(2*n-3)-dual_vars[i] < dual_vars[i_hat])) 
						dual_vars[i_hat]=dist[i][j]-ell*dual_vars[n+1]/(2*n-3)-dual_vars[i];
		//dual_vars[i_hat]=0;
		cout << dual_vars[i_hat] << endl;
	}
	// set lambda
	cout << "lambda:" << endl;
	int sum=n+1; for(int i=1; i<=n-1; i++){ for(int j=i+1; j<=n; j++){
			dual_vars[sum+j-i] = (dist[i][j]-(n-1)*dual_vars[n+1]/(2*n-3)-dual_vars[i]-dual_vars[j])*pow(2.0,2-n);
			cout << dual_vars[sum+j-i] << endl;
		}
	sum += n-i;}	
	double sum_d = 0.0; for(int i=1; i<=n_rows_con_eq;i++) sum_d+=dual_vars[i];
	cout << "-> dual objective function value: " << sum_d << " (dual gap = " << StatsArray[0].LPRootLB*pow(2.0,-n)-sum_d << ")" << endl << endl;
	
	
	/*int max_mu_idx=1, min_mu_idx=1;
	for(int i=2; i<=n; i++){
		if(dual_vars[i]>dual_vars[max_mu_idx]) max_mu_idx = i;
		if(dual_vars[i]<dual_vars[min_mu_idx]) min_mu_idx = i;
	}
	double rebalance = (dual_vars[max_mu_idx]+dual_vars[min_mu_idx])/2;
	for(int i=1; i<=n; i++) sum_d -= dual_vars[i];
	dual_vars[max_mu_idx] = rebalance;
	for(int i=1; i<=n; i++)
		if((i!=max_mu_idx)&&(i!=min_mu_idx))
			dual_vars[i] -= rebalance - dual_vars[min_mu_idx];
	dual_vars[min_mu_idx] = rebalance;
	for(int i=1; i<=n; i++) sum_d += dual_vars[i];*/
	


	double alpha_rate=0.000001; double alpha_max=0.5;
	int counter=0;
	/*int i_hat=1;
	while(counter++<100){
	// dual step: increase nu and decrease each mu_(i_hat) by the largest possible margin
	double opt_alpha_inc, alpha_inc = alpha_rate, sum_n_d_max=0.0;
	while(alpha_inc<=alpha_max){
		double max_LB_beta=0.0;
		sum=n+1; for(int i=1; i<=n; i++){
			for(int j=i+1; j<=n; j++) for(int ell=2; ell<=n-1; ell++){
				double LB_beta = pow(2.0,ell-1)*dual_vars[sum+j-i]+dual_vars[i]+dual_vars[j]+(dual_vars[n+1]+alpha_inc)*ell/(2*n-3)-dist[i][j];
				if(LB_beta>max_LB_beta) max_LB_beta = LB_beta;}
			sum += n-i;
		}
		// check if the dual objective function value got improved
		double sum_n_d = sum_d + alpha_inc - max_LB_beta*n/2;
		if(sum_n_d > sum_n_d_max) {sum_n_d_max = sum_n_d; opt_alpha_inc = alpha_inc;}
		// next alpha
		alpha_inc+=alpha_rate;
	}
	//cout << sum_n_d_max << " | " << sum_d << endl;
	if((sum_n_d_max>sum_d)&&(opt_alpha_inc>alpha_rate)){
		// increase nu by opt_alpha_inc
		dual_vars[n+1] += opt_alpha_inc; 
		double max_LB_beta=0.0;
		sum=n+1; for(int i=1; i<=n; i++){
			for(int j=i+1; j<=n; j++) for(int ell=2; ell<=n-1; ell++){
				double LB_beta = pow(2.0,ell-1)*dual_vars[sum+j-i]+dual_vars[i]+dual_vars[j]+(dual_vars[n+1]+alpha_inc)*ell/(2*n-3)-dist[i][j];
				if(LB_beta>max_LB_beta) max_LB_beta = LB_beta;}
			sum += n-i;
		}
		// decrease mu_(i_hat) by the tightest feasibility bound on beta_(i_hat)
		for(int i=1; i<=n; i++) dual_vars[i] -= max_LB_beta/2;
		sum_d = sum_n_d_max; 
		cout << "-> dual objective function value (nu += " << opt_alpha_inc << "; mu decreased): " << sum_d << " (dual gap = " << StatsArray[0].LPRootLB*pow(2.0,-n)-sum_d << ")" << endl << endl;
		// if there was an increase of nu, the descrease of mu_(i_hat) lead to more space for nu to increase independently
		double max_nu=INF;
		sum=n+1; for(int i=1; i<=n-1; i++){for(int j=i+1; j<=n; j++)
		for(int ell=2; ell<=n-1; ell++){
			double UB_nu = (dist[i][j]-pow(2.0,ell-1)*dual_vars[sum+j-i]-dual_vars[i]-dual_vars[j])*(2*n-3)/ell;
			if(UB_nu < max_nu) max_nu = UB_nu;
		}
		sum += n-i;}
		if(sum_d + max_nu - dual_vars[n+1] > sum_d){
			sum_d += max_nu - dual_vars[n+1];
			cout << "-> dual objective function value (nu += " << (max_nu - dual_vars[n+1]) <<  "): " << sum_d << " (dual gap = " << StatsArray[0].LPRootLB*pow(2.0,-n)-sum_d << ")" << endl << endl;
			dual_vars[n+1] = max_nu;
		}
	} else break;

	}*/


	/*counter=0;
	while(counter++<100){
	// dual step: increase nu and decrease all lambdas by the largest possible margin
	double alpha_inc = alpha_rate;
	double opt_alpha_inc, sum_n_d_max=0.0;
	while(alpha_inc<=alpha_max){
		double beta_sum=0.0;
		sum=n+1; for(int i=1; i<=n-1; i++){ 
			for(int j=i+1; j<=n; j++){
				double max_LB_beta=0.0;
				for(int ell=2; ell<=n-1; ell++){
					double LB_beta = dual_vars[sum+j-i]+pow(2.0,1-ell)*(dual_vars[i]+dual_vars[j])+(dual_vars[n+1]+alpha_inc)*ell*pow(2.0,1-ell)/(2*n-3)-dist[i][j]*pow(2.0,1-ell);
					if(LB_beta>max_LB_beta) max_LB_beta = LB_beta;}
				beta_sum+=max_LB_beta;}
			sum += n-i;}
		// check if the dual objective function value got improved
		double sum_n_d = sum_d + alpha_inc - beta_sum;
		if(sum_n_d > sum_n_d_max) {sum_n_d_max = sum_n_d; opt_alpha_inc = alpha_inc;}
		// next alpha
		alpha_inc+=alpha_rate;
	}
	//cout << sum_n_d_max << " | " << sum_d << endl;
	if((sum_n_d_max>sum_d)&&(opt_alpha_inc>alpha_rate)){
		// increase nu by opt_alpha_inc
		dual_vars[n+1] += opt_alpha_inc;
		sum=n+1; for(int i=1; i<=n-1; i++){ 
			for(int j=i+1; j<=n; j++){
				double max_LB_beta=0.0;
				for(int ell=2; ell<=n-1; ell++){
					double LB_beta = dual_vars[sum+j-i]+pow(2.0,1-ell)*(dual_vars[i]+dual_vars[j])+dual_vars[n+1]*ell*pow(2.0,1-ell)/(2*n-3)-dist[i][j]*pow(2.0,1-ell);
					if(LB_beta>max_LB_beta) max_LB_beta = LB_beta;}
				// decrease lambda_ij by the tightest feasibility bound on beta_ij
				dual_vars[sum+j-i] -= max_LB_beta;}
			sum += n-i;}
		sum_d = sum_n_d_max;
		cout << "-> dual objective function value (nu += " << opt_alpha_inc << "; lambda decreased): " << sum_d << " (dual gap = " << StatsArray[0].LPRootLB*pow(2.0,-n)-sum_d << ")" << endl << endl;
		// if there was an increase of nu, the descrease of lambda_ij lead to more space for nu to increase independently
		double max_nu=INF;
		sum=n+1; for(int i=1; i<=n-1; i++){for(int j=i+1; j<=n; j++)
			for(int ell=2; ell<=n-1; ell++){
				double UB_nu = (dist[i][j]-pow(2.0,ell-1)*dual_vars[sum+j-i]-dual_vars[i]-dual_vars[j])*(2*n-3)/ell;
				if(UB_nu < max_nu) max_nu = UB_nu;
			}
		sum += n-i;}
		//max_nu = floor(max_nu*pow(10,12))*pow(10,-12); // numerical correction
		if(sum_d + max_nu - dual_vars[n+1] > sum_d){
			sum_d += max_nu - dual_vars[n+1];
			cout << "-> dual objective function value (nu += " << (max_nu - dual_vars[n+1]) <<  "): " << sum_d << " (dual gap = " << StatsArray[0].LPRootLB*pow(2.0,-n)-sum_d << ")" << endl << endl;
			dual_vars[n+1] = max_nu;
		}
	} else break;

	}*/

	int c=0; while(c++<10){
	sum=n+1; for(int i=1; i<=n-1; i++){
		for(int j=1; j<=n; j++) 
			if(dual_vars[sum+i-j]>=0) dual_vars[sum+i-j] = 0;
		sum += n-i;}

	counter=0;
	for(int i_hat=1; i_hat<=n; i_hat++) while(counter++<100){
		//if(((i_hat==1)&&(c==2))||((i_hat==n)&&(c==2))) break;
	// dual step: increase mu_(i_hat) and decrease all lambdas by the largest possible margin
	double opt_alpha_inc, alpha_inc = alpha_rate, sum_n_d_max=0.0;
	while(alpha_inc<=alpha_max){
		double beta_sum=0.0;
		sum=n+1; for(int j=1; j<=i_hat-1; j++){
			double max_LB_beta=0.0;
			for(int ell=2; ell<=n-1; ell++){
				double LB_beta = dual_vars[sum+i_hat-j]+pow(2.0,1-ell)*(dual_vars[i_hat]+alpha_inc+dual_vars[j])+dual_vars[n+1]*ell*pow(2.0,1-ell)/(2*n-3)-dist[j][i_hat]*pow(2.0,1-ell);
				if(LB_beta>max_LB_beta) max_LB_beta = LB_beta;}
			beta_sum+=max_LB_beta;
			sum += n-j;
		}
		for(int j=i_hat+1; j<=n; j++){
			double max_LB_beta=0.0;
			for(int ell=2; ell<=n-1; ell++){
				double LB_beta = dual_vars[sum+j-i_hat]+pow(2.0,1-ell)*(dual_vars[i_hat]+alpha_inc+dual_vars[j])+dual_vars[n+1]*ell*pow(2.0,1-ell)/(2*n-3)-dist[i_hat][j]*pow(2.0,1-ell);
				if(LB_beta>max_LB_beta) max_LB_beta = LB_beta;}
			beta_sum+=max_LB_beta;
		}
		// check if the dual objective function value got improved
		double sum_n_d = sum_d + alpha_inc - beta_sum;
		if(sum_n_d > sum_n_d_max) {sum_n_d_max = sum_n_d; opt_alpha_inc = alpha_inc;}
		// next alpha
		alpha_inc+=alpha_rate;
	}
	if((sum_n_d_max>sum_d)&&(opt_alpha_inc>alpha_rate)){
		// increase mu_(i_hat) by opt_alpha_inc
		dual_vars[i_hat] += opt_alpha_inc; 
		sum=n+1; for(int j=1; j<=i_hat-1; j++){
			double max_LB_beta=0.0;
			for(int ell=2; ell<=n-1; ell++){
				double LB_beta = dual_vars[sum+i_hat-j]+pow(2.0,1-ell)*(dual_vars[i_hat]+dual_vars[j])+dual_vars[n+1]*ell*pow(2.0,1-ell)/(2*n-3)-dist[j][i_hat]*pow(2.0,1-ell);
				if(LB_beta>max_LB_beta) max_LB_beta = LB_beta;}
			// decrease lambda_ij by the tightest feasibility bound on beta_ij
			dual_vars[sum+i_hat-j] -= max_LB_beta;
			sum += n-j;
		}
		for(int j=i_hat+1; j<=n; j++){
			double max_LB_beta=0.0;
			for(int ell=2; ell<=n-1; ell++){
				double LB_beta = dual_vars[sum+j-i_hat]+pow(2.0,1-ell)*(dual_vars[i_hat]+dual_vars[j])+dual_vars[n+1]*ell*pow(2.0,1-ell)/(2*n-3)-dist[i_hat][j]*pow(2.0,1-ell);
				if(LB_beta>max_LB_beta) max_LB_beta = LB_beta;}
			// decrease lambda_ij by the tightest feasibility bound on beta_ij
			dual_vars[sum+j-i_hat] -= max_LB_beta;
		}
		sum_d = sum_n_d_max; 
		cout << "-> dual objective function value (mu_"<<i_hat<<" += " << opt_alpha_inc << "; lambda decreased): " << sum_d << " (dual gap = " << StatsArray[0].LPRootLB*pow(2.0,-n)-sum_d << ")" << endl << endl;
			// if there was an increase of mu_i, the descrease of lambda_ij lead to more space for mu_i to increase independently
			double max_mu=INF;
			sum=n+1; for(int j=1; j<=i_hat-1; j++){
				for(int ell=2; ell<=n-1; ell++){
					double mu = dist[j][i_hat]-pow(2.0,ell-1)*dual_vars[sum+i_hat-j]-dual_vars[j]-ell*dual_vars[n+1]/(2*n-3);
					if(mu < max_mu) max_mu = mu;}
				sum += n-j;
			}
			for(int j=i_hat+1; j<=n; j++)
				for(int ell=2; ell<=n-1; ell++){
					double mu = dist[i_hat][j]-pow(2.0,ell-1)*dual_vars[sum+j-i_hat]-dual_vars[j]-ell*dual_vars[n+1]/(2*n-3);
					if(mu < max_mu) max_mu = mu;}
			
			if(sum_d + max_mu - dual_vars[i_hat] > sum_d){
				sum_d += max_mu - dual_vars[i_hat];
				cout << "-> dual objective function value (mu_"<<i_hat<<" += " << (max_mu - dual_vars[i_hat]) <<  "): " << sum_d << " (dual gap = " << StatsArray[0].LPRootLB*pow(2.0,-n)-sum_d << ")" << endl << endl;
				dual_vars[i_hat] = max_mu;
			}
	} else break;

	}

	
	// nu
	cout << "nu:" << endl << dual_vars[n+1] << endl;
	// mu
	cout << "mu:" << endl;
	for(int i_hat=1; i_hat<=n; i_hat++) cout << dual_vars[i_hat] << endl;
	// lambda
	cout << "lambda:" << endl;
	sum=n+1; for(int i=1; i<=n-1; i++){ for(int j=i+1; j<=n; j++)
			cout << dual_vars[sum+j-i] << endl;
	sum += n-i;}
		

	double sum_n_d=0.0; for(int i=1; i<=n_rows_con_eq; i++) sum_n_d+=dual_vars[i];
	cout << "-> dual objective function value " << sum_n_d << " (dual gap = " << StatsArray[0].LPRootLB*pow(2.0,-n)-sum_n_d << ")" << endl << endl;
	
	}

	sum=n+1;
	for(int i=1; i<=n-1; i++){
		for(int j=i+1; j<=n; j++)
			for(int ell=2; ell<=n-1; ell++)
				cout << "2^" << ell << " * " << dual_vars[sum+j-i] << "+2(" << dual_vars[i] 
					 << "+" << dual_vars[j] << ")+2*" << ell << "/(2n-3) * " << dual_vars[n+1]
					 << " | " << 2*dist[i][j]-pow(2.0,ell)*dual_vars[sum+j-i]-2*(dual_vars[i]+dual_vars[j])-2*ell*dual_vars[n+1]/(2*n-3) << endl;
		sum += n-i;
	}

	sum=n+1;
	for(int i=1; i<=n-1; i++){
		for(int j=i+1; j<=n; j++)
			for(int ell=2; ell<=n-1; ell++)
				if(2*dist[i][j]-pow(2.0,ell)*dual_vars[sum+j-i]-2*(dual_vars[i]+dual_vars[j])-2*ell*dual_vars[n+1]/(2*n-3)==0){
					cout << "mu_" << i << " = " << dual_vars[i] << " vs. mu_" << j << " = " << dual_vars[j] << endl; 
				}
		sum += n-i;
	}


	/*int counter=0;
	while(counter++<10){
		cout << setprecision(10) << "dual gaps in step " << counter << ":" << endl;*/
		/*for(int i=1; i<=n-1; i++)
			for(int j=i+1; j<=n; j++){
				for(int ell=2; ell<=n-1; ell++) cout << setprecision(10) << dual_gaps[i][j][ell] << " ";
				cout << endl;
			}*/
		// DUAL STEP
		// find minimum violated Kraft equality
		/*double min=INF; int i_hat=-1;
		for(int i=1; i<=n; i++){
			double sum=0.0;
			for(int j=1; j<=n; j++)
				for(int ell=2; ell<=n-1; ell++){
					if(i<j) sum += pow(2.0,1-ell)*primal_vars[i][j][ell];
					if(i>j) sum += pow(2.0,1-ell)*primal_vars[j][i][ell];
				}
			if(abs(sum-1)<min) {min = abs(sum-1); i_hat=i;}
		}
		cout << endl << "minimum violated Kraft equality: " << i_hat << endl;*/
		// find minimum epsilon
		/*double min_eps=INF; sum=n+1;
		for(int i=1; i<=n-1; i++){
			for(int j=i+1; j<=n; j++)
				for(int ell=2; ell<=n-1; ell++)
					if(primal_vars[i][j][ell]==0){
						double eps = dist[i][j] - dual_vars[sum+j-i]*pow(2.0,ell) - 2*(dual_vars[i]+dual_vars[j]) - dual_vars[n+1]*ell/(2*n-3);
						if(eps<min_eps) min_eps = eps;
					}
			sum += n-i;
		}
		cout << endl << "minimum epsilon for inactive inequalities: " << min_eps << endl;
		for(int i=1; i<=n-1; i++)
			for(int j=i+1; j<=n; j++){
				for(int ell=2; ell<=n-1; ell++) cout << dual_gaps[i][j][ell] << " ";
				//if((i==i_hat)||(j==i_hat)) cout << "*" << endl; else 
				cout << endl;
			}*/
		// find the minimum dual gap
		/*min=INF; int min_j=-1;
		for(int j=1; j<=n; j++)
			for(int ell=2; ell<=n-1; ell++){
				if((i_hat<j)&&(dual_gaps[i_hat][j][ell]>0)&&(dual_gaps[i_hat][j][ell]<min)) {min=dual_gaps[i_hat][j][ell]; min_j=j;}
				if((i_hat>j)&&(dual_gaps[j][i_hat][ell]>0)&&(dual_gaps[j][i_hat][ell]<min)) {min=dual_gaps[j][i_hat][ell]; min_j=j;}
			}
		cout << endl << "minimum dual gap: " << min << endl;*/
		// increase mu_{i_hat} until at least one new dual constraints is active
		/*dual_vars[i_hat] += min/2;
		for(int j=1; j<=n; j++)
			for(int ell=2; ell<=n-1; ell++){
				if(i_hat<j) dual_gaps[i_hat][j][ell] -= min;
				if(i_hat>j) dual_gaps[j][i_hat][ell] -= min;
				//if(i_hat<j) {cout << setprecision(4) << ell << " | " << dual_gaps[i_hat][j][ell]; dual_gaps[i_hat][j][ell] -= min; cout << " | " << dual_gaps[i_hat][j][ell] << endl;}
				//if(i_hat>j) {cout << setprecision(4) << ell << " | " << dual_gaps[j][i_hat][ell]; dual_gaps[j][i_hat][ell] -= min; cout << " | " << dual_gaps[j][i_hat][ell] << endl;}		
			} 
		for(int i=1; i<=n-1; i++)
			for(int j=i+1; j<=n; j++){
				for(int ell=2; ell<=n-1; ell++) cout << dual_gaps[i][j][ell] << " ";
				if((i==i_hat)||(j==i_hat)) cout << "*" << endl; else cout << endl;
			}
		// increase lambda_{i_hat,j} until all dual gaps are non-negative
		sum=0; for(int i=1; i<=i_hat-1; i++) sum += n-i;
		for(int j=1; j<=n; j++) if(j!=min_j){
			min=0.0; for(int ell=2; ell<=n-1; ell++){
				if((i_hat<j)&&(dual_gaps[i_hat][j][ell]*pow(2.0,-ell)<min)) min=dual_gaps[i_hat][j][ell]*pow(2.0,-ell);
				if((i_hat>j)&&(dual_gaps[j][i_hat][ell]*pow(2.0,-ell)<min)) min=dual_gaps[j][i_hat][ell]*pow(2.0,-ell);
			}
			for(int ell=2; ell<=n-1; ell++){
				if(i_hat<j) dual_gaps[i_hat][j][ell] -= min*pow(2.0,ell);
				if(i_hat>j) dual_gaps[j][i_hat][ell] -= min*pow(2.0,ell);
			}
			if(i_hat<j) dual_vars[sum+j-i_hat] += min;
			if(i_hat>j) dual_vars[sum+j-i_hat] += min;
		}*/
		/*for(int j=1; j<=n; j++)
			for(int ell=2; ell<=n-1; ell++){
				if(i_hat<j) {cout << setprecision(4) << ell << " | " << dual_gaps[i_hat][j][ell] << endl; }
				if(i_hat>j) {cout << setprecision(4) << ell << " | " << dual_gaps[j][i_hat][ell] << endl;}		
			}*/ 
		/*cout << endl << "updated lambda" << endl;
		for(int i=1; i<=n-1; i++)
			for(int j=i+1; j<=n; j++){
				for(int ell=2; ell<=n-1; ell++) cout << dual_gaps[i][j][ell] << " ";
				if((i==i_hat)||(j==i_hat)) cout << "*" << endl; else cout << endl;
			}*/
		// increase all dual variables by min_eps
		/*for(int i=1; i<=n_rows_con_eq; i++) dual_vars[i] += min_eps;
		for(int i=1; i<=n_rows_con_eq; i++)
			cout << dual_vars[i] << endl;
		// update all dual gaps
		sum=n+1;
		for(int i=1; i<=n-1; i++){
			for(int j=i+1; j<=n; j++)
				for(int ell=2; ell<=n-1; ell++)
					dual_gaps[i][j][ell] = dist[i][j] - dual_vars[sum+j-i]*pow(2.0,ell) - 2*(dual_vars[i]+dual_vars[j]) - dual_vars[n+1]*ell/(2*n-3);
			sum += n-i;
		}*/
		// update lambda
		/*sum=n+1;
		for(int i=1; i<=n-1; i++){
			for(int j=i+1; j<=n; j++) {
				double min=INF; int opt_ell=-1;
				for(int ell=2; ell<=n-1; ell++)
					if(dual_gaps[i][j][ell]*pow(2.0,-ell)<min){min = dual_gaps[i][j][ell]*pow(2.0,-ell); opt_ell=ell;}
				dual_vars[sum+j-i] += min;//dist[i][j]/4 - min_d[i]/8 - min_d[j]/8;
				for(int ell=2; ell<=n-1; ell++) dual_gaps[i][j][ell] -= min*pow(2.0,ell);
			}
			sum += n-i;
		}*/
		/*cout << endl << "dual variable += epsilon: " << min_eps << endl;
		for(int i=1; i<=n-1; i++)
			for(int j=i+1; j<=n; j++){
				for(int ell=2; ell<=n-1; ell++) cout << dual_gaps[i][j][ell] << " ";
				//if((i==i_hat)||(j==i_hat)) cout << "*" << endl; else 
				cout << endl;
			}
		// PRIMAL STEP
		cout << "primal solution: " << endl; 
		for(int i=1; i<=n-1; i++)
			for(int j=i+1; j<=n; j++){
				double min=INF;
				for(int ell=2; ell<=n-1; ell++){
					if(dual_gaps[i][j][ell]<min) min = dual_gaps[i][j][ell];*/
					//primal_vars[i][j][ell] = 1.0/(n-2);
					//cout << primal_vars[i][j][ell] << " | " << dual_gaps[i][j][ell] << endl;
					/*if(dual_gaps[i][j][ell]==0){//&&(2*dist[i][j]!=min_d[i]+min_d[j]))
						primal_vars[i][j][ell] = 1; 
					//else if(dual_gaps[i][j][ell]==0) 
					//	primal_vars[i][j][ell] = 1.0/(n-2);
					} else
						primal_vars[i][j][ell] = 0;*/
					//cout << setprecision(5) << primal_vars[i][j][ell] << " ";
				//}
				/*for(int ell=2; ell<=n-1; ell++){
					if(dual_gaps[i][j][ell]==min) primal_vars[i][j][ell] = 1; 
					else primal_vars[i][j][ell] = 0; 
					cout << setprecision(5) << primal_vars[i][j][ell] << " ";
				}
				cout << endl;
			}
		cout << endl;
	}*/
}

double Solver::EntropyBound(int n, int taxon, double * &d){
	if(n==taxon){
		double sum=0.0;
		for(int j=2;j<=taxon;j++)
			sum += (d[j]-opt_alpha*alpha_sigma[j])/pow(2.0,alpha_sigma[j]);
		return sum;
	} else {
		// build vector d^alpha
		double* d_alpha = new double[n];
		for(int i=1; i<=n-2; i++) d_alpha[i]=d[i];
		double value=(d[n-1]+d[n])/2-opt_alpha;
		// sort d^alpha
		int idx = 1; for(int i=1; i<=n-2; i++) if(d_alpha[i]>value){idx=i;break;}
		for(int i=n-1; i>idx; i--) d_alpha[i] = d_alpha[i-1]; d_alpha[idx]=value;
		// continue recursion on d^alpha
		double bound = EntropyBound(n-1,taxon,d_alpha);
		delete [] d_alpha;
		// reverse the sorting according to d^alpha
		double sigma_value=alpha_sigma[idx];
		for(int i=idx; i<=n-2; i++) alpha_sigma[i]=alpha_sigma[i+1]; alpha_sigma[n-1]=sigma_value;
		// expand the path-length sequence
		alpha_sigma[n-1]++; alpha_sigma[n]=alpha_sigma[n-1];
		return bound;
	}
}
	
// Assumption: dist is not doubly-stochastic, matrix Tau is well-defined
void Solver::EntropyAnalysis(int taxon){
	//cout << "start entropy analysis for taxon " << taxon << " with Tau matrix " << endl;
	//for(int i=1; i<=taxon; i++){ for(int j=1; j<=taxon; j++) cout << Tau[i][j] << " "; cout << endl;}
	double *d_alpha=new double[n+1]; 
	int *pi=new int[n+1];
	int *pi2=new int[n+1];
	alpha_sigma=new int[n+1];
	int **sigma=new int*[n+1]; for(register unsigned int i=1; i<=n; i++) sigma[i] = new int[n+1];
	double entropybound=0.0;

	for(register unsigned int k=1; k<=taxon; k++){
		// pick k-th row of distance matrix
		for(register unsigned int j=1; j<=n; j++){d_alpha[j]=dist[k][j]; pi[j]=j;}
		// sort k-th row of distance matrix
		for(int i=1; i<=n-1; i++){int min = i;
			for(int j=i+1; j<=n; j++) if(d_alpha[min]>d_alpha[j]) min=j;
			double value = d_alpha[min]; int idx = pi[min];
			while(min > i){d_alpha[min]=d_alpha[min-1]; pi[min]=pi[min-1]; min--;}
			d_alpha[i]=value; pi[i]=idx;}
		// define the optimal solution for f_taxon
		alpha_sigma[1]=0;
		for(int i=1; i<k; i++)
			alpha_sigma[i+1]=Tau[k][i]-1; 
		for(int i=k+1; i<=taxon; i++)
			alpha_sigma[i]=Tau[k][i]-1; 
		// sort alpha_sigma
		for(int i=1; i<=taxon; i++){int min = i;
			for(int j=i+1; j<=taxon; j++) if(alpha_sigma[min]>alpha_sigma[j]) min=j;
			double value = alpha_sigma[min]; int idx = pi2[min];
			while(min > i){alpha_sigma[min]=alpha_sigma[min-1]; pi2[min]=pi2[min-1]; min--;}
			alpha_sigma[i]=value; pi2[i]=idx;}
		//cout << "sorted alpha_sigma for taxon " << taxon << endl;
		//for(int i=1; i<=taxon; i++) cout << alpha_sigma[i] << " "; cout << endl;
		// calculate the entropy lower bound
		entropybound += EntropyBound(n,taxon,d_alpha);
		//cout << "sorted alpha_sigma after entropy calculations " << endl;
		//for(int i=1; i<=n; i++) cout << alpha_sigma[i] << " "; cout << endl;
		//cout << "the entropy lower bound for row " << k << " is: " << entropy_bound << endl;
		// reverse sort
		for(register unsigned int j=1; j<=n; j++) sigma[k][pi[j]] = alpha_sigma[j];
		//cout << "sigma after entropy calculations " << endl;
		//for(int i=1; i<=n; i++) cout << sigma[k][i] << " "; cout << endl;
		/*for(register unsigned int j=1; j<=n; j++) cout << pi[j] << " "; cout << endl;
		for(register unsigned int j=1; j<=n; j++) cout << d_alpha[j] << " "; cout << endl;
		for(register unsigned int j=1; j<=n; j++) cout << alpha_tau[j] << " "; cout << endl;*/
		//for(register unsigned int j=1; j<=n; j++) cout << sigma[k][j] << " "; cout << endl;
	}
	//cout << "the entropy lower bound is: " << entropybound << endl;
	entropybound+=opt_alpha*(3*n-6);
	entropybound/=2;
	if(entropybound>opt_entropybound) opt_entropybound=entropybound;
	opt_entropybound*=precomputedpow[n];
	/*for(register unsigned int i=1; i<=n;i++){
		for(register unsigned int j=1; j<=n;j++)
			cout << sigma[i][j] << " "; 
		cout << endl;
	}*/
	delete [] pi; delete [] pi2; delete [] d_alpha; for(register unsigned int i=1; i<=n;i++) delete [] sigma[i]; delete [] sigma;
}

// solves, for each row i, min{ sum_i dist_i * 2^-tau_i : tau_i path-length sequence}
double Solver::KraftLowerBound(){
	double bound=0.0;
	double *d=new double[n+1];
	for(register unsigned int i=1; i<=n; i++){
		for(register unsigned int j=1; j<=n; j++) d[j]=dist[i][j];
		for(int k=1; k<=n-1; k++){
			int min = k;
			for(int j=k+1; j<=n; j++) if(d[min]>d[j]) min=j;
			double value = d[min];
			while(min > k){d[min]=d[min-1]; min--;}
			d[k]=value;
		}
		for(register unsigned int j=2; j<=n-1; j++)	bound+=d[j]/pow(2.0,j);
		bound+=d[n]/pow(2.0,n-1);
	}
	return bound*pow(2.0,n);
}


double Solver::PardiLowerBound(int taxon){
	double acc=0; for(register unsigned int k=taxon+1; k<=n; k++) acc+=beta[k];
	acc=precomputedpow[n-1]*acc + len;
	return acc;
}

double Solver::PardiGammaLowerBound(int taxon){
	double acc=len; for(register unsigned int k=taxon+1; k<=n; k++) acc+=bgamma[k];
	return acc;
}

double Solver::LagrangianDualLowerBound(int taxon){
	double LB1 = DualValuesAtLeaf[taxon][0]*(2*n-3)*precomputedpow[n];
	for(register unsigned int i=1; i<=n; i++) LB1 += DualValuesAtLeaf[taxon][i]*precomputedpow[n-1];
	
	double temp,temp1;
	int temp2;
	for(register unsigned int i=1; i<=taxon-1; i++) 
		for(register unsigned int j=i+1; j<=taxon; j++) {   
			temp = 2*dist[i][j]-DualValuesAtLeaf[taxon][i]-DualValuesAtLeaf[taxon][j];
			temp1 = (temp - 2*DualValuesAtLeaf[taxon][0]*Tau[i][j])*precomputedpow[n-Tau[i][j]];
			temp2 = MIN(n-taxon+Tau[i][j],n-1);			
			for(register unsigned int k = Tau[i][j]+1; k <= temp2; k++) 
				if ((temp - 2*DualValuesAtLeaf[taxon][0]*k)*precomputedpow[n-k] < temp1) temp1 = (temp - 2*DualValuesAtLeaf[taxon][0]*k)*precomputedpow[n-k];				
			LB1 += temp1;
		}		

	for(register unsigned int i=1; i<=taxon-1; i++){ 
		for(register unsigned int j=taxon+1; j<=n; j++){
			temp = 2*dist[i][j]-DualValuesAtLeaf[taxon][i]-DualValuesAtLeaf[taxon][j];
			temp1 = (temp - 2*DualValuesAtLeaf[taxon][0]*2)*precomputedpow[n-2];
			temp2 = MIN(n-taxon+Tau[i][0],n-1);	// MODIFIED 2/6/2009		
			for(register unsigned int k = 3; k <= temp2; k++)  // MODIFIED 2/6/2009
				if ((temp - 2*DualValuesAtLeaf[taxon][0]*k)*precomputedpow[n-k] < temp1) temp1 = (temp - 2*DualValuesAtLeaf[taxon][0]*k)*precomputedpow[n-k];
			LB1 += temp1;
		}
	}
	for(register unsigned int i=taxon; i<=n-1; i++){ 
		for(register unsigned int j=i+1; j<=n; j++){
			temp = 2*dist[i][j]-DualValuesAtLeaf[taxon][i]-DualValuesAtLeaf[taxon][j];
			temp1 = (temp - 2*DualValuesAtLeaf[taxon][0]*2)*precomputedpow[n-2];
			temp2 = MIN(n-taxon+Tau[0][0],n-1); // MODIFIED 2/6/2009	
			for(register unsigned int k = 3; k <= n-1; k++) // MODIFIED 2/6/2009
				if ((temp - 2*DualValuesAtLeaf[taxon][0]*k)*precomputedpow[n-k] < temp1) temp1 = (temp - 2*DualValuesAtLeaf[taxon][0]*k)*precomputedpow[n-k];
			LB1 += temp1;
		}
	}
	return LB1;
}


double Solver::FastLinearProgrammingLowerBound(int taxon){
	// Eq 13
	for(register unsigned int i=1; i<=taxon-1; i++){
		for(register unsigned int j=i+1; j<=taxon; j++){
			for(register unsigned int k=2; k<=Tau[i][j]-1; k++) path[i][j][k].setUB(0.0); 
			for(register unsigned int k=(n-taxon)+Tau[i][j]+1; k<=n-1; k++) path[i][j][k].setUB(0.0); 
		}
	}
	
	// Eq 14
	for(register unsigned int i=1; i<=taxon; i++){
		for(register unsigned int j=taxon+1; j<=n; j++){
			for(register unsigned int k=(n-taxon)+Tau[taxon][0]+1; k<=n-1; k++) path[i][j][k].setUB(0.0); 
		}
	}

	// Eq 15
	for(register unsigned int i=taxon+1; i<=n-1; i++){
		for(register unsigned int j=i+1; j<=n;   j++){
			for(register unsigned int k=(n-taxon)+Tau[0][0]+1; k<=n-1; k++) path[i][j][k].setUB(0.0); 
		}
	}

	problem.minim("ld");
	double LB=problem.getObjVal(); 
	DualValuesAtLeaf[taxon][0]=Third.getDual(); 
	for(register unsigned int i=1; i <= n; i++) DualValuesAtLeaf[taxon][i]=Kraft[i].getDual();

	// Eq 13
	for(register unsigned int i=1; i<=taxon-1; i++){
		for(register unsigned int j=i+1; j<=taxon; j++){
			for(register unsigned int k=2; k<=Tau[i][j]-1; k++) path[i][j][k].setUB(1.0); 
			for(register unsigned int k=(n-taxon)+Tau[i][j]+1; k<=n-1; k++) path[i][j][k].setUB(1.0); 
		}
	}
	
	// Eq 14
	for(register unsigned int i=1; i<=taxon; i++){
		for(register unsigned int j=taxon+1; j<=n; j++){
			for(register unsigned int k=(n-taxon)+Tau[taxon][0]+1; k<=n-1; k++) path[i][j][k].setUB(1.0); 
		}
	}

	// Eq 15
	for(register unsigned int i=taxon+1; i<=n-1; i++){
		for(register unsigned int j=i+1; j<=n;   j++){
			for(register unsigned int k=(n-taxon)+Tau[0][0]+1; k<=n-1; k++) path[i][j][k].setUB(1.0); 
		}
	}

	return LB; 
}

/*double Solver::LinLog(int taxon){
	double sum=0.0;
	for(register unsigned int i=1; i<=taxon; i++){
		for(register unsigned int j=i+1; j<=taxon; j++){
			sum+=precomputedpow[n-Tau[i][j]]*(1-pow(2.0,-dist[i][j]));
		}
	}
	return sum;
}*/

bool Solver::Pruned(int taxon, int selector, double bound){
	double LowerBound=0; 
	switch(selector){ 
		case 1: LowerBound=PardiGammaLowerBound(taxon);
				if(LowerBound > Optimum.tree_len) return true; 
				break;
		case 2:	EntropyAnalysis(taxon);
				//if(opt_entropybound > Optimum.tree_len) return true;
				if(!obj_rescaling){
					StatsArray[Core].DualCounter++; 
					if(eorder!=1) LowerBound=LagrangianDualLowerBound(taxon);
					if(LowerBound > Optimum.tree_len) return true;
				}
				StatsArray[Core].LPCounter++; 
				if(eorder==1) LowerBound = bound;
				else LowerBound=FastLinearProgrammingLowerBound(taxon);
				if(LowerBound > Optimum.tree_len) return true;
				break;	
    }
	return false;  	 
}

void Solver::Explore(int start, int end, int taxon){
	/*#pragma omp critical
    {
    	cout << "Core "<<Core<<": " << taxon << " taxa inserted: Explore in (" << start << "," << end << ")" << endl;
    }*/
	double tmp_insertion_bound;
	for(register unsigned int e=start; e<end; e++){
		tmp_insertion_bound=tree[e].prio;
		/*#pragma omp critical
    	{
    		cout << "add taxon "<<taxon<<endl;
    	}*/
	 	Add(taxon,e);
	 	ComputeTopologicalMatrix(taxon);
		StatsArray[Core].nodecounter++;
		DFS(taxon,tmp_insertion_bound);
		Restore(e);		
		/*#pragma omp critical
    	{
    		cout << "remove taxon "<<taxon<<endl;
    	}*/
		tree[e].prio = tmp_insertion_bound;	
	}	
}
	   		

void Solver::insertPartition(int taxon, double sep, int n_sep){
	if(n_sep>0)
		if(ThereAreIdleCores(OverallCoreNumber) || Queue.size()<OverallCoreNumber){
			bool inserted_properly=false; 
			int shift=0; if(floor(sep)==ceil(sep) && (n_sep+1)*sep < (int)tree.size()) shift=1;
			if(floor(n_sep*sep)+1<ceil((n_sep+1)*sep)){
				while(!inserted_properly && !Out_of_Time())
					inserted_properly=Push(floor(n_sep*sep)+1,ceil((n_sep+1)*sep)+shift,taxon);
				/*#pragma omp critical
    			{
    				cout << "n_sep = " << n_sep << ", sep = "<<sep<<": "<< floor(n_sep*sep)+1 << " < " << ceil((n_sep+1)*sep)<<"; leads to "<<inserted_properly << endl;
    			}*/
			}
			if(inserted_properly && n_sep>1) insertPartition(taxon,sep,n_sep-1);
			else Explore(0,floor(n_sep*sep)+1,taxon+1);
		} else Explore(0,ceil((n_sep+1)*sep),taxon+1);
	else
		Explore(0,1,taxon+1);
}

void Solver::RandomizeTaxaOrder(){
	vector <int> new_vertices;
	vector <int> v;
	// create a random order of taxa
	for(int i=1;i<=n;i++)v.push_back(i);
	while(v.size()){
		srand(time(NULL));
		int idx=rand()%((int)v.size());
		int num = v[idx];
		swap(v[idx],v[(int)v.size()-1]);
		v.pop_back();
		new_vertices.push_back(num);
	}
	// reorder taxa according to the order 'new_vertices'
	double **new_dist=new double*[n+1]; for(register unsigned int i=0; i<=n;i++) new_dist[i] = new double[n+1];
	for(register unsigned int i=0; i<new_vertices.size()-1; i++){
		new_dist[i+1][i+1]=0;
		for(register unsigned int j=i+1; j<new_vertices.size(); j++){
			new_dist[i+1][j+1]=dist[new_vertices[i]][new_vertices[j]];
			new_dist[j+1][i+1]=new_dist[i+1][j+1];
		}
	}
	for(register unsigned int i=0; i<=n; i++) delete [] dist[i]; delete [] dist;
	dist=new_dist;
	new_vertices.clear();
	v.clear();
}

bool Solver::compare_edge_prio(EDGE e1, EDGE e2){return (e1.prio<e2.prio);}

int Solver::LBEdgeOrder(int taxon){
	int size=(int)tree.size();
	vector <struct EDGE> InsertedLowerBound; 
	double tmp_future_insertion_bound;
	for(register unsigned int e=0; e<size; e++){
		EDGE edge0;
		edge0.i=tree[e].i;
		edge0.j=tree[e].j;
		tmp_future_insertion_bound=tree[e].prio;
		Add(taxon+1,e); 
	 	ComputeTopologicalMatrix(taxon+1);
		// determine lower bounds after insertion of taxon+1 on edge e
		edge0.prio=LagrangianDualLowerBound(taxon+1);
	 	if(edge0.prio<=Optimum.tree_len) edge0.prio=FastLinearProgrammingLowerBound(taxon+1);
	 	if(edge0.prio>Optimum.tree_len) edge0.prio=MAX_DOUBLE; // insertion of taxon+1 on edge e prunes the search-tree
	 	InsertedLowerBound.push_back(edge0);
		Restore(e);	
		tree[e].prio = tmp_future_insertion_bound;	
	}
	// sort all lower bounds
	stable_sort(InsertedLowerBound.begin(),InsertedLowerBound.end(),compare_edge_prio);
	// reorder tree accordingly
	int idx=size;
	old_tree[taxon].clear();
	for(register unsigned int e=0; e<size; e++){
		EDGE edge0;
		edge0.i=tree[e].i;
		edge0.j=tree[e].j;
		edge0.prio=tree[e].prio;
		old_tree[taxon].push_back(edge0);
		tree[e].i=InsertedLowerBound[e].i;
		tree[e].j=InsertedLowerBound[e].j;
		tree[e].prio=InsertedLowerBound[e].prio;
		if((idx==size)&&(InsertedLowerBound[e].prio==MAX_DOUBLE))idx=e;
	}
 	InsertedLowerBound.clear();
 	return idx;
}

void Solver::RandomizeEdgeOrder(int taxon){
	int size=(int)tree.size();
	vector <struct EDGE> new_edges;
	vector <int> v;
	// create a random order of edges
	for(int i=0;i<size;i++)v.push_back(i);
	while(v.size()){
		srand(time(NULL));
		int idx=rand()%((int)v.size());
		int num = v[idx];
		swap(v[idx],v[(int)v.size()-1]);
		v.pop_back();
		// define the order of edges of tree according to the order v
		EDGE edge0;
		edge0.i=tree[num].i;
		edge0.j=tree[num].j;
		new_edges.push_back(edge0);
	}
	// save a copy of the current tree to restore it
	old_tree[taxon].clear();
	for(register unsigned int e=0; e<size; e++){
		EDGE edge0;
		edge0.i=tree[e].i;
		edge0.j=tree[e].j;
		old_tree[taxon].push_back(edge0);
	}
	// reorder tree
	for(register unsigned int e=0; e<size; e++){
		tree[e].i=new_edges[e].i;
		tree[e].j=new_edges[e].j;
	}
	new_edges.clear();
	v.clear();
}

void Solver::RestoreEdgeOrder(int taxon){
	int size=(int)tree.size();
	// restore the non-random order of the edges of tree
	for(int k=0;k<size;k++){
		tree[k].i=old_tree[taxon][k].i;
		tree[k].j=old_tree[taxon][k].j;
		tree[k].prio=old_tree[taxon][k].prio;
	}
}

void Solver::DFS(int taxon, double bound){ 
	if(Out_of_Time()){StatsArray[Core].outoftime=true; return;}
	if(taxon==n){
		StatsArray[Core].int_sol_counter++;
		#pragma omp critical
		{
			if(len < Optimum.tree_len){
				UpdateOptimum(); 
			}
		}
		NNI(n,true);              
	}
    else if (!Pruned(taxon,boundtype,bound)){ 
    	// the relative size of the tree; might exclude edges
    	/*cout << "current tree:" << endl;
    	for(int e=0; e<(int)tree.size();e++)
    		cout << "(" << tree[e].i << "," << tree[e].j << ")" << endl;*/
		int rel_size=(int)tree.size();
		if(eorder==1) rel_size=LBEdgeOrder(taxon);
		if(eorder==2) RandomizeEdgeOrder(taxon);
	   	if(rel_size>=1){
	   		double sep = 1.0;
	   		if(rel_size>=n_p+1){ // default: >= 4
	   			sep = rel_size/(n_p*1.0); // sep = uniform size of partition sets (>= 2 by default)
	   			/*#pragma omp critical
    			{
    				cout << "Calculated: rel_size = " << rel_size << ", sep = "<<sep << endl;
    			}*/
	   			if(sep>=s_p) insertPartition(taxon,sep,n_p-1);
	   			/*#pragma omp critical
    			{
    				cout << "finished inserting partitions" << endl;
    			}*/
	   		} else if(s_p==1) insertPartition(taxon,sep,rel_size-1);
	   		if((s_p>sep)&&(s_p<=floor(rel_size/2))){
	   			int n_pe = floor(rel_size/(s_p*1.0));
	   			sep = rel_size/(n_pe*1.0);
	   			insertPartition(taxon,sep,n_pe-1);
	   		}
	   	}	
		if((eorder==1)||(eorder==2)) RestoreEdgeOrder(taxon);	
	}		
	return;
}

void Solver::UpdateOptimum(){
	lock_guard<mutex> lg(m); //blocking access to the optimum data...
	Optimum.core=Core;
	Optimum.tree_len=len;
	Optimum.primal_bound_list.push_back(Optimum.tree);
	Optimum.OverallNumber_of_Primal_Solutions++;
	Optimum.tree.clear();
	struct EDGE edge;
	for(register unsigned int e=0; e<(int)tree.size(); e++){
		edge.i=tree[e].i;
		edge.j=tree[e].j;
		Optimum.tree.push_back(edge);
	} 
	cout<<"* "
			<<setw(numWidth2)<< setfill(separator)<<Optimum.OverallNumber_of_Primal_Solutions
			<<setw(numWidth+PRECISION) << setfill(separator)<<std::setprecision(PRECISION)<<Optimum.tree_len<<std::setprecision(PRECISION)
			<<setw(numWidth+PRECISION) << setfill(separator)<<std::scientific<<abs(Optimum.tree_len-StatsArray[0].RootLB)/Optimum.tree_len*100
			<<setprecision(8)
			<<setw(numWidth) << setfill(separator)<<std::fixed<<omp_get_wtime() - StatsArray[0].start_time
			<<setw(numWidth) << setfill(separator)<<Queue.size()
		    <<setw(numWidth) << setfill(separator)<<"no"
			<<endl;		
}

void Solver::UpdateOptimum2(bool print_line){
	lock_guard<mutex> lg(m); //blocking access to the optimum data...
	Optimum.tree_len=len;
	Optimum.primal_bound_list.push_back(Optimum.tree);
	Optimum.OverallNumber_of_Primal_Solutions++;
	Optimum.tree.clear();
	EDGE edge;
	register unsigned int tsize=tree.size();
	for(register unsigned int e=0; e<tsize; e++){
		edge.i=tree[e].i; 
		edge.j=tree[e].j;
		Optimum.tree.push_back(edge);
	}	
	if(print_line){
		cout<<"* "
			<<setw(numWidth2)<< setfill(separator)<<Optimum.OverallNumber_of_Primal_Solutions
			<<setw(numWidth+PRECISION) << setfill(separator)<<std::setprecision(PRECISION)<<Optimum.tree_len<<std::setprecision(PRECISION)
			<<setw(numWidth+PRECISION) << setfill(separator)<<std::scientific<<abs(Optimum.tree_len-StatsArray[0].RootLB)/Optimum.tree_len*100
			<<setprecision(8)
			<<setw(numWidth) << setfill(separator)<<std::fixed<<omp_get_wtime() - StatsArray[0].start_time
			<<setw(numWidth) << setfill(separator)<<Queue.size()
			<<setw(numWidth) << setfill(separator)<<"yes"
			<<endl;
	}		
}


/*NNI BLOCK*/
void Solver::NNISwap(int node1, int node2, int node3, int node4){
	// cout<<"swapping nodes "<<node1<<" and "<<node2<<" with "<<node3<<" and "<<node4<<endl;
	register unsigned int tsize=tree.size();
	for(register unsigned int e=0; e<tsize; e++){
		if(tree[e].i==node1 && tree[e].j==node2) {tree[e].j=node4; if(tree[e].j < tree[e].i){int temp = tree[e].j; tree[e].j=tree[e].i; tree[e].i=temp;} break;}
		else if(tree[e].i==node2 && tree[e].j==node1) {tree[e].i=node4; if(tree[e].j < tree[e].i){int temp = tree[e].j; tree[e].j=tree[e].i; tree[e].i=temp;} break;}
	}	
	for(register unsigned int e=0; e<tsize; e++){
		if(tree[e].i==node3 && tree[e].j==node4) {tree[e].j=node2; if(tree[e].j < tree[e].i){int temp = tree[e].j; tree[e].j=tree[e].i; tree[e].i=temp;} break;}
		else if(tree[e].i==node4 && tree[e].j==node3) {tree[e].i=node2; if(tree[e].j < tree[e].i){int temp = tree[e].j; tree[e].j=tree[e].i; tree[e].i=temp;} break;}
	}	
}


void Solver::NNI(int taxon,bool flag){
	// cout<<"In NNI..."<<endl;
	register unsigned int tsize=tree.size();
	vector<EDGE> temp_tree; 
	EDGE f; for(register unsigned int e=0; e<tsize; e++){f.i=tree[e].i; f.j=tree[e].j; temp_tree.push_back(f);}
	for(register unsigned int e=0; e<tsize; e++){
		//if e is internal
		if(tree[e].i >n && tree[e].j>n){
			register unsigned int node[5]; 
			register unsigned int c=1;
			for(register unsigned int k=1; k<=3; k++) if(M[tree[e].i-n][k]!=tree[e].j) {node[c]=M[tree[e].i-n][k]; c++;}
			for(register unsigned int k=1; k<=3; k++) if(M[tree[e].j-n][k]!=tree[e].i) {node[c]=M[tree[e].j-n][k]; c++;}
			NNISwap(tree[e].i,node[2],tree[e].j,node[3]);
			ComputeTopologicalMatrix(taxon);
			if(Optimum.tree_len > len){UpdateOptimum2(flag);} 
			tree.clear(); for(register unsigned int s=0; s<tsize; s++){f.i=temp_tree[s].i; f.j=temp_tree[s].j; tree.push_back(f);}
			NNISwap(tree[e].i,node[2],tree[e].j,node[4]);
			ComputeTopologicalMatrix(taxon);
			if(Optimum.tree_len > len){UpdateOptimum2(flag);} 
			tree.clear(); for(register unsigned int s=0; s<tsize; s++){f.i=temp_tree[s].i; f.j=temp_tree[s].j; tree.push_back(f);}
			ComputeTopologicalMatrix(taxon);
		}
	}	
}

/*END NNI BLOCK*/

void Solver::NJtree(){
	cout<<"Running the NJ algorithm...";
	double **NJdist;
    NJdist=new double*[2*n-1]; 
    for(int i=0; i<=2*n-2; i++) NJdist[i]=new double[2*n-1];
	NJdist[0][n]=n;NJdist[0][n+1]=0;
	for(register unsigned int i=1;i<=n-1;i++){
		NJdist[0][i]=i;
		for(register unsigned int j=i+1;j<=n;j++){
			NJdist[i][j]=dist[i][j];
			NJdist[j][i]=dist[j][i];
		}
	}
	tree.clear();

	struct EDGE e; for(register unsigned int k=1;k<=n;k++){e.i=k;e.j=n+1;tree.push_back(e);}
	
	// cout<<"printing star tree..."<<endl;
	// for(int e=0; e<tree.size(); e++) cout<<e<<": "<<tree[e].i<<" "<<tree[e].j<<endl;	

	struct MINPOS{double value; int a; int i; int b; int j;	int idx;};
	MINPOS min;
	for(register unsigned int k=4; k<=n; k++){
		min.value = INF;
		for(register unsigned int i=1; i<=n-k+3; i++){
			for(register unsigned int j=i+1; j<=n-k+4; j++){
 				if (NJdist[(int)NJdist[0][i]][(int)NJdist[0][j]] < min.value){
					min.value=NJdist[(int)NJdist[0][i]][(int)NJdist[0][j]];
					min.a=NJdist[0][i]; min.i=i;
					min.b=NJdist[0][j]; min.j=j;
				}
			}	
		}
		for(register unsigned int e=0; e<(int)tree.size(); e++)
			if ((tree[e].i == min.a && tree[e].j == n+1) || (tree[e].i == n+1 && tree[e].j == min.a)) {min.idx = e; break;}
		for(register unsigned int e=0; e<(int)tree.size(); e++)
			if ((tree[e].i == min.b && tree[e].j == n+1) || (tree[e].i == n+1 && tree[e].j == min.b)) {tree.erase(tree.begin()+e); break;}

		struct EDGE edge1, edge2, edge3; 	
		edge1.i=min.b;      edge1.j=n+k-2;
		edge2.i=tree[min.idx].i;  edge2.j=n+k-2;
		edge3.i=tree[min.idx].j;  edge3.j=n+k-2;

		tree.erase(tree.begin()+min.idx);
		tree.push_back(edge1);
		tree.push_back(edge2);
		tree.push_back(edge3);		

		// cout<<"printing current tree..."<<endl;
		// for(int e=0; e<tree.size(); e++) cout<<e<<": "<<tree[e].i<<" "<<tree[e].j<<endl;	

		for(register unsigned int i=min.i; i<=n-k+4; i++) NJdist[0][i]=NJdist[0][i+1];
		for(register unsigned int i=min.j-1; i<=n-k+4; i++) NJdist[0][i]=NJdist[0][i+1];
		NJdist[0][n-k+3]=n+k-2;
		for(register unsigned int i=1; i<=n-k+2; i++){
			NJdist[(int)NJdist[0][i]][n+k-2]=0.5*(NJdist[min.a][(int)NJdist[0][i]]+NJdist[min.b][(int)NJdist[0][i]]-min.value);
			NJdist[n+k-2][(int)NJdist[0][i]]=NJdist[(int)NJdist[0][i]][n+k-2];
		}
	}
	ComputeTopologicalMatrix(n);
	cout<<"done! BME length of the NJ-tree:   "<<"\x1b[92m"<<len<<"\x1b[0m"<<endl;

    for(int i=0; i<=2*n-2; i++) delete [] NJdist[i]; delete [] NJdist;	
}

void Solver::SWAA(){
	cout<<"Running the SWA algorithm...";
	tree.clear();
	SetStartingTree();
	ComputeTopologicalMatrix(3);
	for(register unsigned int k=4; k<=n; k++){ 
		double min=INF; 
		register unsigned int pos=0; 
		register unsigned int tsize=tree.size();
		for(register unsigned int e=0; e<tsize; e++){
			Add(k,e); 
		    ComputeTopologicalMatrix(k); 
		    // double LB=PardiLowerBound(k);
		    double LB=PardiLowerBound(k);
			if(min > LB){pos=e; min=LB;}
			Restore(e);			
	 	}
        Add(k,pos); 
        cout<<setprecision(PRECISION);
        //cout<<"Adding taxon "<<k<<" in edge "<<pos<<". Current length: "<<len<<endl;
	}
	ComputeTopologicalMatrix(n); 
	cout<<"done! BME length of the SWA tree: "<<"\x1b[92m"<<len<<"\x1b[0m"<<endl;
}

void Solver::SWAA2(){
	cout<<"Running the 2nd SWA algorithm...";
	tree.clear();
	SetStartingTree();
	ComputeTopologicalMatrix(3);
	for(register unsigned int k=4; k<=n; k++){ 
		double min=INF; 
		register unsigned int pos=0; 
		register unsigned int tsize=tree.size();
		for(register unsigned int e=0; e<tsize; e++){
			Add(k,e); 
		    ComputeTopologicalMatrix(k); 
		    // double LB=PardiLowerBound(k);
		    double LB=PardiGammaLowerBound(k);
			if(min > LB){pos=e; min=LB;}
			Restore(e);			
	 	}
        Add(k,pos); 
        cout<<setprecision(PRECISION);
        //cout<<"Adding taxon "<<k<<" in edge "<<pos<<". Current length: "<<len<<endl;
	}
	ComputeTopologicalMatrix(n); 
	cout<<"done! BME length of the 2nd SWA tree: "<<"\x1b[92m"<<len<<"\x1b[0m"<<endl;
}

void Solver::ShortcutTree_and_reorder_Taxa(){
	//Assumption: ComputeTopologicalMatrix(n) was already applied for the current 'tree'
	double **new_dist=new double*[n+1]; for(register unsigned int i=0; i<=n;i++) new_dist[i] = new double[n+1];
	vector<int> visited; visited.push_back(1);
	vector<int> unvisited; for(int i=2;i<=n;i++) unvisited.push_back(i);
	struct MINPOS{double value; int vertex; int idx;}; MINPOS min;
	int current_taxon=1;
	while(unvisited.size()){
		min.value = MAX_DOUBLE;
		for (register unsigned int v=1; v<=n; v++)
			if((current_taxon!=v)&&(Tau[current_taxon][v]<min.value))
				for(register unsigned int i=0; i<unvisited.size(); i++)
					if(unvisited[i]==v){min.value=Tau[current_taxon][v];min.vertex=v;min.idx=i;break;}
		current_taxon=min.vertex;
		visited.push_back(current_taxon);
		unvisited.erase(unvisited.begin()+min.idx);
	}
	for(register unsigned int i=0; i<visited.size(); i++) cout << visited[i] << endl;
	for(register unsigned int i=0; i<visited.size()-1; i++){
		new_dist[i+1][i+1]=0;
		for(register unsigned int j=i+1; j<visited.size(); j++){
			new_dist[i+1][j+1]=dist[visited[i]][visited[j]];
			new_dist[j+1][i+1]=new_dist[i+1][j+1];
		}
	}
	visited.clear();
	unvisited.clear();
	for(register unsigned int i=0; i<=n; i++) delete [] dist[i]; delete [] dist;
	dist=new_dist;
}

void Solver::SetPrimal_and_Taxa_Order(){
	cout<<"Computing some primal bounds for the instance..."<<endl;
	if(!readtree){	
		// Greedy SAS
		SWAA();
		StatsArray[0].RootUB=len; 
		Optimum.tree_len=len;
		Optimum.OverallNumber_of_Primal_Solutions++;
		Optimum.tree.clear();
		EDGE edge;
		register unsigned int tsize=tree.size();
		for(register unsigned int e=0; e<tsize; e++){
			edge.i=tree[e].i;
			edge.j=tree[e].j;
			Optimum.tree.push_back(edge);
		} 
		Optimum.primal_bound_list.push_back(Optimum.tree);
		double current_val=len; 
		if(torder==3) ShortcutTree_and_reorder_Taxa();
		// Greedy Agglo
		NJtree();
		if(Optimum.tree_len > len){
			UpdateOptimum2(false);
			cout<<"+The NJ-tree has a length shorter than the SWA's one. Updating primal bound..."<<endl;
		}else{
			cout<<"+Keeping the SWA tree as primal bound."<<endl;
		}
		if(torder==2) ShortcutTree_and_reorder_Taxa();
		//NNI...
		cout<<"Calculate all NNIs for the current primal bound..."<<endl;
		current_val = Optimum.tree_len;
		while(true){
			NNI(n,false);
			if(current_val != Optimum.tree_len){
				cout<<"+New NNI local optimum: "<<"\x1b[92m"<<Optimum.tree_len<<"\x1b[0m"<<endl;
				// Optimum.primal_bound_list.push_back(Optimum.tree);
				current_val=Optimum.tree_len;
			} else break;
		}
		cout<<"Done with NNI updates."<<endl;
		cout<<"Overall number of primal bound updates: "<<Optimum.OverallNumber_of_Primal_Solutions<<endl;
		
		if(torder==4) RandomizeTaxaOrder();
	} else {
		cout<<"Calculate all NNIs for the current primal bound..."<<endl;
		double current_val = Optimum.tree_len;
		while(true){
			NNI(n,false);
			if(current_val != Optimum.tree_len){
				cout<<"+New NNI local optimum: "<<"\x1b[92m"<<Optimum.tree_len<<"\x1b[0m"<<endl;
				// Optimum.primal_bound_list.push_back(Optimum.tree);
				current_val=Optimum.tree_len;
			} else break;
		}
		cout<<"Done with NNI updates."<<endl;
		cout<<"Overall number of primal bound updates: "<<Optimum.OverallNumber_of_Primal_Solutions<<endl;
		
		//Best FastME solution was read and stored as a phylogeny in Optimum.tree

		// Calculate the estimated evolutionary distance between all pairs of taxa
		ComputeWeightedOptimalTopologicalMatrix();
		// Reorder dist based on the estimated evolutionary distances Wtau (taxa order stored in pi)
		PhylogeneticDiversityTaxaOrder();
		// Reorder the primal solution
		tree.clear();
		struct EDGE edge; double input_length=0;
		for(register unsigned int e=0; e<(int)Optimum.tree.size(); e++){
			if(Optimum.tree[e].i<=n){
				for(register unsigned int k=1; k<=n;k++)if(pi[k]==Optimum.tree[e].i) edge.i=k; 
			}else edge.i=Optimum.tree[e].i;
			if(Optimum.tree[e].j<=n){
				for(register unsigned int k=1; k<=n;k++)if(pi[k]==Optimum.tree[e].j) edge.j=k;
			}else edge.j=Optimum.tree[e].j;
			edge.weight=Optimum.tree[e].weight*precomputedpow[n];
			input_length+=edge.weight;
			tree.push_back(edge);
		}
		delete [] pi;
		// Calculate the length of the phylogeny based on its topology and the distance matrix
		ComputeTopologicalMatrix(n);
		//cout << "input tree length = " << input_length << endl;
		//cout << "tree length derived from input topology = " << len << endl;
		Optimum.tree_len=len;
		Optimum.tree.clear();
		for(register unsigned int e=0; e<(int)tree.size(); e++){
			edge.i=tree[e].i;
			edge.j=tree[e].j;
			edge.weight=tree[e].weight;
			Optimum.tree.push_back(edge);
		}
		Optimum.OverallNumber_of_Primal_Solutions++;
		Optimum.primal_bound_list.push_back(Optimum.tree);
	}
}

void Solver::SetLower_Bounds(int boundtype){
	cout<<"Computing both Pardi's and the LP relaxation at the root node of the search tree...";

	// Pardi
	len=0; tree.clear();
	SetStartingTree();
	ComputeTopologicalMatrix(3);
	StatsArray[0].PardiRootLB=PardiLowerBound(3);       //Setting Pardi's root LB
	StatsArray[0].PardiGammaRootLB=PardiGammaLowerBound(3); 

	if(boundtype!=1){
		// Manifold
		problem2.setMsgLevel(0);
		XPRSsetintcontrol(problem2.getXPRSprob(), XPRS_THREADS,  OverallCoreNumber);
		XPRSsetintcontrol(problem2.getXPRSprob(), XPRS_LPTHREADS,OverallCoreNumber);
		XPRSsetintcontrol(problem2.getXPRSprob(), XPRS_BARTHREADS,OverallCoreNumber);

		problem2.minim("ld");
	
		XPRSsetintcontrol(problem2.getXPRSprob(), XPRS_THREADS,   1);
		XPRSsetintcontrol(problem2.getXPRSprob(), XPRS_LPTHREADS, 1);
		XPRSsetintcontrol(problem2.getXPRSprob(), XPRS_BARTHREADS,1);
		problem2.setMsgLevel(0);

		cout<<"done."<<endl;
		StatsArray[0].ManifoldRootLB=problem2.getObjVal(); //Setting Manifold root LB
		cout<<endl<<endl;  

		// Kraft
		StatsArray[0].KraftRootLB=KraftLowerBound();

		// LP
		problem.setMsgLevel(3);
		XPRSsetintcontrol(problem.getXPRSprob(), XPRS_THREADS,  OverallCoreNumber);
		XPRSsetintcontrol(problem.getXPRSprob(), XPRS_LPTHREADS,OverallCoreNumber);
		XPRSsetintcontrol(problem.getXPRSprob(), XPRS_BARTHREADS,OverallCoreNumber);

		problem.minim("ld");
	
		XPRSsetintcontrol(problem.getXPRSprob(), XPRS_THREADS,   1);
		XPRSsetintcontrol(problem.getXPRSprob(), XPRS_LPTHREADS, 1);
		XPRSsetintcontrol(problem.getXPRSprob(), XPRS_BARTHREADS,1);
		problem.setMsgLevel(0);

	/*for(register unsigned int i=1; i<=n-1; i++){
		cout << i << ": ";
		for(register unsigned int j=i+1; j<=n; j++){
			double sum = 0.0;
			for(register unsigned int k=2; k<=n-1; k++){
				cout << std::setprecision(6) << abs(path[i][j][k].getSol()) << " ";
			}
			cout << "| ";
		}
		cout << endl;
	}*/


		cout<<"done."<<endl;
		StatsArray[0].LPRootLB=problem.getObjVal(); //Setting LP root LB
		cout<<endl<<endl;  

		//PrimalDualLB();

	}
	tree.clear();
	SetStartingTree(); 	
	Push(0,3,3);
}



bool Solver::Push(int start, int end, int taxon){
    NODE node; 
	//Recopying partial tree
    node.partial_tree.clear();
	EDGE edge0;
	for(register unsigned int e=0; e<(int)tree.size(); e++){
		edge0.i=tree[e].i;
		edge0.j=tree[e].j;
		if(eorder==1) edge0.prio=tree[e].prio;
		node.partial_tree.push_back(edge0);
	}
	//Storing taxon, start edge and end edge
	node.start_edge=start;
    node.end_edge=end;
    node.taxon=taxon;
    node.empty=false;
    bool inserted_properly=false;
    #pragma omp critical
    {
		if(Queue.size()<OverallCoreNumber){
			switch(qrule){
				// FIFO
				case 0: node.prio=1; break;
				// LIFO
				case 1: if(Queue.size()) node.prio=Queue.top().prio+1; else node.prio=1; break;
				// non-increasing #edges
				case 2: node.prio=(int)tree.size(); break;
				// non-decreasing #edges
				case 3: node.prio=(-1)*((int)tree.size()); break;
				// random
				case 4: srand(time(NULL)); node.prio=rand()%100; break;
			}
			//cout << "Core "<<Core<<": "<< "push taxon "<<node.taxon<<" in ("<<node.start_edge<<","<<node.end_edge<<")"<<endl;
			Queue.push(node);
			inserted_properly=true;
		}	
	}
	return inserted_properly;
}


NODE Solver::RetrieveNode_and_Delete_from_Queue() {
	NODE node;
	node.empty=true; 
	node.partial_tree.clear();
    #pragma omp critical
    {	
      if(!Queue.empty()){
   		node = Queue.top();
   		Queue.pop();
 	   }
 	}
 	return node;
}


void Solver::SetNode_and_Run(NODE *node){
	int taxon = node->taxon;
	tree.clear();
	//Reset of the decoding matrices.
	for(register unsigned int i=0; i<=n; i++)    for(register unsigned int j=0; j<=n; j++) Tau[i][j]=0;
	for(register unsigned int i=0; i<=n-2;  i++) for(register unsigned int j=0; j<7;  j++) M[i][j]=0;
	for(register unsigned int i=0; i<n-1; i++) mem[i]=0;	

	struct EDGE edge;
	for(register unsigned int e=0; e<(int)node->partial_tree.size(); e++){
		edge.i=node->partial_tree[e].i;
		edge.j=node->partial_tree[e].j;
		edge.prio=node->partial_tree[e].prio;
		tree.push_back(edge);
	}	
	
	//clearing the tree of node 
	node->partial_tree.clear();
	Explore(node->start_edge,node->end_edge, node->taxon+1);
}
/*Stat rosa pristina nomine 
  Nomina nuda tenemus */