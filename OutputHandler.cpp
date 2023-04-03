class OutputHandler{
private:
	ofstream out;
	int **M;
	int **Tau;  
	int *mem;
	int n;
public:
	OutputHandler(int taxa){
		n=taxa;
		M=new int*[n-1];  for(register unsigned int i=0; i<=n-2;  i++) M[i]=new int[7]; 
		Tau=new int*[n+1];  for(register unsigned int i=0; i<=n;  i++)   Tau[i]=new int[n+1]; 
	  mem=new int[n-1];
	}
	~OutputHandler(){
	  for(register unsigned int i=0; i<=n-2; i++)  delete [] M[i]; delete [] M; 
	  for(register unsigned int i=0; i<=n; i++)    delete [] Tau[i]; delete [] Tau; 
	  delete [] mem; 
	}
	void PrintOptimalTree(vector <struct EDGE>tree){
		for(register unsigned int i=0; i<(int)tree.size(); i++){
			cout<<tree[i].i<<"\t"<<tree[i].j<<endl;
		}
	}
	
	void PrintOptimalTreeCVSFormat(vector <struct EDGE>tree){
		cout<<"Taxon,Taxon,Weight"<<endl;
		for(register unsigned int i=0; i<(int)tree.size(); i++){
			cout<<tree[i].i;
			cout<<",";
			cout<<tree[i].j;
			cout<<",";
			cout<<1;
			cout<<endl;
		}
	}

	void ComputeAdjacents(vector <struct EDGE>tree){
		int tsize=tree.size();
		for(register unsigned int i=1; i<=n-2;i++) mem[i]=1;
		for(register unsigned int e=0; e<tsize; e++){
			int temp1=tree[e].i; int temp2=tree[e].j;	
			if (temp1 > n) {M[temp1-n][mem[temp1-n]]=temp2; mem[temp1-n]++;}
			if (temp2 > n) {M[temp2-n][mem[temp2-n]]=temp1; mem[temp2-n]++;}
		}
	}	

	void RecursiveTraversal(int node, int father,int depth){
		for(int i=0; i<=depth; i++) out<<"  ";
		out<<"\"name\" : "<<node<<","<<endl;
		for(int i=0; i<=depth; i++) out<<"  ";
		out<<"\"parent\" : "; father == -1 ? out<<"\"null\"" : out<<father;
		if(node <= n){out<<endl; return;}
		else{
			out<<","<<endl;
			for(int i=0; i<=depth; i++) out<<"  ";
			out<<"\"children\" : ["<<endl;
			for(int i=1; i<=3; i++){
				if(M[node-n][i]!=father){	
					for(int i=0; i<=depth; i++) out<<"  ";
					out<<"  {"<<endl;
					RecursiveTraversal(M[node-n][i], node,depth+1);
					for(int i=0; i<=depth; i++) out<<"  ";
					out<<"  }"; 
					if(i<=2) out<<","<<endl; 
					else out<<endl;
				}	
			}	
			for(int i=0; i<=depth; i++) out<<"  ";
			out<<"]"<<endl;
		}
	}

	int** partition_leaves;
	int* n_leaves;
	int current_leafset;

	// get the leafset of the subtree rooted in 'node' not including 'father'
	void getLeafset(int node, int father){
		if(node<=n) partition_leaves[current_leafset][++n_leaves[current_leafset]]=node;
		else for(int child=1; child<=3; child++) if(M[node-n][child]!=father) getLeafset(M[node-n][child], node);
	}

	int scale_factor;

	// Assumption: leaf <= n
	double external_edge_weight(int leaf){
		bool found_father = false; int father=n;
		while(!found_father){father++; for(int child=1; child<=3; child++) if(M[father-n][child]==leaf) found_father=true;}
		// (leaf) --- (father) is an edge in the tree
		partition_leaves = new int*[3]; for(int i=1;i<=2;i++) partition_leaves[i]=new int[n];
		n_leaves = new int[3]; n_leaves[1]=0; n_leaves[2]=0;
		current_leafset=0; for(int child=1; child<=3; child++) if(M[father-n][child]!=leaf) {current_leafset++; getLeafset(M[father-n][child],father);}
		// Print partition
		/*cout << endl;
		for(int i=1;i<=n_leaves[1];i++) cout << partition_leaves[1][i] << ", "; cout << endl;
		for(int i=1;i<=n_leaves[2];i++) cout << partition_leaves[2][i] << ", "; cout << endl;
		cout << endl;*/
		//for(int i=1; i<=n; i++){for(int j=1; j<=n; j++) cout << dist[i][j] << ", "; cout << endl;}
		// calculate edge weight
		double weight=0.0;
		for(int i=1; i<=2; i++)for(int j=1;j<=n_leaves[i];j++) weight+=dist[leaf][partition_leaves[i][j]]*pow(2,scale_factor-Tau[leaf][partition_leaves[i][j]]+1);
		for(int i=1;i<=n_leaves[1];i++)for(int j=1;j<=n_leaves[2];j++) weight-=dist[partition_leaves[1][i]][partition_leaves[2][j]]*pow(2,scale_factor-Tau[partition_leaves[1][i]][partition_leaves[2][j]]+1);
		//cout << "edge weight for leaf " << leaf << ": " << weight << endl;
		delete [] n_leaves;
		delete [] partition_leaves[1]; delete [] partition_leaves[2]; delete [] partition_leaves;
		return weight;
	}

	double internal_edge_weight(int node1, int node2){
		bool edge_exists=false;
		for(int child=1; child<=3; child++) if(M[node1-n][child]==node2) edge_exists=true;
		if(edge_exists){
			// (node1) --- (node2) is an edge in the tree
			partition_leaves = new int*[5]; for(int i=1;i<=4;i++) partition_leaves[i]=new int[n];
			n_leaves = new int[5]; n_leaves[1]=0; n_leaves[2]=0; n_leaves[3]=0; n_leaves[4]=0;
			current_leafset=0;
			for(int child=1; child<=3; child++) if(M[node1-n][child]!=node2) {current_leafset++; getLeafset(M[node1-n][child],node1);}
			for(int child=1; child<=3; child++) if(M[node2-n][child]!=node1) {current_leafset++; getLeafset(M[node2-n][child],node2);}
			// Print partition
			/*cout << endl;
			cout << "Partition by edge (" << node1 << "," << node2 << "):" << endl;
			for(int i=1;i<=n_leaves[1];i++) cout << partition_leaves[1][i] << ", "; cout << endl;
			for(int i=1;i<=n_leaves[2];i++) cout << partition_leaves[2][i] << ", "; cout << endl;
			for(int i=1;i<=n_leaves[3];i++) cout << partition_leaves[3][i] << ", "; cout << endl;
			for(int i=1;i<=n_leaves[4];i++) cout << partition_leaves[4][i] << ", "; cout << endl;
			cout << endl;*/
			// calculate edge weight
			double weight=0.0;
			for(int i=1;i<=n_leaves[1];i++){
				for(int j=1;j<=n_leaves[3];j++) weight+=dist[partition_leaves[1][i]][partition_leaves[3][j]]*pow(2,scale_factor-Tau[partition_leaves[1][i]][partition_leaves[3][j]]+1);
				for(int j=1;j<=n_leaves[4];j++) weight+=dist[partition_leaves[1][i]][partition_leaves[4][j]]*pow(2,scale_factor-Tau[partition_leaves[1][i]][partition_leaves[4][j]]+1);
			}
			for(int i=1;i<=n_leaves[2];i++){
				for(int j=1;j<=n_leaves[3];j++) weight+=dist[partition_leaves[2][i]][partition_leaves[3][j]]*pow(2,scale_factor-Tau[partition_leaves[2][i]][partition_leaves[3][j]]+1);
				for(int j=1;j<=n_leaves[4];j++) weight+=dist[partition_leaves[2][i]][partition_leaves[4][j]]*pow(2,scale_factor-Tau[partition_leaves[2][i]][partition_leaves[4][j]]+1);
			}
			for(int i=1;i<=n_leaves[1];i++)for(int j=1;j<=n_leaves[2];j++) weight-=dist[partition_leaves[1][i]][partition_leaves[2][j]]*pow(2,scale_factor-Tau[partition_leaves[1][i]][partition_leaves[2][j]]+1);
			for(int i=1;i<=n_leaves[3];i++)for(int j=1;j<=n_leaves[4];j++) weight-=dist[partition_leaves[3][i]][partition_leaves[4][j]]*pow(2,scale_factor-Tau[partition_leaves[3][i]][partition_leaves[4][j]]+1);
			//cout << "edge weight for internal edge: " << weight << endl;
			return weight;
		}
	}

	void ComputeD3Hierarchy(vector <struct EDGE>tree, string namefile){
		cout<<"\nGenerating a visual representation of the optimal solution...";
		ComputeAdjacents(tree);
		// for(int i=1; i<=n-2; i++){
		// 	cout<<i+n<<":\t";
		// 	for(int j=1; j<=3; j++) cout<<M[i][j]<<"\t"; cout<<endl;
		// }
		out.open("treeData.txt");
		out<<"\n\n\n"<<endl;
		out<<"var treeData =[{"<<endl;
		RecursiveTraversal(2*n-2,-1,0);
		out<<"}];"<<endl;
		out.close();
		stringstream outputnamefile; 
		outputnamefile << "visual_output_"<<namefile<<"_"<<n<<".html";
		string command1 = "cat Resources/VisualizationScripts/00.txt > "+outputnamefile.str()+"; cat treeData.txt >> "+outputnamefile.str()+"; cat Resources/VisualizationScripts/01.txt >> "+outputnamefile.str();
		system(command1.c_str());
		cout<<"done"<<endl;
		cout<<"The graphical representation of the optimal phylogeny is ready!\nIf the browser does not open automatically, please click on the file "+outputnamefile.str()<<endl;
		system("rm treeData.txt");
		string command2 = "open "+outputnamefile.str(); 
		system(command2.c_str());
	}	


	double ComputeTopologicalMatrix(vector <struct EDGE>tree, int taxon){
		int tsize=tree.size();
		for(int i=0; i<=n; i++) for(int j=0; j<=n; j++) Tau[i][j]=0;
		for(register unsigned int i=1; i<=n-2;i++) mem[i]=1;
		for(register unsigned int e=0; e<tsize; e++){
			int temp1=tree[e].i; int temp2=tree[e].j;	
			if (temp1 > n) {M[temp1-n][mem[temp1-n]]=temp2; mem[temp1-n]++;}
			if (temp2 > n) {M[temp2-n][mem[temp2-n]]=temp1; mem[temp2-n]++;}
		}	
		register int CurrentNode,VisitedNode; 
		double len=0;
		Tau[0][0] = 0; // Modified 2/6/2009
		for(register unsigned int j=0; j<tsize; j++){
			if(tree[j].i<=taxon){
				int i = tree[j].i;
				int father=tree[j].j; 
				
				Tau[i][i]=0; // Modified 2/6/2009 - rimodificato il 29/9/2020
				Tau[i][0]=0; // modificato il 29/9/2020

				M[father-n][4]=1; M[father-n][5]=i; M[father-n][6]=0; 
				CurrentNode=father; 
				while(true){
					if (M[CurrentNode-n][6]<3){
						M[CurrentNode-n][6]++; VisitedNode=M[CurrentNode-n][M[CurrentNode-n][6]];
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
								 if(obj_rescaling) len+=dist[i][VisitedNode]*precomputedpow[n]/precomputedpowRes[Tau[i][VisitedNode]];              //pow(2.0,(double)Tau[i][VisitedNode]/n);
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
	    return len;	     	
	}




	void OutputPrimalBoundList(vector <vector <struct EDGE> > tree_list, string namefile){
		stringstream outputnamefile; 
		outputnamefile << namefile << "_primal_bound_list_" << n <<".txt";
		out.open((outputnamefile.str()).c_str());

		cout << tree_list.size() <<endl;
		for(int k=1; k<tree_list.size(); k++){
			vector <struct EDGE> T = tree_list[k]; 
			out<<"Printing primal solution "<<k<<"..."<<endl;
			double len=0.0;
			for(register unsigned int i=0; i<(int)T.size(); i++){
				out<<T[i].i<<"\t"<<T[i].j<<"\t"<<T[i].weight<<endl;
				len+=T[i].weight;
			}
			out<<"Input tree length = " << len <<"\n"<< endl;
			out<<"Printing the relative topological matrix..."<<endl;
			double val=ComputeTopologicalMatrix(T,n);
			for(int i=1; i<=n; i++){
				for(int j=1; j<=n; j++) out << Tau[i][j] <<" ";
				out<<endl;
			}
			out<<"length: "<<std::setprecision(PRECISION)<<val<<endl;
			out<<endl;			
		}
		out.close();

		scale_factor=n; 
		double bme_length=0.0;
		for(int i=1;i<=n;i++) bme_length+=external_edge_weight(i);
		for(int i=n+1;i<=2*n-2;i++)for(int j=i+1;j<=2*n-2;j++) bme_length+=internal_edge_weight(i,j);
		//bme_length*=pow(2,n);
		//cout << "BME length of the final primal solution = " << bme_length/4 << endl;
	}	

	void OutputSpectralAnalysis(vector <vector <struct EDGE> > tree_list, double **dist){
		cout<<setprecision(10)<<endl<<endl;		
		cout<<"Printing data for spectral analysis..."<<endl;
		cout<<endl;
		// cout<<"n = "<<n<<";"<<endl;
		cout<<"dist"<<n<<"={"<<endl;
		for(int i=1; i<=n; i++){
			cout<<"{";
			for(int j=1; j<=n-1; j++) cout<<dist[i][j]<<",";
			cout<<dist[i][n]<<"}";
			if(i<n) cout<<","<<endl; 	
		}
		cout<<"};"<<endl<<endl;


		for(int k=1; k<tree_list.size(); k++){
			vector <struct EDGE> T = tree_list[k]; 
			double val=ComputeTopologicalMatrix(T,n);
			cout<<"T"<<n<<"n"<<k<<"={"<<endl;
			for(int i=1; i<=n; i++){
				cout<<"{";
				for(int j=1; j<=n-1; j++) cout<<Tau[i][j]<<",";
				cout<<Tau[i][n]<<"}";
				if(i<n) cout<<","<<endl; 	
			}
			cout<<"};"<<endl<<endl;		
		}		
	}	


};


/*Stat rosa pristina nomine 
  Nomina nuda tenemus */