#include <map>

class InputHandler{
private:
	ifstream in;
	int* alpha_sigma;
	map<string,int> taxa_names;
public:
	InputHandler(){}
	~InputHandler(){}
	bool ActivateEntropyAnalysis;
	void ReadDistanceMatrix(char * const, string &, int, double** &, int &, int);
	void ReadPrimalSolution(char * const, string &, int &);
	void DoubleStochasticRescaling(int n, double **d, double** &dist);
	void Print(int **,int,int);
	void Print(double **,int,int);
	void PrintForMathematica(double **,int,int);
	void CallMosel(string, string, double** &, int);
	void CallConcorde(string, string, double** &, int);
	void InputSettler(int, char * const [], int &, bool &, bool &, int &, int &, bool &, int &, int &, int &, int &, int &);
	void HowTo();
	void StabilityChecker(int, double **);
	void MachineEpsilon();
	double EntropyLB(int, double, double * &);
	void EntropyAnalysis(int, double ** &, bool);
	void TriangleInequalityChecker(double **dist, int n);
};

void InputHandler::MachineEpsilon(){while ((1+EPS) != 1){TOLL = EPS; EPS /=2;} cout << "Machine epsilon: "<< TOLL << endl;} 
void InputHandler::StabilityChecker(int n, double **d){
	double m = INF;
	int acc=0;
	double smallest=INF;
	for(int i=1; i<=n-1; i++){
		for(int j=i+1; j<=n; j++){
			if(m > log2(d[i][j]/EPS) && log2(d[i][j]/EPS) >0){
				m=log2(d[i][j]/EPS);	
				// cout<<"d("<<i<<","<<j<<")="<<d[i][j]<<"\t"<<"EPS: "<<EPS<<"\tValue: "<<m<<endl;
			}
			if(d[i][j] < 1/pow(10.0,10.0)) acc++;
			if(smallest > d[i][j]) smallest=d[i][j];
		}
	}
	cout<<"The smallest entry in the matrix has value "<<smallest<<"."<<endl;
	cout<<"The matrix contains "<<acc<<" number of entries having value smaller than 10e-10."<<endl;
	cout<<"Numerical stability issues for the considered instance concerns path-lengths longer than or equal to "<<floor(m)<<"."<<endl;	
}

void InputHandler::DoubleStochasticRescaling(int n, double **d, double** &dist){
	cout<<"Running the double stochastic rescaling algorithm...";
	// d = 2^-dist
	for(register unsigned int i=1;i<=n;++i)for(register unsigned int j=1;j<=n;++j)if(i!=j)d[i][j]=pow(2.0,-dist[i][j]);
	double RescalingTOLL = MAX_TOLL;
	double** old_d;
	while(true){
		old_d=d;
		for(register unsigned int j=1;j<=n;++j){
			double sum=0.0;for(register unsigned int l=1;l<=n;++l)sum+=d[l][j];
			for(register unsigned int i=1;i<=n;++i)d[i][j]/=sum;}
		for(register unsigned int i=1;i<=n;++i){
			double sum=0.0;for(register unsigned int l=1;l<=n;++l)sum+=d[i][l];
			for(register unsigned int j=1;j<=n;++j)d[i][j]/=sum;}
		double sum=0.0;
		for(register unsigned int i=1;i<=n;++i)
			for(register unsigned int j=1;j<=n;++j)
				sum+=pow(old_d[i][j]-d[i][j],2);
		if(sum>RescalingTOLL) old_d=d; else break;
	}
	cout<<"done!"<<endl;
	for(register unsigned int i=1;i<=n-1;i++)
		for(register unsigned int j=i+1;j<=n;j++){
			d[i][j]=(-log(d[i][j])/log(2));
			d[j][i]=d[i][j];}	
	old_d=NULL;
}

double InputHandler::EntropyLB(int n, double alpha, double * &d){
	if(n==4){
		return (d[2]/2+(d[3]+d[4])/4-1.5*alpha);
	} else {
		// build vector d^alpha
		double* d_alpha = new double[n];
		for(int i=2; i<=n-2; i++) d_alpha[i]=d[i];
		double value=(d[n-1]+d[n])/2-alpha;
		// sort d^alpha
		int idx = -1; for(int i=2; i<=n-2; i++) if(d_alpha[i]>value){idx=i;break;}
		if(idx == -1) idx=n-1;
		for(int i=n-1; i>idx; i--) d_alpha[i] = d_alpha[i-1]; d_alpha[idx]=value;
		// continue recursion on d^alpha
		double bound = EntropyLB(n-1,alpha,d_alpha);
		delete [] d_alpha;
		// reverse the sorting according to d^alpha
		double sigma_value=alpha_sigma[idx];
		for(int i=idx; i<=n-2; i++) alpha_sigma[i]=alpha_sigma[i+1]; alpha_sigma[n-1]=sigma_value;
		// expand the path-length sequence
		alpha_sigma[n-1]++; alpha_sigma[n]=alpha_sigma[n-1];
		return bound;
	}
}

void InputHandler::EntropyAnalysis(int n, double ** &dist, bool DSrescaling){
	cout<<"\nRunning entropic analysis..."<<endl;	
	//Searching for the maximum entry of D. If the max d_ij is such that 2^-dij <= MAX_TOLL then no rescaling possible
	double maxentry = -1;
	for(int i=1; i<=n-1; i++) for(int j=i+1; j<=n; j++) if(maxentry < dist[i][j]) maxentry = dist[i][j]; 
	if(pow(2.0,-maxentry) > MAX_TOLL){
		double **d=new double*[n+1]; for(register unsigned int i=0; i<=n; i++) {d[i] = new double[n+1]; d[i][i]=0;}
		DoubleStochasticRescaling(n,d,dist);
		//now d is the double stochastic form of dist
		//cout<<"Printing rescaled matrix..."<<endl;
		//Print(d,n,n);
		//for(register unsigned int i=1; i<=n; i++){double sum=0.0; for(register unsigned int j=1; j<=n; j++) sum+=pow(2.0,-d[i][j]); cout << sum << endl;}
		//cout<<"Printing original matrix..."<<endl;
		//Print(dist,n,n);
		double alpha=0.0; double max_alpha = 1.0; double TOLL = 0.000001;
		double *d_alpha=new double[n+1]; 
		int *pi=new int[n+1];
		alpha_sigma=new int[n+1];
		int **sigma=new int*[n+1]; for(register unsigned int i=1; i<=n; i++) sigma[i] = new int[n+1];

		cout << "Running the auxillary function minimization...";
		while(alpha<=max_alpha){
			double entropybound=0.0;
			for(register unsigned int k=1; k<=n; k++){
				// pick k-th row of rescaled distance matrix d
				for(register unsigned int j=1; j<=n; j++){d_alpha[j]=dist[k][j]; pi[j]=j;}
				// sort k-th row of d
				for(int i=1; i<=n-1; i++){int min = i;
					for(int j=i+1; j<=n; j++) if(d_alpha[min]>d_alpha[j]) min=j;
					double value = d_alpha[min]; int idx = pi[min];
					while(min > i){d_alpha[min]=d_alpha[min-1]; pi[min]=pi[min-1]; min--;}
					d_alpha[i]=value; pi[i]=idx;}
				// check for submodularity guarantee
				double a = (d_alpha[n-3]+d_alpha[n-4])/2-d_alpha[n-5];
				if(a>max_alpha) max_alpha=a;
				// define the optimal solution for f_4
				alpha_sigma[1]=0; alpha_sigma[2]=1; alpha_sigma[3]=2; alpha_sigma[4]=2;
				// calculate the entropy lower bound
				entropybound += EntropyLB(n,alpha,d_alpha);
				// reverse sort
				for(register unsigned int j=1; j<=n; j++) sigma[k][pi[j]] = alpha_sigma[j];
			}
			// here shift K calculated using sigma. However, K can (in theory) be derived independently from sigma.
			K=0.0; for(int i=1; i<=n; i++)for(int j=1; j<=n; j++) K+=(dist[i][j]-d[i][j])/pow(2.0,sigma[i][j]+1);
			entropybound+=alpha*(3*n-6);
			entropybound/=2;
			if(entropybound>opt_entropybound){opt_entropybound=entropybound;opt_alpha=alpha;}
			alpha += TOLL;
		}
		cout<<"done!"<<endl;
		//cout << "lower bound for alpha (="<< opt_alpha <<") to ensure submodularity (unique solution): " << max_alpha << endl;
		delete [] pi; delete [] d_alpha; for(register unsigned int i=1; i<=n;i++) delete [] sigma[i]; delete [] sigma;

		double third_eq = 0; 
		for(int i=1; i<=n; i++)for(int j=i+1; j<=n; j++) third_eq += d[i][j]*pow(2.0,-d[i][j]); 	
		double ratio = third_eq*2/(2*n-3); 
		cout<<"The entropic ratio for the rescaled distance matrix is: "<<ratio<<endl;
		if(DSrescaling){
			cout<<"Keeping the rescaled matrix for the global search...";
			for(int i=0; i<=n; i++) for(int j=0; j<=n; j++) dist[i][j]=d[i][j];
			opt_entropybound-=K;
			opt_entropybound*=pow(2.0,n);
			cout<<"done."<<endl;
		}	
		else {cout<<"Assuming no rescaling of the input distance matrix. "<<endl;
			opt_entropybound*=pow(2.0,n);
		}	
		for(register unsigned int i=0; i<n+1;i++) delete [] d[i]; delete [] d;
	}		
	else cout<<"No rescaling possible for the considered matrix. Proceeding without rescaling. "<<endl;
	cout<<"Entropy analysis completed.\n"<<endl;	
}


void InputHandler::CallMosel(string namefile, string result, double** &dist, int n){
	 cout<<"\nPreparing data for Mosel..."<<endl;	
	 ofstream myfile; 
	 myfile.open(result+"Resources/MoselTSPSolver/input.txt");
  	 myfile << "N: "<<n<<endl;
  	 myfile << "dist:["<<endl;
  	 for(register unsigned int i=1; i<=n; i++){
  		for(register unsigned int j=1; j<=n; j++) myfile<<dist[i][j]<<" ";
  		myfile <<endl;
  	 }
  	 myfile << "]"<<endl;
 	 myfile.close();
  	 string cmd="mosel exec "+result+"Resources/MoselTSPSolver/TSP.mos DATAFILE="+result+"Resources/MoselTSPSolver/input.txt VERBOSE=1 SECPATH="+result+"Resources/MoselTSPSolver/";
  	 system(cmd.c_str());
	 ifstream in1; in1.open("TSPout.txt",std::ios::in); 
	 if(in1.is_open()){
		cout<<"Reading taxa order..."; 
		double **d=new double*[n+1]; for(register unsigned int i=1; i<=n; i++) d[i] = new double[n+1];
	   	for(register unsigned int i=1;i<=n;i++)for(register unsigned int j=1; j<=n; j++) d[i][j]=dist[i][j];
	   	int *olist=new int[n+1];
	    int t=1; 
	    while(in1 >> olist[t]) t++;
	    // cout<<"\nMosel optimum "<<endl; 
	    // for(register unsigned int i=1;i<=n;i++) cout<< olist[i]<<endl;
	    // cout<<endl;
	    if(t-1 == n){
	   		// for(register unsigned int i=1;i<=n;i++) in >> olist[i];
			for(register unsigned int i=1;i<=n;i++)for(register unsigned int j=1; j<=n; j++) dist[i][j]=d[olist[i]][olist[j]];
		}	
		else cout<<"Proceeding without circular order."<<endl;	
		for(register unsigned int i=1; i<n+1;i++) delete [] d[i]; delete [] d; delete[] olist;
		system("rm -f TSPout.txt");
		in1.close();
		cout<<"done!"<<endl;
	}
	else{
	 	cout<<"Proceeding without circular order."<<endl;
	}
}

void InputHandler::CallConcorde(string namefile, string result, double** &dist, int n){
	 cout<<"\nPreparing data for Concorde..."<<endl;	
	 ofstream myfile;
	 myfile.open (result+"Resources/MoselTSPSolver/input_concorde.txt");
  	 myfile<<"NAME: "<<namefile<<endl;
  	 myfile<<"TYPE: TSP"<<endl;
   	 myfile<<"DIMENSION: "<<n<<endl;
   	 myfile<<"EDGE_WEIGHT_TYPE: EXPLICIT"<<endl;
   	 myfile<<"EDGE_WEIGHT_FORMAT: FULL_MATRIX"<<endl;
   	 myfile<<"EDGE_WEIGHT_SECTION"<<endl;
   	 int max_precision=0;  
  	 for(register unsigned int i=1; i<=n; i++){
  		for(register unsigned int j=1; j<=n; j++){
  			stringstream s; 
  			s<<dist[i][j];
  			string s1 = s.str();
  			size_t pos = s1.find_last_of("."); 
  			string s2 = s1.substr(pos+1);
  			int len = s2.length();
  			if(len > max_precision) max_precision=len;
  		}
  	 }
  	 cout<<"Max precision: "<<max_precision<<endl;
  	 for(register unsigned int i=1; i<=n; i++){
  		for(register unsigned int j=1; j<=n; j++){
  			int num = round(dist[i][j]*pow(10,max_precision)); 
  			myfile<< num <<"\t";
  		}	
  		myfile <<endl;
  	 }  	 
	 myfile<<"EOF"<<endl;
	 myfile.close();
  	 string cmd="concorde "+result+"Resources/MoselTSPSolver/input_concorde.txt";
	 system(cmd.c_str());
	 system("rm -f Oinput_concorde.mas");
	 system("rm -f Oinput_concorde.pul");
	 system("rm -f Oinput_concorde.sav");
	 system("rm -f input_concorde.mas");
	 system("rm -f input_concorde.pul");
	 system("rm -f input_concorde.sav");
	 ifstream in1; in1.open("input_concorde.sol",std::ios::in); 
	 if(in1.is_open()){
		cout<<"Reading taxa order..."; 
		double **d=new double*[n+1]; for(register unsigned int i=1; i<=n; i++) d[i] = new double[n+1];
	   	for(register unsigned int i=1;i<=n;i++)for(register unsigned int j=1; j<=n; j++) d[i][j]=dist[i][j];
	   	int *olist=new int[n+1];
	    int t=0; 
	    while(in1 >> olist[t]){olist[t]++; t++;}
	    // cout<<"\nConcorde optimum "<<endl; 
	    // for(register unsigned int i=1;i<=n;i++) cout<< olist[i]<<endl;
	    // cout<<endl;
	    if(t-1 == n){
			for(register unsigned int i=1;i<=n;i++)for(register unsigned int j=1; j<=n; j++) dist[i][j]=d[olist[i]][olist[j]];
		}	
		else cout<<"Proceeding without circular order."<<endl;	
		for(register unsigned int i=1; i<n+1;i++) delete [] d[i]; delete [] d; delete[] olist;
		// system("rm -f input_concorde.sol");
		in1.close();
		cout<<"done!"<<endl;
	}
	else{
	 	cout<<"Proceeding without circular order."<<endl;
	}		 
}

void InputHandler::TriangleInequalityChecker(double **dist, int n){
	for(int i=1; i<=n; i++){
		for(int j=1; j<=n; j++){
			for(int k=1; k<=n; k++){
				if(i!=j && i!=k && j!=k){
					if( dist[i][j] > dist[i][k] + dist[k][j] ){
						cout<<"\x1b[91m"<<"failure: "<<"\x1b[0m"<<"d("<<i<<","<<j<<") = "<<dist[i][j]<<" > d("<<i<<","<<k<<") = "<<dist[i][k]<<" + d("<<k<<","<<j<<") = "<<dist[k][j]<<endl;
						return;
					}	
				}
			}
		}
	}
	cout<<"the property is satisfied."<<endl;
}


void InputHandler::ReadDistanceMatrix(char * const fullnamefile, string &namefile, int val, double** &dist, int &n, int readingPattern){ 
	if(n < 4){cout<<"\nWarning: the number of taxa must be a positive integer greater than or equal to 4. Please, try again."<<endl; exit(0);}
	in.open(fullnamefile,std::ios::in); 
	if(in.is_open()==false){cout << "Warning: Unable to open file. Please, try again."; exit(0);}
	string fullpath_namefile=fullnamefile;
	size_t pos = fullpath_namefile.find_last_of("/\\"); 
	namefile = fullpath_namefile.substr(pos+1);

	cout<<"Preprocessing file "<<"\x1b[92m"<<namefile<<"\x1b[0m...";

	int N; string line;
	if(readingPattern==1) in >> N;
	if(readingPattern==2){getline(in,line); ostringstream ss; ss << line; N=atoi(ss.str().c_str());}
	if(val==MAX_INT) n=N; else{if(n > N) n=N; else n=val;}	
	if(n<=N && n >=4){
		cout<<"\nTaxa in the current instance: \x1b[91m"<<n<<"\x1b[0m"<<endl;
		double **d=new double*[N];   for(register unsigned int i=0; i<N;i++)     d[i] = new double[N];
		      dist=new double*[n+1]; for(register unsigned int i=0; i<=n;i++) dist[i] = new double[n+1];
		if(readingPattern==1)
			for(register unsigned int i=0;i<N;i++) for(register unsigned int j=0; j<N; j++) in >> d[i][j];
		if(readingPattern==2)
			for(register unsigned int i=0;i<N;i++) {
				getline(in,line); ostringstream ss; ss << line;
				int p=0; while(ss.str()[p]!=' ') p++; //cout << ss.str()[p++]; cout << endl;
				taxa_names.insert(make_pair(ss.str().substr(0,p),i+1)); //cout << ss.str().substr(0,p) << endl;
				//cout << "insert " << ss.str().substr(0,p) << " at pos " << taxa_names[ss.str().substr(0,p)] << endl;
				for(register unsigned int j=0; j<N; j++){
					while(ss.str()[++p]==' ') continue;
					string entry; while((p<ss.str().length())&&(ss.str()[p]!=' ')) entry.append(1,ss.str()[p++]);
					d[i][j] = atof(entry.c_str());
				} 
			}
		for(register unsigned int i=1;i<=n;i++)for(register unsigned int j=1; j<=n; j++) dist[i][j]=d[i-1][j-1];
		for(register unsigned int i=0; i<N;i++) delete [] d[i]; delete [] d;
	}	
	else{
		 cout<<"\nWarning!\nThe specified input number of taxa must be a positive integer greater than or equal to 4.\nIf this number exceeds the number of taxa contained in the input file, it will be automatically set to the number of taxa in such a file.\nPlease, try again by keeping in mind this constaint."<<endl; exit(0); 
	}	
	in.close();	

	cout<<"Checking if the distance satisfies the triangle inequality property...";
	TriangleInequalityChecker(dist,n);
	
	if(torder==1){
		//trying to reading the associate order
		stringstream sm2; sm2<<namefile; string sm22=sm2.str(); 
		sm22.erase(sm22.end()-4, sm22.end()); 
		string sm33=sm22;
		sm22.erase(sm22.begin(), sm22.begin()+3); 
		stringstream sm3; sm3<<"instances/orders/"<<sm33<<n<<"order.dat"; sm22=sm3.str();
		cout<<"Searching for taxa order file "<<sm22<<"...";
		in.open(sm22.c_str(),std::ios::in); 
		if(in.is_open()){
			cout<<"File found! Reading taxa order..."; 
			double **d=new double*[n+1]; for(register unsigned int i=1; i<=n; i++) d[i] = new double[n+1];
    		for(register unsigned int i=1;i<=n;i++)for(register unsigned int j=1; j<=n; j++) d[i][j]=dist[i][j];
    		int *olist=new int[n+1];
    		for(register unsigned int i=1;i<=n;i++) in >> olist[i];
			for(register unsigned int i=1;i<=n;i++)for(register unsigned int j=1; j<=n; j++) dist[i][j]=d[olist[i]][olist[j]];
			for(register unsigned int i=1; i<n+1;i++) delete [] d[i]; delete [] d; delete[] olist;
			in.close();
			cout<<"done!"<<endl;
		} 
		else{
			cout<<"No taxa order found."<<endl;
		 	cout<<"Calling TSP Solver...";
		 	std::array<char, 128> buffer; 
	     	std::string result;
	     	std::unique_ptr<FILE, decltype(&pclose)> pipe(popen("pwd", "r"), pclose);
	     	if (!pipe){throw std::runtime_error("Absolute path finder failed!"); exit(1);}
	     	while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) result += buffer.data();
	     	result.erase(remove(result.begin(), result.end(), '\n'), result.end());
	     	result+="/";
	     	//result contains the absolute path
	   	 	switch(torder){
	   			case 1: CallMosel(namefile,result,dist,n); break;
	   			//case 2: CallConcorde(namefile,result,dist,n); break;
	  	 	}		 
		}
	}

	cout<<"The processing of "<<"\x1b[92m"<<namefile<<"\x1b[0m is completed."<<endl<<endl;
}

/* Assumptions: - ReadDistanceMatrix was applied (-> map taxa_names is well-defined)
                - fullnamefile is in newick format (with root taxon at the end of the file)
                - all taxa names in fullnamefile are indices of the map taxa_names */
void InputHandler::ReadPrimalSolution(char * const fullnamefile, string &namefile, int &n){ 
	if(n < 4){cout<<"\nWarning: the number of taxa must be a positive integer greater than or equal to 4. Please, try again."<<endl; exit(0);}
	string mod_fullnamefile = fullnamefile;
	mod_fullnamefile.replace(std::string(fullnamefile).length()-7,7,"mat.txt.nwk");
	in.open(mod_fullnamefile,std::ios::in); 
	if(in.is_open()==false){cout << "Warning: Unable to open file. Please, try again."; exit(0);}
	string fullpath_namefile=mod_fullnamefile;
	size_t pos = fullpath_namefile.find_last_of("/\\"); 
	namefile2 = fullpath_namefile.substr(pos+1);

	cout<<"Preprocessing file "<<"\x1b[92m"<<namefile2<<"\x1b[0m..."<<endl;

	double length=0.0;
	vector <struct EDGE> tree; 	
	string line; getline(in,line);
	ostringstream ss; ss << line;
	int idx = 0; int internal = n+1;
	do{
		// last entry opens a new edge
		// -> take last (-1,-1,-1) from tree and push back half an edge (taxon,-1,w)
		if((ss.str()[idx-1]=='(')&&(ss.str()[idx]!='(')){
			struct EDGE edge = tree.back(); tree.pop_back();
			int idx2=idx+1; while(ss.str()[idx2]!=':') idx2++;
			edge.i = taxa_names[ss.str().substr(idx,idx2-idx)];
			//cout << ss.str().substr(idx,idx2-idx) << " is at pos " << edge.i << endl;
			idx=idx2+1; while(ss.str()[idx2]!=',') idx2++;
			edge.weight = stod(ss.str().substr(idx,idx2-idx));
			tree.push_back(edge); idx=idx2+1;
		}
		// last entry created half an edge and new entry doesn't opens new edge
		// -> take last (A,-1,w) from tree and push back half an edge (B,-1,w'')
		//    while push back edges (A,B,w) and (taxon,B,w') to Optimum.tree
		if((ss.str()[idx-1]==',')&&(ss.str()[idx]!='(')){
			struct EDGE e1 = tree.back(); tree.pop_back(); e1.j = internal; 
			struct EDGE e2; e2.j = internal;
			int idx2=idx+1; while(ss.str()[idx2]!=':') idx2++;
			e2.i = taxa_names[ss.str().substr(idx,idx2-idx)];
			//cout << ss.str().substr(idx,idx2-idx) << " is at pos " << e2.i << endl;
			idx=idx2+1; while((ss.str()[idx2]!=')')&&(ss.str()[idx2]!=',')) idx2++; 
			// 2nd while-condition violated -> close to end of file
			e2.weight = stod(ss.str().substr(idx,idx2-idx));
			if(ss.str()[idx2+1]==';') e2.j=e1.i;
			else {Optimum.tree.push_back(e1); length += e1.weight;}
			Optimum.tree.push_back(e2); length += e2.weight;
			if(ss.str()[idx2]==',') idx=idx2; // read last "third of tree"
			// (A,B,w) and (taxon,B,w') created and in Optimum.tree
			// -> last entry in tree: (-1,-1,-1) or (X,-1,w)
			while((ss.str()[idx2]==')')&&(ss.str()[idx2+1]==':')){
				struct EDGE edge = tree.back(); tree.pop_back(); struct EDGE edge2; edge2.i=-1;
				if(edge.i != -1) edge2 = edge; // edge2 = (X,-1,w)
				edge.i = internal; edge.j=-1;edge.weight=-1; // edge = (B,-1,-1)
				idx2++; idx=idx2+1; while((ss.str()[idx2]!=',')&&(ss.str()[idx2]!=')')) idx2++;
				edge.weight = stod(ss.str().substr(idx,idx2-idx)); // edge = (internal,-1,w)
				internal++;
				if(edge2.i != -1){ // we have ss.str()[idx2]==')' before end of file
					if(ss.str()[idx2+1]==';'){
						edge2.j=internal-1; Optimum.tree.push_back(edge2); length += edge2.weight;
					} else {
						// we have two weighted half-edges -> create new internal vertex to connect them to
						edge.j = internal; edge2.j = internal;
						Optimum.tree.push_back(edge); Optimum.tree.push_back(edge2);
						length += edge.weight + edge2.weight;
					}
					if(ss.str()[idx2]==',') idx=idx2; // read last "third of tree"
				} else { // we have ss.str()[idx2]==','
					tree.push_back(edge); idx=idx2+1; 
				}
			}
		}
		if(ss.str()[idx]=='('){
			struct EDGE edge; edge.i = -1; edge.j = -1; edge.weight = -1;
			tree.push_back(edge); idx++;
		}
		if(ss.str()[idx]==','){
			struct EDGE edge; edge.i = internal; edge.j = -1; edge.weight = -1;
			tree.push_back(edge); idx++; internal++;
		}
	} while(tree.size()>=1);
	Optimum.tree_len = length;

	in.close();	

	if(Optimum.tree.size()!=(2*n-3)) cout<<"\x1b[91m Error: Read "<<Optimum.tree.size()<<" out of the required total of "<< (2*n-3)<<" edges.\x1b[0m"<<endl;
	
	//for(int i=0;i<Optimum.tree.size();i++)cout<<"("<<Optimum.tree[i].i<<","<<Optimum.tree[i].j<<")"<<endl;
	cout<<"The processing of "<<"\x1b[92m"<<namefile<<"\x1b[0m is completed."<<endl<<endl;
}

void InputHandler::Print(int **C,int row,int col){cout<<endl;for(register unsigned int i=1;i<=row;i++){for(register unsigned int j=1; j<=col;j++) cout << C[i][j]<< "\t"; cout <<endl;}cout<<endl;}
void InputHandler::Print(double **C,int row,int col){cout<<endl; for(register unsigned int i=1;i<=row;i++){for(register unsigned int j=1; j<=col;j++) cout << C[i][j]<< "\t"; cout <<endl;}cout<<endl;}

void InputHandler::PrintForMathematica(double **C,int row,int col){
	cout<<"Printing distance matrix in Mathematica format..."<<endl; 
	cout<<"Dist={";
	for(register unsigned int i=1;i<=row;i++){
		cout<<"{";
		for(register unsigned int j=1; j<=col;j++){
			if(j==col) cout << C[i][j]<< "}";
			else cout<<C[i][j]<<",";
		}
		if(i<row) cout<<","<<endl;
	}
	cout<<"};"<<endl;
}


void InputHandler::InputSettler(int argc, char * const argv[], int & taxa_to_consider, bool &DSrescaling, bool &readtree, int &cores_to_use, int &boundtype, bool &obj_rescaling, int &torder, int &eorder, int &qrule, int &n_p, int &s_p){
	ActivateEntropyAnalysis=false;
	if(argc>=2){
		for(register unsigned int i=2; i<argc; i++){
			string s=argv[i]; 
			std::size_t found=s.find("-newick=");
			if(found!=std::string::npos){
				if(s.length()<=8){
					cout<<"Syntax error with the parameter -newick"<<endl<<endl;
					HowTo();
				}
				if(s.substr(found+8)=="true" || s.substr(found+8)=="false"){
					if(s.substr(found+8)=="true") {readtree=true; readingPattern=2;}
				}
				else{cout<<"Bad value for parameter -newick. Aborting."<<endl; exit(0);}
			}
			else{ 
			found=s.find("-taxa=");
			if(found != std::string::npos){
				if(s.length()<=6){
			 		cout<<"Syntax error with the parameter -taxa"<<endl<<endl;
			  		HowTo();
			  	}	
			  	if(!readtree)
					try{taxa_to_consider = stoi(s.substr(found+6));}
					catch(std::exception& e){std::cout << "Error converting taxa value. Aborting." << endl; exit(0);}
				else
					cout<<"\x1b[91m Error: -newick=true, -taxa can't be specified.\x1b[0m"<<endl<<endl;
			} 
		    else{
		    	found=s.find("-boundtype=");
		    	if(found!=std::string::npos){
			    	if(s.length()<=11){
					   	cout<<"Syntax error with the parameter -boundtype"<<endl<<endl;
					   	HowTo();
					}	    		
					else if(s.substr(found+11)=="Pardi") boundtype=1;
						else if(s.substr(found+11)=="Res") obj_rescaling=true;
							 else{cout<<"Bad value for parameter -boundtype. Aborting."<<endl; exit(0);}
				} 
		    	else{
			    	found=s.find("-DSrescaling=");
			    	if(found!=std::string::npos){
			    		if(s.length()<=13){
					    	cout<<"Syntax error with the parameter -DSrescaling"<<endl<<endl;
					    	HowTo();
					    }
						if(s.substr(found+13)=="true" || s.substr(found+13)=="false"){
							if(s.substr(found+13)=="true") DSrescaling=true;
						}	
						else{cout<<"Bad value for parameter -DSrescaling. Aborting. "<< s.substr(found+11) <<endl; exit(0);}
					} 
	        	 	else{
				    	found=s.find("-cores=");
				    	if(found!=std::string::npos){
					    	if(s.length()<=7){
					    		cout<<"Syntax error with the parameter -cores"<<endl<<endl;
					    		HowTo();
					    	}
		    				try{cores_to_use = stoi(s.substr(found+7));}
							catch(std::exception& e){std::cout << "Error converting cores value. Aborting." << endl; exit(0);}
				    		if(cores_to_use<=0) cores_to_use=MAX_INT;
						} 
						else{
					    	found=s.find("-ordertaxa=");
					    	if(found!=std::string::npos){
					    		if(s.length()<=11){
					    			cout<<"Syntax error with the parameter -ordertaxa"<<endl<<endl;
					    			HowTo();
					    		} 
			    				try{torder = stoi(s.substr(found+11));}
								catch(std::exception& e){std::cout << "Error converting order value. Aborting." << endl; exit(0);}
								if((torder<0)||(torder>4)){torder=0; std::cout << "Reset to default." << endl;}
							} 
		        	 		else{
		        	 			found=s.find("-orderedges=");
								if(found!=std::string::npos){
									if(s.length()<=12){
					    			cout<<"Syntax error with the parameter -orderedges"<<endl<<endl;
					    			HowTo();
					    			} 
			    					try{eorder = stoi(s.substr(found+12));}
									catch(std::exception& e){std::cout << "Error converting order value. Aborting." << endl; exit(0);}
									if((eorder<0)||(eorder>2)){eorder=0; std::cout << "Reset to default." << endl;}
								}
								else{
									found=s.find("-queueselect=");
									if(found!=std::string::npos){
					    				if(s.length()<=13){
					    					cout<<"Syntax error with the parameter -queueselect"<<endl<<endl;
					    					HowTo();
					    				} 
			    						try{qrule = stoi(s.substr(found+13));}
										catch(std::exception& e){std::cout << "Error converting the queue selection value. Aborting." << endl; exit(0);}
										if((qrule<0)||(qrule>4)){qrule=0; std::cout << "Reset to default." << endl;}
									} 
									else{ 
		        	 					found=s.find("-numpartition=");
				    					if(found!=std::string::npos){
					    					if(s.length()<=14){
					    						cout<<"Syntax error with the parameter -numpartition"<<endl<<endl;
					    						HowTo();
					    					}
		    								try{n_p = stoi(s.substr(found+14));}
											catch(std::exception& e){std::cout << "Error converting cores value. Aborting." << endl; exit(0);}
				    						if(n_p<=1) n_p=2;
										} 
										else{
											found=s.find("-minpartition=");
				    						if(found!=std::string::npos){
					    						if(s.length()<=14){
					    							cout<<"Syntax error with the parameter -minpartition"<<endl<<endl;
					    							HowTo();
					    						}
		    									try{s_p = stoi(s.substr(found+14));}
												catch(std::exception& e){std::cout << "Error converting cores value. Aborting." << endl; exit(0);}
				    							if(s_p<=0) s_p=1;
											}
											else{ 
												if(strcmp(argv[i],"-p")==0){
		        	 								output_for_experiment=true;
		        	 							}
		        	 							else{
		        	 								found=s.find("-entropyanalysis=");
							    					if(found!=std::string::npos){
							    						if(s.length()<=17){
									    					cout<<"Syntax error with the parameter -entropanalysis"<<endl<<endl;
									    					HowTo();
									    				}
									    				if(s.substr(found+17)=="true" || s.substr(found+17)=="false"){
															if(s.substr(found+17)=="true") ActivateEntropyAnalysis=true;	
														} 
														else{cout<<"Bad value for parameter -entropyanalysis. Aborting."<<endl; exit(0);}
													}
														else{
															found=s.find("-maxruntime=");
															if(found!=std::string::npos){
								    							if(s.length()<=12){
										    						cout<<"Syntax error with the parameter -maxruntime"<<endl<<endl;
										    						HowTo();
										    					}
																try{
																	double temp = stoi(s.substr(found+12));
																	if (temp > 0) user_max_time = temp;
																}
																catch(std::exception& e){std::cout << "Error converting taxa value. Aborting." << endl; exit(0);}
															}
														}
													
												}		
											}
										}
									}
		        	 			}
		        	 		}
	        	 		}
	        	 	}	
	        	}		
	    	}} 
	    }
	}
}

void::InputHandler::HowTo(){
	cout<<"Synopsis: BME-Solver distance_file [optional parameters]"<<endl<<endl;
	cout<<"List of optional parameters"<<endl;
	cout<<"-taxa=[an integer K]"
		<<left<<setw(17)<<"\t"<<"Enables the processing of the first K taxa in an instance of N >= K taxa.\n"
		<<left<<setw(33)<<"\t"<<"For example, the option -taxa 10 consider the first 10 taxa in the considered instance;\n"
		<<left<<setw(33)<<"\t"<<"The option requires that K>=4; if K exceeds N then the solver sets automatically K=N.\n"
		<<left<<setw(33)<<"\t"<<"Default value: all taxa in the given instance. Deactivated if -newick=true."
		<<endl
		<<endl
	    <<"-boundtype=[Pardi, LP/Lagrangian, Res]"
	    <<left<<setw(1) <<"\t"<<"Specifies the type of lower bound to use during the search.\n"
	    <<left<<setw(33)<<"\t"<<"For example, the option -boundtype=Pardi enables the use of Pardi's bound.\n"
	    <<left<<setw(33)<<"\t"<<"The option -boundtype=LP/Lagrangian enables the use of the linear programming relaxation.\n"
	    <<left<<setw(33)<<"\t"<<"The option -boundtype=Res enable the use of the linear programming with rescaling of the objective function.\n"
	    <<left<<setw(33)<<"\t"<<"Default value: LP/Lagrangian."
	    <<endl
	    <<endl
	    <<"-DSrescaling=[true or false]"
	    <<left<<setw(9) <<"\t"<<"Enables or disables the double-stochastic rescaling of the distance matrix.\n"
	    <<left<<setw(33)<<"\t"<<"For example, the option -DSrescaling=true activates the rescaling; -DSrescaling=false deactivates it.\n"
	    <<left<<setw(33)<<"\t"<<"Default value: false."
	    <<endl
	    <<endl
	    <<"-cores=[a positive integer]"
	    <<left<<setw(9) <<"\t"<<"Specifies the number of computing cores to use.\n"
	    <<left<<setw(33)<<"\t"<<"For example, the option -cores 8 enables allows the use of 8 parallel threads.\n"
	    <<left<<setw(33)<<"\t"<<"Default value: all of the computing cores available."
	    <<endl
	    <<endl
	    <<"-ordertaxa=[0, 1, 2, 3, 4]"
	    <<left<<setw(9) <<"\t"<<"Specifies the type of circular order to use.\n"
	    <<left<<setw(33)<<"\t"<<"0 is none; 1 - TSP; 2 - NJ; 3 - Greedy SAS; 4 - random\n"
	    <<left<<setw(33)<<"\t"<<"Default value: 0. Deactivated if -newick=true."		
	    <<endl
	    <<endl
	    <<"-orderedges=[0, 1, 2]"
	    <<left<<setw(17) <<"\t"<<"Specifies the order of edges for branchings.\n"
	    <<left<<setw(33)<<"\t"<<"0 - insertion time; 1 - lower bound; 2 - random\n"
	    <<left<<setw(33)<<"\t"<<"Default value: 0."	    		    	    
	    // <<"-p"		  <<left<<setw(15)<<"\t"<<"print output for experiments";
	    <<endl
	    <<endl
	    <<"-queueselect=[0, 1, 2, 3, 4]"
	    <<left<<setw(9) <<"\t"<<"Specifies the type of queue order to use.\n"
	    <<left<<setw(33)<<"\t"<<"0 - FIFO; 1 - LIFO; 2 - Max. #edges; 3 - Min. #edges; 4 - random\n"
	    <<left<<setw(33)<<"\t"<<"Default value: 3."		
	    <<endl
	    <<endl
	    <<"-numpartition=[a positive integer]"
	    <<left<<setw(1) <<"\t"<<"Specifies the number of partition sets of edges to use to create jobs for the queue.\n"
	    <<left<<setw(33)<<"\t"<<"For example, the option -numpartition=5 partitions edges of a job into 5 sets of similar size.\n"
	    <<left<<setw(33)<<"\t"<<"In case the chosen number of partition sets exceeds the number of edges, the latter defines the number of partition sets.\n"		
	    <<left<<setw(33)<<"\t"<<"Default value: 3."		
	    <<endl
	    <<endl
	    <<"-minpartition=[a positive integer]"
	    <<left<<setw(1) <<"\t"<<"Specifies the minimum size of partition sets of edges to use to create jobs for the queue.\n"
	    <<left<<setw(33)<<"\t"<<"For example, the option -minpartition=3 partitions edges of a job into sets of similar size such that no individual set contains less than 3 edges.\n"
	    <<left<<setw(33)<<"\t"<<"In case the chosen minimum size of partition sets exceeds half the number of edges, the latter defines the minimum size of partition sets.\n"		
	    <<left<<setw(33)<<"\t"<<"Default value: 2."		
	    <<endl
	    <<endl
	    // <<"-distort=[true or false]"
	    // <<left<<setw(9) <<"\t"<<"Enables or disables the distortion of the distance matrix.\n"
	    // <<left<<setw(33)<<"\t"<<"For example, the option -distort=true activates the distortion; -rescaling=false deactivates it.\n"
	    // <<left<<setw(33)<<"\t"<<"Default value: false."	    		    	    
	    // <<"-p"		  <<left<<setw(15)<<"\t"<<"print output for experiments";	
	    <<"-entropyanalysis=[true or false]"
	    <<left<<setw(1) <<"\t"<<"Enables or disables the entropic analysis of the input distance matrix.\n"
	    <<left<<setw(33)<<"\t"<<"Default value: -DSrescaling. Can't be disabled if -DSrescaling=true."	    		    	    
	    // <<"-p"		  <<left<<setw(15)<<"\t"<<"print output for experiments";
	    <<endl
	    <<endl
	    <<"-newick=[true or false]"
	    <<left<<setw(17) <<"\t"<<"Enables or disables the reading of a phylogeny in newick format.\n"
	    <<left<<setw(33)<<"\t"<<"Default value: false. When set to true assumes that the distance matrix is read from [..]mat.txt and the phylogeny from [..]mat.txt.nwk\n"	    		    	    
	    <<left<<setw(33)<<"\t"<<"The taxa order is calculated from the phylogenetic diversity scores of the given phylogeny."	    		    	    
	    // <<"-p"		  <<left<<setw(15)<<"\t"<<"print output for experiments";
	    <<endl
	    <<endl
	    <<"-maxruntime=[a positive value]"
	    <<left<<setw(9) <<"\t"<<"Specifies the time limit for the global search.\n"
	    <<left<<setw(33)<<"\t"<<"For example -maxruntime=300 sets the time limit to 5 minutes.\n"
	    <<left<<setw(33)<<"\t"<<"Default value: -1 (no time limit)."		
	    <<endl
	    <<endl;	     
		exit(0);
}
/*Stat rosa pristina nomine 
  Nomina nuda tenemus */