#include <iostream> 
#include <fstream> 
#include <string> 
#include <vector> 
#include <cmath> 
#include <list>
#include <sstream>
#include <random>

using namespace std;

class systemOfAgents{

	private:
			double lamda;
			double mu;
			double c;
			double T;
			double alpha;
			int N;
			double gamma;
			double Beta;
			double fV;
			double fI;
			
			std::uniform_real_distribution<> dis;
			std::mt19937 gen;
			std::random_device rd;
			
			double countV;
			double countI;
			double countIv;
			
			vector<vector<int> > EdgesG;
			vector<vector<int> > neighborsG;
			
			vector<vector<int> > EdgesE;
			
			
			vector<int> V;
			vector<double> win;
			vector<int> Inf;
			vector<int> NewHealth;
			
			void initEdgesG(){
				ifstream f("EdgesG.txt");
				int NEdges;
				f >> NEdges;
		
				vector<int> e(2);
				
				for (int i = 0; i < NEdges; i++){
					for (int j = 0; j < 2; j++){
						f >> e[j];
					}
					EdgesG.push_back(e);
				}
				f.close();
			}
			
			void initEdgesE(){
				ifstream f("EdgesE.txt");
				int NEdges;
				f >> NEdges;
		
				vector<int> e(2);
				
				for (int i = 0; i < NEdges; i++){
					for (int j = 0; j < 2; j++){
						f >> e[j];
					}
					EdgesE.push_back(e);
				}
				f.close();
			}
			
			void initNeighborsG(){
				ifstream inputfile("NeighborsG.txt");
				string s;
	
				while (getline(inputfile,s)){
					string a(s);
					std::stringstream stream(a);
					vector<int> d;
					int b;
					while (stream >> b) {
						d.push_back(b);
					}
					
					neighborsG.push_back(d);
				}
				inputfile.close();
				N = neighborsG.size();
				
				vector<int> A(N);
				V = A;
				Inf = A;
				
				vector<double> B(N);
				win = B;
			}
			
			void initParameters(){
				vector<double> Parameters(10);
				ifstream para("Parameters.txt");
	
				for (int i = 0; i < 10; i++){
					para >> Parameters[i];
				}
				
				para.close();
				
				fV = Parameters[0];
				T = Parameters[1];
				c = Parameters[2];
				mu = Parameters[3];
				fI = Parameters[4];
				gamma = Parameters[5];
				lamda = Parameters[6];
				Beta = Parameters[7];
				alpha = Parameters[8];
			}
			
			void initDistribution(){
				std::mt19937 a(rd());
				gen = a;
				std::uniform_real_distribution<> b(0, 1);
				dis = b;
			}
			
			void updatePayoff(){
				for (int i = 0; i < N; i++){
					if(V[i] == 1){
						win[i] += -c;
						if(neighborsG[i][0] != 0){
							for (int j = 0;j<neighborsG[i].size();j++){
								win[i] += -T*Inf[neighborsG[i][j]-1]*gamma*(1-alpha)/neighborsG[i].size();
							}
							win[i] += -T*countI/N*alpha*gamma;
						}
					}
					else{
						if(neighborsG[i][0] != 0){
							for (int j = 0;j<neighborsG[i].size();j++){
								win[i] += (-T*Inf[neighborsG[i][j]-1]*V[neighborsG[i][j]-1]*gamma-T*Inf[neighborsG[i][j]-1]*(1-V[neighborsG[i][j]-1]))*(1-alpha)/neighborsG[i].size();
							}
							win[i] += (-T*gamma*countIv/N-T*(countI-countIv)/N)*alpha;
						}
					}
					
				}
			}
			
			void initVaccination(){
				countV = 0.0;

				for (int i = 0; i < N; i++){
					if(dis(gen) < fV){
						V[i] = 1;
						countV += 1;
					}
					else{
						V[i] = 0;
					}
					win[i] = 0;
				}
			}
			
			void initializeInfections(){
				countI = 0.0;
				countIv = 0.0;
				vector<int> A(N);
				Inf = A;
				NewHealth = A;
				for(int i = 0; i < N; i++){
					if(V[i] == 0 && dis(gen) < fI){
						Inf[i] = 1;
						NewHealth[i] = 1;
						countI += 1;
					}
					else if(V[i] == 1 && dis(gen) < gamma*fI){
						Inf[i] = 1;
						NewHealth[i] = 1;
						countI += 1;
						countIv += 1;
					}
					else{
						Inf[i] = 0;
						NewHealth[i] = 0;
					}
				}

			}
			
			void updateStrategies(){
				vector<int> newStrategies(V);

				int m = 0;
				for (int i = 0; i < N; i++){
					if(neighborsG[i][0] != 0){
						uniform_int_distribution<int> uni(0,neighborsG[i].size()-1);
						m = neighborsG[i][uni(gen)]-1;
						
						if(dis(gen) <= 1.0/(1.0+exp(-Beta*(win[m]-win[i])))){
							newStrategies[i] = V[m];
						}
					}
				}
				countV = 0.0;
				for(int i = 0; i < N; i++){
					V[i] = newStrategies[i];
					win[i] = 0;
					countV += V[i];
				}
				if(countV/N > 1.0){
					cout << "Problem" << endl;
					cout << countV << endl;
				}
			}
					
		public:
			void init(){
				initParameters();
				initNeighborsG();
				initEdgesG();
				initEdgesE();
				initDistribution();
				initVaccination();
				initializeInfections();
			}
			
			double getV(){
				return countV/N;
			}
			
			void updateGame(){
				updatePayoff();
				updateStrategies();
			}
			
			double getI(){
				return countI/N;
			}
			
			void updateHealth(){
				for(int i = 0; i < EdgesE.size(); i++){
					if(Inf[EdgesE[i][0]] == 1 && Inf[EdgesE[i][1]] == 0 && NewHealth[EdgesE[i][1]] == 0){
						if(V[EdgesE[i][1]] == 0 && V[EdgesE[i][0]]==0){
							if(dis(gen) < lamda){
								NewHealth[EdgesE[i][1]] = 1;
							}	
						}
						else{
							if(dis(gen) < lamda*gamma){
								NewHealth[EdgesE[i][1]] = 1;
							}
						}
					}
					else if(Inf[EdgesE[i][0]] == 0 && Inf[EdgesE[i][1]] == 1 && NewHealth[EdgesE[i][0]] == 0){
						if(V[EdgesE[i][1]] == 0 && V[EdgesE[i][0]]==0){
							if(dis(gen) < lamda){
								NewHealth[EdgesE[i][0]] = 1;
							}	
						}
						else{
							if(dis(gen) < lamda*gamma){
								NewHealth[EdgesE[i][0]] = 1;
							}
						}
					}
				}
				for(int i = 0; i < N; i++){
					if(Inf[i] == 1 && dis(gen) < mu){
						NewHealth[i] = 0;
					}
				}
				countI = 0;
				countIv = 0;
				for(int i = 0; i < N; i++){
					Inf[i] = NewHealth[i];
					countI += Inf[i];
					if(V[i] == 1){
						countIv += Inf[i];
					}
				}
				if(countI/N > 1.0){
					cout << "Problem I" << endl;
					cout << countI << endl;
				}
				if(countIv > countI){
					cout << "Problem Iv" << endl;
					cout << countIv << endl;
					cout << countI << endl;
				}
			}
			
			int getN(){
				return N;
			}	
};

class simulateSystemOfAgents{
		private:
			systemOfAgents sys;
			int simulationStepsGame;
			vector<double> V;
			vector<double> Inf;
			
			void updateGame(){
				sys.updateGame();
				V.push_back(sys.getV());
			}
				
			void simulateGame(){
				for(int i = 0; i < simulationStepsGame; i++){
					updateGame();
				}
			}
			
			void updateHealth(){
				sys.updateHealth();
				Inf.push_back(sys.getI());
			}
			
			double getVMean(){
				double VTotalEnd(0.0);
				for (int i = simulationStepsGame/2; i < simulationStepsGame+1; i++){
					VTotalEnd += V[i];
				}
				
				return VTotalEnd/(simulationStepsGame+1-simulationStepsGame/2);
			}
			
			double getIMean(){
				double ITotalEnd(0.0);
				for (int i = simulationStepsGame/2; i < simulationStepsGame+1; i++){
					ITotalEnd += Inf[i];
				}
				
				return ITotalEnd/(simulationStepsGame+1-simulationStepsGame/2);
			}
				
		public:
			void init(){
				simulationStepsGame = 2000;
				sys.init();
				vector<double> A;
				V = A;
				vector<double> B;
				Inf = B;
				V.push_back(sys.getV());
				Inf.push_back(sys.getI());
			}
			
			void simulate(){
				for(int i = 0; i < simulationStepsGame; i++){
					updateGame();
					updateHealth();
				}
			}
			
			void writeOutput(){
				ofstream myfile;
				myfile.open ("output.txt");
				double a(getVMean());
				double b(getIMean());
				/*
				if(a > 1 or b > 1){
					cout << "a: " << a << endl;
					cout << "b: " << b << endl;
				}
				if(a<0.00000001 and b>0.00000001){
					myfile << "0.0" << " " << b << endl;
				}
				else if(a>0.00000001 and b<0.00000001){
					myfile <<  a << " " << "0.0" << endl;
				}
				else if(a<0.00000001 and b<0.00000001){
					myfile << "0.0" << " " << "0.0" << endl;
				}
				else{
					
				}*/
				myfile << fixed << a << " " << b << endl;
				myfile.close();
			}
		
};

int main(int argc,char* argv[]) {
	simulateSystemOfAgents A;
	//A.init();
	//A.simulate();
	//A.writeOutput();
	return 0;
}
