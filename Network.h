/** Last updated 2010-02-03 15:36 **/
/** Fermi updating + Epidemic Spreading SIR model */

#ifndef _Network_h
#define _Network_h

#include<iostream>
#include<fstream>
#include<stdlib.h>
#include "PublicFunction.h"
#include<math.h>
using namespace std;

/** Cooperator、Defector矩阵 **/
const int C=0;//vaccination

const int D=1;//non-vaccination

const int OP_VAC =0; //0 no vaccination opinion update, =1 vaccination opinion update

const int Iseeds = 30; // number of infected individuals

int NVac = 1250;//Number of initially vaccinated individuals

int NVacB = 1250;//Number of initially vaccine believers

const double EPITIME = 10000;

double TR = 0.46;//disease transmission rate

double RR = 1.0/3.0;//recovery rate MAKE SURE R_0 is right

double mu =0.01; //strategy mutation
double Time=0;

double ef=0.8; //1-ef vaccine efficiency

double theta = 0;//additional cost for nonbelivers

double c=0.05;//cost of vaccination

double bias = 0.1; //bias of not vaccinated. [-.5, +.5]

double beta=1; //selection intensity


/** 输出实时错误并退出 */
void handleError(char errorText[])
{
	cout<<"... run-time error ..."<<endl;
	cout<<errorText;
	cout<<"... now exiting to system ..."<<endl;
	exit(1);
}

/***** 通用的网络 。开始*****/
class Network
{
public:
	/** 构造函数，参数n为节点总数 */
	Network(long n);
	~Network();

	/** 初始化邻居表（生成边） */
	virtual void initNeighbors(void)=0;
	
	/** 初始化所有节点的状态、Payoff */
	void initStates(void);

	/** 计算某节点的新状态 */
	inline int getNewState(long x);

	void StraUpdate(void);

	void EpiUpdate(ostream &out);

	void VDegree(ostream &out);//correlation between vaccination behavior with degree

	void snap(ostream &out);//output the system snapshot

	inline long getNumberOfNodes(void) { return NumberOfNodes; }

	inline long getNumberOfCs(void) { return NumberOfCs; } //number of vaccinated individuals

	inline long getNumberOfIs(void) { return NumberOfIs; }//number of infected individuals

	inline long getNumberOfSs(void) { return NumberOfSs; }//number of Susceptiable individuals

	inline long getNumberOfRs(void) { return NumberOfRs; }//number of Recovered individuals

	inline long getNumberOfVs(void){ return NumberOfVs; }// number of individuals of failed vaccination
	inline long getNumberOfVi(void){ return NumberOfVi; }// number of individuals of failed vaccination
	inline long getNumberOfVr(void){ return NumberOfVr; }// number of individuals of failed vaccination

	inline long getNumberOfV0(void){ return NumberOfV0; }//number of vaccinated believers
	inline long getNumberOfV1(void){ return NumberOfV1; }//number of vaccinated nonbelievers
	inline long getNumberOfV2(void){ return NumberOfV2; }//number of unvaccinated believers
	inline long getNumberOfV3(void){ return NumberOfV3; }//number of unvaccinated nonbelievers
	
	
	
	 friend ostream& operator << (ostream& out, const Network& net);
	
protected:

	long NumberOfNodes;//totoal number of nodes

	long NumberOfV0; //0-0 vaccinated, believers
	long NumberOfV1;//0-1 vaccinated, nonbelievers
	long NumberOfV2;// 1-0 unvaccianted believers
	long NumberOfV3;//1-1 unvaccinated nonbelievers

	long NumberOfCs; //number of vaccinated indiviudals  ::state 3
	long NumberOfBs;//number of believers
	long NumberOfVs; //number of vaccinated indiviudals  ::state 3
	long NumberOfVi;//number of infected vaccinated individuals:: state 4
	long NumberOfVr;//number of recoverred and vaccinated indiviudals ::state 5
	long NumberOfIs;//Number of infected susecptibles :: state 0
	long NumberOfSs;//Number of Susceptiable ones :: state 1
	long NumberOfRs;//Number of Recovery ones:: state 2
//	long NumberOfFs; //number of individuals with failed vaccination

	/** Neighbors[0] ~ Neighbors[NumberOfNeighbors[i]-1]表示节点i的所有邻居 */
	long ** Neighbors;
	
	/** NumberOfNeighbors[i]表示节点i的所有邻居的个数 */
	long * NumberOfNeighbors;

	/** 节点的状态，State0为当前时刻，State1为下一时刻 */
	int * State0, * State1;// 0: Vaccinated 1: NV 

	int * VacB0, * VacB1; //0 current, 1 updated vaccination opinon, 0, believers, 1, nonbelievers

	int * Dstate; // 0: Susceptiable; 1: Infected; 2: Recovery; 3: Vaccinated

	/** 节点的Payoff */
	double * Payoff;

	double * AveragePayoff;

};

Network::Network(long n)
{	
	NumberOfNodes=n;
	
	if(NumberOfNodes<=0) {
		Neighbors=0;
		NumberOfNeighbors=0;
		State0=0;
		Dstate=0;
		State1=0;
		VacB0=0;
		VacB1=0;
		Payoff=0;
		AveragePayoff=0;
		return;
	}
	
	Neighbors=new long * [NumberOfNodes];
	if(Neighbors==0) handleError("Neighbors==0");
	
	NumberOfNeighbors=new long [NumberOfNodes];
	if(NumberOfNeighbors==0) handleError("NumberOfNeighbors==0");
	
	State0=new int [NumberOfNodes];
	if(State0==0) handleError("State0==0");

	Dstate=new int [NumberOfNodes];
	if(Dstate==0) handleError("DState==0");
	
	State1=new int [NumberOfNodes];
	if(State1==0) handleError("State1==0");

	VacB0=new int [NumberOfNodes];
	if(VacB0==0) handleError("VacB0==0");

	VacB1=new int [NumberOfNodes];
	if(VacB1==0) handleError("VacB1==0");
	
	Payoff=new double [NumberOfNodes];
	if(Payoff==0) handleError("Payoff==0");
		
			
	AveragePayoff=new double [NumberOfNodes];
	if(AveragePayoff==0) handleError("AveragePayoff==0");
	
}

Network::~Network()
{
	if(Neighbors!=0) {
		for(int i=0;i<NumberOfNodes;i++) {
			if(Neighbors[i]!=0) {
				delete[] Neighbors[i];
				Neighbors[i]=0;
			}
		}
		delete[] Neighbors;
		Neighbors=0;
	}
	
	if(NumberOfNeighbors!=0) {
		delete[] NumberOfNeighbors;
		NumberOfNeighbors=0;
	}
	
	if(State0!=0) {
		delete[] State0;
		State0=0;
	}
	if(Dstate!=0) {
		delete[] Dstate;
		Dstate=0;
	}

	if(State1!=0) {
		delete[] State1;
		State1=0;
	}

	if(VacB0!=0) {
		delete[] VacB0;
		VacB0=0;
	}

	if(VacB1!=0) {
		delete[] VacB1;
		VacB1=0;
	}

	if(Payoff!=0) {
		delete[] Payoff;
		Payoff=0;
	}
	
		
	if(AveragePayoff!=0) {
		delete[] AveragePayoff;
		AveragePayoff=0;
	}

	NumberOfNodes=0;
	NumberOfCs=0; //number of vaccinated indiviudals  ::state 3
	NumberOfVs=0;//number of infected vaccinated individuals:: state 4
	NumberOfVi=0;//number of recoverred and vaccinated indiviudals ::state 5
	NumberOfVr=0;
	NumberOfIs=0;//Number of infected susecptibles :: state 0
	NumberOfSs=0;//Number of Susceptiable ones :: state 1
	 NumberOfRs=0;
	
}

/** 初始化所有节点的状态、Payoff */
void Network::initStates(void)
{
	long * node=new long [NumberOfNodes];
	if(node==0) handleError("node==0");

	long * Bnode=new long [NumberOfNodes];
	if(Bnode==0) handleError("Bnode==0");
	
	long i, nC2, nV2;
	for(i=0;i<NumberOfNodes;i++) {
		State0[i]=0;
		VacB0[i]=0;
		node[i]=i;
		Bnode[i]=i;
	}
	
	NumberOfCs=NumberOfNodes;
	NumberOfBs=NumberOfNodes;
	
	nC2=NVac;//Number of vaccinated individuals
	nV2=NVacB;//Number of believers
	
	while(NumberOfCs!=nC2) {
		i=rand()%NumberOfCs;
		State0[node[i]]=1;
		node[i]=node[--NumberOfCs];
	}
	while(NumberOfBs!=nV2) {
		i=rand()%NumberOfBs;
		VacB0[Bnode[i]]=1;
		Bnode[i]=Bnode[--NumberOfBs];
	}
	
NumberOfV0=0;
NumberOfV1=0;
NumberOfV2=0;
NumberOfV3=0;


	for(i=0;i<NumberOfNodes;i++)
		{	
		
		
			if(State0[i]==0)
				{
					Payoff[i] = -c;
					Dstate[i] = 3;

					if(VacB0[i]==1)
						{
						Payoff[i] -= theta;

						}
					/*if(getRand01()<ef)
						{
							Dstate[i] =0;
							NumberOfFs++;
						}
					else Dstate[i] = 3;*/
				}
			else {Payoff[i] = 0; Dstate[i] = 0;}
			if(State0[i]==0 && VacB0[i]==0)
				NumberOfV0++;
			if(State0[i]==0 && VacB0[i]==1)
				NumberOfV1++;
			if(State0[i]==1 && VacB0[i]==0)
				NumberOfV2++;
			if(State0[i]==1 && VacB0[i]==1)
				NumberOfV3++;

                    //if(State0[i]!=0 &&State0[i]!=1)
					//	cout<<State0[i]<<endl;
					//if(VacB0[i]!=0 &&VacB0[i]!=1)
					//	cout<<VacB0[i]<<endl;
			
				
		}

		//if(NumberOfV0+NumberOfV1+NumberOfV2+NumberOfV3<400)
		//cout<<222222222<<"\t"<<NumberOfV0+NumberOfV1+NumberOfV2+NumberOfV3<<endl;
	
	//NumberOfSs = NumberOfNodes - NumberOfCs - Iseeds + NumberOfFs;
	//NumberOfIs =Iseeds;
	//NumberOfRs = 0;

	//NumberOfCs; //number of vaccinated indiviudals  ::state 3
	NumberOfCs=NumberOfV0+NumberOfV1;
	NumberOfVs=NumberOfCs;
	//NumberOfV0=NumberOfBs;
	NumberOfVi=0;//number of infected vaccinated individuals:: state 4
	NumberOfVr=0;//number of recoverred and vaccinated indiviudals ::state 5
	NumberOfIs=0;//Number of infected susecptibles :: state 0
	NumberOfSs=NumberOfNodes - NumberOfVs;//Number of Susceptiable ones :: state 1
	 NumberOfRs=0;
	//cout<<NumberOfV0+NumberOfV2<<endl;

	 
	delete[] node;
}


/** Mutation in strategy is considered **/
/** Bias in taking vaccine is also considered **/
int Network::getNewState(long x) {
	long y=Neighbors[x][getRand(NumberOfNeighbors[x])];
	double p=1+exp((Payoff[x]-Payoff[y])*beta);
	p=1/p;

	if(getRand01()>p) return State0[x];
			return State0[y];
/*	
	if (getRand01()< mu)
		{
		if (getRand01()< 0.5+bias) return D;
		else return C;
		}
	else {    	
			if(getRand01()>p) return State0[x];
			return State0[y];
		}
*/


}

/** 更新所有节点的状态、Payoff，并计算最新的C的个数 */

void Network::StraUpdate(void)
{
	long i;

	if(OP_VAC==0) //no vaccination opinon update
		{

			for(i=0;i<NumberOfNodes;i++)
		{
			State1[i] = getNewState(i);
			//NumberOfCs += State0[i]-State1[i];
		}

	int * temp;
	temp=State0;
	State0=State1;
	State1=temp;
	
	//calculating the payoff
	//NumberOfFs=0;
	for(i=0;i<NumberOfNodes;i++)
		{	
			if(State0[i]==0) 
				{ 
					Payoff[i] = -c; 
					//if (getRand01()<ef) { Dstate[i] = 0; NumberOfFs++;} // Vaccine effectiveness is considered, w.p. ef vaccination fails.
					//else Dstate[i] = 3;
					 Dstate[i] = 3;
					if(VacB0[i]==1)
						{
						Payoff[i] -= theta;

						}
				}
			else{ Payoff[i] = 0; Dstate[i] = 0;}
		}

		
		}

	if(OP_VAC==1)//vaccination opinon updated simutaneously with vaccination decisions
		{
			for(i=0;i<NumberOfNodes;i++)
				{

				long y=Neighbors[i][getRand(NumberOfNeighbors[i])];
				double p=1+exp((Payoff[i]-Payoff[y])*beta);
				p=1/p;

					//cout<<Payoff[y]<<endl;

													
					//if(VacB0[y]!=0 &&VacB0[y]!=1)
						//cout<<VacB0[y]<<endl;

				if(getRand01()<p)
					{

					State1[i] = State0[y];
					VacB1[i] = VacB0[y];
					//NumberOfCs += State0[i]-State1[i];
					//NumberOfBs += VacB0[i]-VacB1[i];
					
					}
				else 
					{

					State1[i] = State0[i];
					VacB1[i] = VacB0[i];
					//NumberOfCs += State0[i]-State1[i];
					//NumberOfBs += VacB0[i]-VacB1[i];
					

					}

				

				}

			
for(i=0;i<NumberOfNodes;i++)
				{

				State0[i]=State1[i];
				VacB0[i] =VacB1[i];


}
				
			/*int * temp;
			temp=State0;
			State0=State1;
			State1=temp;

			int * temp1;
			temp1=VacB0;
			VacB0=VacB1;
			VacB1=temp1;
			*/

			//calculating the payoff
	//NumberOfFs=0;
	for(i=0;i<NumberOfNodes;i++)
		{	
			if(State0[i]==0) 
				{ 
					Payoff[i] = -c; 
					//if (getRand01()<ef) { Dstate[i] = 0; NumberOfFs++;} // Vaccine effectiveness is considered, w.p. ef vaccination fails.
					//else Dstate[i] = 3;
					 Dstate[i] = 3;
					if(VacB0[i]==1)
						{
						Payoff[i] -= theta;

						}
				}
			else{ Payoff[i] = 0; Dstate[i] = 0;}
		}



		}
	

	
	
	
	//delete[] temp;

	NumberOfV0=0;
	NumberOfV1=0;
	NumberOfV2=0;
	NumberOfV3=0;
		for(i=0;i<NumberOfNodes;i++)
		{
					if(State0[i]==0 && VacB0[i]==0)
						NumberOfV0++;
					if(State0[i]==0 && VacB0[i]==1)
						NumberOfV1++;
					if(State0[i]==1 && VacB0[i]==0)
						NumberOfV2++;
					if(State0[i]==1 && VacB0[i]==1)
						NumberOfV3++;
				//if(State0[i]!=0 &&State0[i]!=1)
				//		cout<<State0[i]<<endl;
				//if(VacB0[i]!=0 &&VacB0[i]!=1)
				//	cout<<VacB0[i]<<endl;
					
					
		}

      NumberOfCs=NumberOfV0+NumberOfV1;
	NumberOfVs=NumberOfCs;
	NumberOfVi=0;//number of infected vaccinated individuals:: state 4
	NumberOfVr=0;//number of recoverred and vaccinated indiviudals ::state 5
	NumberOfIs=0;//Number of infected susecptibles :: state 0
	NumberOfSs=NumberOfNodes - NumberOfVs;//Number of Susceptiable ones :: state 1
	 NumberOfRs=0;
	//if(NumberOfV0+NumberOfV1+NumberOfV2+NumberOfV3<400)
	//cout<<33333333<<"\t"<<NumberOfV0+NumberOfV1+NumberOfV2+NumberOfV3<<endl;
		//cout<<NumberOfV0+NumberOfV1<<endl;
		//cout<<NumberOfV0+NumberOfV1+NumberOfV2+NumberOfV3<<endl;

}


void Network::EpiUpdate(ostream &out)
{	
	NumberOfVs=NumberOfCs;
	NumberOfVi=0;//number of infected vaccinated individuals:: state 4
	NumberOfVr=0;//number of recoverred and vaccinated indiviudals ::state 5
	NumberOfIs=0;//Number of infected susecptibles :: state 0
	NumberOfSs=NumberOfNodes - NumberOfVs;//Number of Susceptiable ones :: state 1
	 NumberOfRs=0;
	long * Inode =  new long [NumberOfNodes]; //for infected individuals book keeping
	if(Inode==0) handleError("Inode==0");
	long * Cnode = new long [NumberOfNodes]; // for individuals' states that can change
	if(Cnode==0) handleError("Cnode==0");

	long i;
	long i1,Cn;//Cn is the total number of S+I
	i1=0;
	for(i=0;i<NumberOfNodes;i++)
		{

					Inode[i1]=i;
					Cnode[i1]=i;
					i1++;
			//if(Dstate[i]!=2 && Dstate[i]!=5)
				//{
					//Inode[i1]=i;
					//Cnode[i1]=i;
					//i1++;
				//}

		}
	Cn=i1;
	long i2;
	i2=0;
	//cout<<33<<endl;
	while(i2<Iseeds) {

		//i=getRand(i1+1);
		i=getRand(i1);
		if (Dstate[Inode[i]]==0)
			{
		Dstate[Inode[i]]=1;//selected for infected seeds
		NumberOfSs--;
		NumberOfIs++;
		i2++;
		Payoff[Inode[i]]-=1;//payoff for infection
		//Inode[i]=Inode[i1--];
		Inode[i] = Inode[--i1];}
				if (Dstate[Inode[i]]==3)
			{
		Dstate[Inode[i]]=4;//selected for infected seeds
		i2++;
		NumberOfVs--;
		NumberOfVi++;
		Payoff[Inode[i]]-=1;//payoff for infection
		//Inode[i]=Inode[i1--];
		Inode[i] = Inode[--i1];}
		
		
		}

	//cout<<34<<endl;
	//cout<<NumberOfIs<<endl;
	double time1;
	time1=0;
	double * Crate = new double [NumberOfNodes];
	double * CrateArea = new double [Cn];
	for (i=0;i<NumberOfNodes;i++) Crate[i]=0;

	for(i=0;i<NumberOfNodes;i++)
		{
			if(Dstate[i]==0)
				{   
					for(i2=0;i2<NumberOfNeighbors[i];i2++)
						{ 
							if(Dstate[Neighbors[i][i2]]==1||Dstate[Neighbors[i][i2]]==4)
								Crate[i] += TR;
						}
					
				}
					if(Dstate[i]==3)
				{   
					for(i2=0;i2<NumberOfNeighbors[i];i2++)
						{ 
							if(Dstate[Neighbors[i][i2]]==1||Dstate[Neighbors[i][i2]]==4)
								Crate[i] += TR*(1-ef);
						}
					
				}
			if(Dstate[i]==1||Dstate[i]==4) Crate[i]=RR;	
		}
	//for(i=0;i<NumberOfNodes;i++) cout<<Crate[i]<<endl;
	while(time1<EPITIME)
		{
			CrateArea[0]=Crate[Cnode[0]];
			for(i=1;i<Cn;i++)
				CrateArea[i] = CrateArea[i-1] + Crate[Cnode[i]];

			if(CrateArea[Cn-1]<0) cout<<CrateArea[Cn-1]<<endl;
			time1+= - log(1-getRand01())/CrateArea[Cn-1];
			Time+= - log(1-getRand01())/CrateArea[Cn-1];
			//cout<<22<<endl;
			i2=binSearch(CrateArea, Cn, getRand01()*CrateArea[Cn-1]);
			if(i2<0) i2=-i2-1;
			//cout<<23<<endl;
			switch(Dstate[Cnode[i2]])
				{
				cout<<100<<endl;
					case 0: //S->I
						Dstate[Cnode[i2]] = 1;
						Payoff[Cnode[i2]] -= 1;
						NumberOfSs--;
						NumberOfIs++;
						Crate[Cnode[i2]]=RR;
						for(i=0;i<NumberOfNeighbors[Cnode[i2]];i++)
							{
								if(Dstate[Neighbors[Cnode[i2]][i]]==0)
								Crate[Neighbors[Cnode[i2]][i]] += TR;
								if(Dstate[Neighbors[Cnode[i2]][i]]==3)
								Crate[Neighbors[Cnode[i2]][i]] += TR*(1-ef);
							}
						//cout<<i2<<"\t"<<Crate[Neighbors[Cnode[i2]][i]] <<"S->I"<<NumberOfIs<<endl;
						break;
					case 1: //I->R
						Dstate[Cnode[i2]]=2;
						//cout<<Cn<<"Cn"<<endl;
						for(i=0;i<NumberOfNeighbors[Cnode[i2]];i++)
							{
								if(Dstate[Neighbors[Cnode[i2]][i]]==0)
										Crate[Neighbors[Cnode[i2]][i]] -= TR;
								if(Dstate[Neighbors[Cnode[i2]][i]]==3)
										Crate[Neighbors[Cnode[i2]][i]] -= TR*(1-ef);
								//if(Crate[Neighbors[Cnode[i2]][i]]<0) cout<< Neighbors[Cnode[i2]][i]<<"neg"<<Crate[Neighbors[Cnode[i2]][i]]<<"I"<<NumberOfIs<<endl;
							}
						NumberOfIs--;
						NumberOfRs++;
						Cnode[i2]=Cnode[--Cn];
						//cout<<i2<<"I->R"<<NumberOfIs<<endl;
						break;
						case 3: //V->Iv
						Dstate[Cnode[i2]]=4;
						Payoff[Cnode[i2]] -= 1;
						Crate[Cnode[i2]]=RR;
												for(i=0;i<NumberOfNeighbors[Cnode[i2]];i++)
							{
								if(Dstate[Neighbors[Cnode[i2]][i]]==0)
										Crate[Neighbors[Cnode[i2]][i]] += TR;
								if(Dstate[Neighbors[Cnode[i2]][i]]==3)
										Crate[Neighbors[Cnode[i2]][i]] += TR*(1-ef);
								//if(Crate[Neighbors[Cnode[i2]][i]]<0) cout<< Neighbors[Cnode[i2]][i]<<"neg"<<Crate[Neighbors[Cnode[i2]][i]]<<"I"<<NumberOfIs<<endl;
							}
												NumberOfVs--;
												NumberOfVi++;
												break;
						case 4:
							Dstate[Cnode[i2]]=5;
								for(i=0;i<NumberOfNeighbors[Cnode[i2]];i++)
							{
								if(Dstate[Neighbors[Cnode[i2]][i]]==0)
										Crate[Neighbors[Cnode[i2]][i]] -= TR;
								if(Dstate[Neighbors[Cnode[i2]][i]]==3)
										Crate[Neighbors[Cnode[i2]][i]] -= TR*(1-ef);
								//if(Crate[Neighbors[Cnode[i2]][i]]<0) cout<< Neighbors[Cnode[i2]][i]<<"neg"<<Crate[Neighbors[Cnode[i2]][i]]<<"I"<<NumberOfIs<<endl;
							}
												NumberOfVi--;
												NumberOfVr++;
											
													Cnode[i2]=Cnode[--Cn];
														break;
							
						
						

				}
			//for(i=0;i<Cn;i++) cout<<Crate[i]<<" ";
		
			//out<<time1<<"\t"<<NumberOfVs<<"\t"<<NumberOfVi<<"\t"<<NumberOfVr<<"\t"<<NumberOfSs<<"\t"<<NumberOfIs<<"\t"<<NumberOfRs<<endl;
			if(NumberOfIs==0 && NumberOfVi==0) break;
			
		
		}
	//out<<Time<<"\t"<<NumberOfCs<<"\t"<<NumberOfSs<<"\t"<<NumberOfIs<<"\t"<<NumberOfRs<<"\t"<<NumberOfFs<<endl;
	delete[] Inode;
	delete[] Cnode;
	delete[] Crate;
	delete[] CrateArea;
	//cout<<time1<<endl;
}


void Network::snap(ostream &snapshot)
{	
	int i, j, M1,N1;
	M1=N1=100;
	for(i=0;i<M1;i++){
		for(j=0;j<N1;j++){
			snapshot<<i<<"\t"<<j<<"\t"<<Dstate[i*N1+j]<<endl;
			}
		}

}
void Network::VDegree(ostream &out)
{
	long *ks, *pk,*vk;
	ks = new long [NumberOfNodes];
	if(ks==0) handleError("ks==0");
	pk = new long [NumberOfNodes];
	if(pk==0) handleError("pk==0");
	vk = new long [NumberOfNodes];
	if(vk==0) handleError("vk==0");
	
	long i,j,k;
		for(i=0;i<NumberOfNodes;i++)
	{
		ks[i]=0;
		pk[i]=0;
		vk[i]=0;

	}

	for(j=0;j<10000;j++)
		{
		StraUpdate();

	for(i=0;i<NumberOfNodes;i++)
	{	
		ks[NumberOfNeighbors[i]]++;
		if(State0[i]==C) pk[NumberOfNeighbors[i]]++;
		for(k=0;k<NumberOfNeighbors[i];k++)
			if(State0[Neighbors[i][k]]==C) vk[NumberOfNeighbors[i]]++;
	}

		}
	out<<"Degree\tNumber\tVacLevel\tfracVN"<<endl;
	for(i=0;i<NumberOfNodes;i++)
		if(ks[i]!=0)
		out<<i<<"\t"<<ks[i]/10000<<"\t"<<pk[i]/(ks[i]+Zero)<<"\t"<<(vk[i]+Zero)/i/ks[i]<<endl;
	delete[] ks;
	delete[] pk;
	
}

/** 输出网络 */
ostream& operator << (ostream& out, const Network& net)
{
	out<<net.NumberOfNodes<<endl;
	long i,j;
	for(i=0;i<net.NumberOfNodes;i++) {
		out<<i<<"\t"<<net.NumberOfNeighbors[i]<<endl;
		for(j=0;j<net.NumberOfNeighbors[i];j++)
			out<<net.Neighbors[i][j]<<" ";
		out<<endl;
	}
	return out;
}

/***** 通用的网络 。结束*****/

#endif
