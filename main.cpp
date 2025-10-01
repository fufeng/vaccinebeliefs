/** 更新于：2008-09-02 12:20 */

#include "Lattice4.h" // 采用何种网络
#include <math.h>  
int main(int argc, char ** argv)
{
	Lattice *network; // full graph

	int M=50;//length
	int N=50; // width
	/** 总的演化次数：先演化TE0次，再演化TE1-TE0次 */
	int TE0=3000, TE1=4000;

	int TN=1; // 生成网络的总次数
	int TS=50; // 初始化节点状态的总次数

	/** 上述三个总次数的计数值 */
	long te, tn, ts;
	
	time_t timeNow;
	struct tm tmNow;
	
	double rowB,rowB1, rowB2, rowB3, rowB4;
	ofstream outRowB,outStat;
	time(&timeNow);

	srand((unsigned) time(&timeNow));

	outRowB.open("rowB.txt");
	if(outRowB==0) handleError("outRowB==0");
	outRowB<<"beta\tc\tvac\tVB\tVNB\tNVB\tNVNB"<<endl;

	outStat.open("snapshot.txt");
	if(outStat==0) handleError("outStat==0");
	

	/*outVDegree.open("Degree.txt");
	if(outVDegree==0) handleError("outVDegree==0");*/
	
	/*outRowT.open("Time.txt");
	if(outRowT==0) handleError("outRowT==0");
	outRowB<<"t\tvac\tinf"<<endl;*/

	for(beta=1;beta<10.1;beta*=100)
		{
			for(c=0;c<=1;c+=0.02)
				{	
					rowB=0;rowB1=0;rowB2=0;rowB3=0;rowB4=0;
					for(tn=0;tn<TN;tn++)
						{	
							network=new Lattice(M,N);
							for(ts=0;ts<TS;ts++)
								{
									network->initStates(); 
									for(te=0;te<TE0;te++)
										{	network->EpiUpdate(outRowB);
											network->StraUpdate();
											if(network->getNumberOfSs()<=0) {
												break;
												}
												
											
										}

											for(te=TE0;te<TE1;te++)
												{	if(network->getNumberOfSs()<=0) {
													rowB+=network->getNumberOfCs()*(TE1-te);
													rowB1+=network->getNumberOfV0()*(TE1-te);
													rowB2+=network->getNumberOfV1()*(TE1-te);
													rowB3+=network->getNumberOfV2()*(TE1-te);
													rowB4+=network->getNumberOfV3()*(TE1-te);

													//if(network->getNumberOfV0()+network->getNumberOfV1()==400)
														//cout<<777777<<endl;

													
													//cout<<1111111111<<endl;
													//cout<<network->getNumberOfCs()<<endl;
													
													//rowB1+=network->getNumberOfSs()*(TE1-TE0);
													//rowB2+=Iseeds*(TE1-TE0);
													//rowB3+=Iseeds*(TE1-TE0);
													//rowB4+=network->getNumberOfFs()*(TE1-TE0);

														break;
													}
													network->EpiUpdate(outRowB);
													
													//rowB1+=network->getNumberOfSs();
													//rowB2+=network->getNumberOfIs();
													//rowB3+=network->getNumberOfRs();
													//rowB4+=network->getNumberOfFs();
													network->StraUpdate();

													//cout<<network->getNumberOfV2()<<endl;

													rowB+=network->getNumberOfCs();
													rowB1+=network->getNumberOfV0();
													rowB2+=network->getNumberOfV1();
													rowB3+=network->getNumberOfV2();
													rowB4+=network->getNumberOfV3();
													
												}

										}
										
								delete network;	



				}

					//cout<<rowB4<<endl;
			
	rowB = rowB/TS/TN/N/(TE1-TE0)/M;
	rowB1 = rowB1/TS/TN/N/(TE1-TE0)/M;
	rowB2 = rowB2/TS/TN/N/(TE1-TE0)/M;
	rowB3 = rowB3/TS/TN/N/(TE1-TE0)/M;
	rowB4 = rowB4/TS/TN/N/(TE1-TE0)/M;
	//cout<<rowB4<<endl;
	outRowB<<beta<<"\t"<<c<<"\t"<<rowB<<"\t"<<rowB1<<"\t"<<rowB2<<"\t"<<rowB3<<"\t"<<rowB4<<endl;
	if (1-rowB<0.05)
		{
		NVac = (int) ceil(0.95*M*N);}
	else if (rowB<0.05){
	NVac = (int) ceil(0.05*M*N);}
	else {NVac = (int) ceil(rowB*M*N);}
                    
                    
                    double tempB = rowB1+rowB3;
                    if (1-tempB<0.05)
                    {
                        NVacB = (int) ceil(0.95*M*N);}
                    else if (tempB<0.05){
                        NVacB = (int) ceil(0.05*M*N);}
                    else {NVacB = (int) ceil(tempB*M*N);}
                    
                    
						}
            
            
									for(c=1;c>=-.01;c-=0.02)
				{	
					rowB=0;rowB1=0;rowB2=0;rowB3=0;rowB4=0;
					for(tn=0;tn<TN;tn++)
						{	
							network=new Lattice(M,N);
							for(ts=0;ts<TS;ts++)
								{
									network->initStates(); 
									for(te=0;te<TE0;te++)
										{	network->EpiUpdate(outRowB);
											network->StraUpdate();
											if(network->getNumberOfSs()<=0) {
												break;
												}
												
											
										}

											for(te=TE0;te<TE1;te++)
												{	if(network->getNumberOfSs()<=0) {
													rowB+=network->getNumberOfCs()*(TE1-te);
													rowB1+=network->getNumberOfV0()*(TE1-te);
													rowB2+=network->getNumberOfV1()*(TE1-te);
													rowB3+=network->getNumberOfV2()*(TE1-te);
													rowB4+=network->getNumberOfV3()*(TE1-te);

														break;
													}
													network->EpiUpdate(outRowB);
													//rowB+=network->getNumberOfCs();
													//rowB1+=network->getNumberOfSs();
													//rowB2+=network->getNumberOfIs();
													//rowB3+=network->getNumberOfRs();
													//rowB4+=network->getNumberOfFs();
													network->StraUpdate();

													rowB+=network->getNumberOfCs();
													rowB1+=network->getNumberOfV0();
													rowB2+=network->getNumberOfV1();
													rowB3+=network->getNumberOfV2();
													rowB4+=network->getNumberOfV3();
													
												}

										}
										
								delete network;	

								}
			
	rowB = rowB/TS/TN/N/(TE1-TE0)/M;
	rowB1 = rowB1/TS/TN/N/(TE1-TE0)/M;
	rowB2 = rowB2/TS/TN/N/(TE1-TE0)/M;
	rowB3 = rowB3/TS/TN/N/(TE1-TE0)/M;
	rowB4 = rowB4/TS/TN/N/(TE1-TE0)/M;
	outRowB<<beta<<"\t"<<c<<"\t"<<rowB<<"\t"<<rowB1<<"\t"<<rowB2<<"\t"<<rowB3<<"\t"<<rowB4<<endl;
                    
	if (1-rowB<0.05)
		{
		NVac = (int) ceil(0.95*M*N);}
	else if (rowB<0.05){
	NVac = (int) ceil(0.05*M*N);}
	else {NVac = (int) ceil(rowB*M*N);}
                    
                    double tempB = rowB1+rowB3;
                    if (1-tempB<0.05)
                    {
                        NVacB = (int) ceil(0.95*M*N);}
                    else if (tempB<0.05){
                        NVacB = (int) ceil(0.05*M*N);}
                    else {NVacB = (int) ceil(tempB*M*N);}
                    
	
	
							

						}

				}


	outRowB.close();
	return 0;
}
