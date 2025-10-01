/** 更新于：2007-11-08 21:20 */

#ifndef _Lattice_h
#define _Lattice_h

#include "Network.h"

/***** lattice 网络 。开始*****/
class Lattice : public Network
{
public:
	/** 构造函数的2个参数
	  m  -- 行数
	  n -- 列数m*n=total number of nodes
	*/
	Lattice(long m, long n) : Network(m*n), M(m), N(n) {
		initNeighbors();
	}
	virtual void initNeighbors(void);
private:
	long M;
	long N;
};

void Lattice::initNeighbors(void)
{
	long i, j,node;
	
	/**	分配节点的邻居*/
	for(node=0;node<NumberOfNodes;node++) {
		NumberOfNeighbors[node]=4;
		Neighbors[node] = new long [4];
		if(Neighbors[node]==0) handleError("Neighbors[node]==0");
	}

	for(i=0;i<M;i++){
		for(j=0;j<N;j++){
			Neighbors[i*N+j][0]=((M+i-1)%M)*N+j;
			Neighbors[i*N+j][1]=((i+1)%M)*N+j; 
			Neighbors[i*N+j][2]=i*N+(N+j-1)%N;
			Neighbors[i*N+j][3]=i*N+(j+1)%N;
			}

	}

//for(n0=0;n0<8;n0++)
	//cout<<Neighbors[0][n0]<<endl;
	
	
}
/***** Lattice网络 。结束*****/

#endif
