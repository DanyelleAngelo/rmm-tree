#include <iostream>
#include <assert.h>
#include "binary-rmm-tree.h"
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/bits.hpp>
#include <sdsl/bp_support_sada.hpp>
#include <algorithm>
#include <vector>
#include <math.h>


using namespace std;
using namespace sdsl;
long long int bNotDivisible=0;

RMMTree_Bin::RMMTree_Bin(int_vector<1> &bv, int sizeBlock,int w){
	this->bv = bv;
	this->b_rank1 = rank_support_v<1>(&bv);
	this->b_rank0 = rank_support_v<0>(&bv);
	this->b_rank10 = rank_support_v<10,2>(&bv);
	this->b_sel1 = select_support_mcl<1>(&bv);
	this->b_sel0 = select_support_mcl<0>(&bv);
	this->b_sel10 = select_support_mcl<10,2>(&bv);
	this->sizeBlock = sizeBlock;
	this->w =  w;
	this->size = bv.size();
	this->numberLeaves = ceil((double)this->size/sizeBlock);
	this->numberNodes = (2*this->numberLeaves) -1;
	this->height = cLog_2(this->numberLeaves);
	this->tree.resize(this->numberNodes);
}

uint16_t RMMTree_Bin::reverse_16(uint16_t x){
	uint16_t y;
	unsigned char *q = (unsigned char *)&y;
	unsigned char *p = (unsigned char *)&x;
	q[1] = BitReverseTable256_bin[p[0]];
	q[0] = BitReverseTable256_bin[p[1]];
	return y;
}

long long int RMMTree_Bin::bitsread(uint64_t idx){
	if(int(idx+16) >= size){
		return bNotDivisible;
	}
	uint64_t word = bv.data()[idx>>6];
	auto x = reverse_16( (word >> (idx & 0x3f)) & bits::lo_set[16]);
	return x;
}

unsigned long long RMMTree_Bin::fLog_2(unsigned long long  n){
	return  log2(n);
}

unsigned long long RMMTree_Bin::cLog_2(unsigned long long  n){
	return  ceil(log2(n));
}

long long int RMMTree_Bin::getNumberLeaves(){
	return numberLeaves;
}

long long int RMMTree_Bin::min(long long int a , long long int b){
	return (a < b )? a:b;
}

long long int RMMTree_Bin::leafInTree(long long int k){
	long long int t = 1 << height;
	if(k < (2*numberLeaves) - t)return t-1+k;
	else return t-1-numberLeaves+k;
}

long long int RMMTree_Bin::numLeaf(long long int v){
	long long int t = 1 << height;

	if(v >= t-1 )return v - t  + 1 ;
	else return v - t  + numberLeaves + 1 ;
}

void RMMTree_Bin::buildingTableC(){
	Node_bin node;
	node.excess = 0;
	node.excessMax = 0 - w;
	node.excessMin = w;
	node.numberExcessMin = 0;
	for(long long int i=0; i < (1<<w); i++){
		tableC.push_back(node);

		for(long long int j=0; j<w;j++){
			tableC[i].excess += i & (1 << (w-j-1)) ? 1 : -1;

			if(tableC[i].excess > tableC[i].excessMax){
				tableC[i].excessMax = tableC[i].excess;
			}
			if(tableC[i].excess < tableC[i].excessMin){
				tableC[i].excessMin = tableC[i].excess;
				tableC[i].numberExcessMin = 1;
			}else if(tableC[i].excess == tableC[i].excessMin){
				tableC[i].numberExcessMin+=1;
			}
		}
  	}
	  
	if(size % (sizeBlock) !=0){
		long long int restofDivision = size - (sizeBlock * (size/sizeBlock));
		long long int bitsStoredInC =  (restofDivision/w)*w;
		long long int i = (sizeBlock * (size/sizeBlock)) + bitsStoredInC;
		bNotDivisible = tableC.size();
		tableC.push_back(node);

		tableC[bNotDivisible].excess += (bv[i])? 1 : -1;
		tableC[bNotDivisible].excessMax= tableC[bNotDivisible].excess;
		tableC[bNotDivisible].excessMin = tableC[bNotDivisible].excess;
		tableC[bNotDivisible].numberExcessMin = 1;
		
		for(i=i+1; i < size; i++){
			tableC[bNotDivisible].excess += (bv[i] ==1)? 1 : -1;
			if(tableC[bNotDivisible].excessMax < tableC[bNotDivisible].excess)tableC[bNotDivisible].excessMax = tableC[bNotDivisible].excess;
			if(tableC[bNotDivisible].excessMin > tableC[bNotDivisible].excess){
				tableC[bNotDivisible].excessMin = tableC[bNotDivisible].excess;
				tableC[bNotDivisible].numberExcessMin =1;
			}
			else if(tableC[bNotDivisible].excessMin == tableC[bNotDivisible].excess){
				tableC[bNotDivisible].numberExcessMin++;
			}

		}
	}
}

void RMMTree_Bin::buildingTree(){
    buildingTableC();
	buildingLeaves();
	buildingInternalNodesRoot();
}

void RMMTree_Bin::buildingLeaves(){
	long long int x,v;
	
	for(long long int k=0;k<numberLeaves;k++){
		v = leafInTree(k);
		tree[v].excess = 0;
		tree[v].excessMax = 0 - w;
		tree[v].excessMin =w;
		tree[v].numberExcessMin =0;
		
		for(long long int p = (k*(sizeBlock/w))+1;p<=((k+1)*sizeBlock)/w;p++){
			x = bitsread((p-1)*w);
			if(tree[v].excess + tableC[x].excessMax > tree[v].excessMax){
				tree[v].excessMax = tree[v].excess + tableC[x].excessMax;
			}
			if(tree[v].excess + tableC[x].excessMin < tree[v].excessMin){
				tree[v].excessMin = tree[v].excess + tableC[x].excessMin;
				tree[v].numberExcessMin = tableC[x].numberExcessMin;
			}else if(tree[v].excess + tableC[x].excessMin == tree[v].excessMin){
				tree[v].numberExcessMin = tree[v].numberExcessMin + tableC[x].numberExcessMin;
			}
			tree[v].excess += tableC[x].excess;
		}
	}
}

void RMMTree_Bin::buildingInternalNodesRoot(){
    for(long long int i=numberNodes - numberLeaves-1;i>=0;i--){
		long long int vL = (2*i)+1;
		long long int vR = (2*i)+2;

		tree[i].excess = tree[vL].excess + tree[vR].excess;

		long long int excessTemp = tree[vL].excess + tree[vR].excessMax;
		tree[i].excessMax = (tree[vL].excessMax >= excessTemp) ? tree[vL].excessMax : excessTemp;

		excessTemp = tree[vL].excess + tree[vR].excessMin;
		if(tree[vL].excessMin < excessTemp){
			tree[i].excessMin = tree[vL].excessMin;
			tree[i].numberExcessMin = tree[vL].numberExcessMin;
		}else if(tree[vL].excessMin > excessTemp){
			tree[i].excessMin = excessTemp;
			tree[i].numberExcessMin = tree[vR].numberExcessMin;
		}else{
			tree[i].excessMin = tree[vL].excessMin;
			tree[i].numberExcessMin = tree[vL].numberExcessMin + tree[vR].numberExcessMin;
		}
	}
}

void RMMTree_Bin::printNode(vector<Node_bin> vector, long long int i){
	cout << "v[" << i<< "].e = " << vector[i].excess << "; "
		 << "v[" << i<< "].m = " << vector[i].excessMin << "; "
		 << "v[" << i<< "].M = " << vector[i].excessMax << "; "
	     << "v[" << i << "].n = " << vector[i].numberExcessMin << ".\n"
		 << endl;
}

void RMMTree_Bin::printTree(){
	long long int v,k;
	cout << " ----- Root ----- \n";
	printNode(tree, 0);
	cout << " ----- Internal nodes ----- \n";
	for(v=1;v<numberNodes-numberLeaves;v++){
		cout << " Nó " << v << "\n";
		printNode(tree, v);
	}

	cout << " ----- Folhas -----" << endl;
	for(;v<numberNodes;v++){
		k=numLeaf(v);
		cout <<  k << "-th folha " << " - nó " << v << ": área de cobertura: B[" << k*sizeBlock << "," << (k+1)*sizeBlock -1<< "]\n";
		printNode(tree,v);
	}
}

void RMMTree_Bin::printTableC(){
	cout << "---- Tabela C, para acelerar a construção da RMM-tree ----" << "\n\n";
	for(long long int i=0;i<(1 << w);i++)printNode(tableC,i);
}

void RMMTree_Bin::printInfoTree(){
    cout << "Tamanho de bloco: " << sizeBlock << '\n'
		 << "Quantidade de folhas: " << numberLeaves << '\n'
		 << "Quantidade de nós: " << numberNodes << '\n'
		 << "Altura da árvore: " << height << '\n'
		 << endl;
}

long long int RMMTree_Bin::fwdBlock(long long int i,int d,int &dr){
	long long int p;
	long long int fb = ceil((double)(i+1)/w);
	long long int lb = ceil((double)(i+2)/sizeBlock)* (sizeBlock/w); 
	
	for(long long int j=i+1;j<=(fb*w)-1  && j<size;j++){
		dr += (bv[j] == 1)? 1 : -1;
		if(dr == d)return j;
	}
	
	for(p=fb+1; p<=lb;p++){
		long long int x = bitsread((p-1)*w);

		if(dr + tableC[x].excessMin <= d && dr + tableC[x].excessMax >= d){
			break;
		}
		dr += tableC[x].excess;
	}
	
	if(p > lb)return lb*sizeBlock;
	
	for(long long int j= (p-1)*w; j <= (p*w)-1  && j<size;j++){
		dr += (bv[j] ==1)? 1:-1;
		if(dr == d)return j;
	}
    return size;
}

long long int RMMTree_Bin::bwdBlock(long long int i,int d,int &dr){
	long long int p,x,j;
	long long int fb = i/w;
	long long int lb = (i/sizeBlock) * (sizeBlock/w);
	
	for(j=i;j>=fb*w;j--){
		dr += (bv[j] == 1)? -1 : 1;
		if(dr == d)return j-1;
	}
	
	for(p=fb-1;p>=lb;p--){
		x = bitsread(p*w);
		if( (dr - tableC[x].excess + tableC[x].excessMin <= d)&& (dr - tableC[x].excess + tableC[x].excessMax >= d) ){
			break;
		}
		dr -= tableC[x].excess;
	}

	if(p < lb)return (lb*w)-1;

	for(j=(p+1)*w-1; j>=p*w; j--){
		dr += (bv[j] == 1)? -1 : 1;
		if(dr == d)return j-1;
	}
	
	return size;
}

long long int RMMTree_Bin::minBlock(long long int i,long long int j, int &d){
	long long int p,x, m=w;
	long long int fb = ceil((double)i/w)+1;
	long long int lb = j/w;
	long long int lim = min(j,(fb*w)-1);

	for(p=i; p <= lim;p++){
		d += (bv[p] == 1)? 1 : -1;
		if(d < m)m=d;
	}
	
	if(j<=lim)return m;
	
	for(p=fb+1; p <=lb;p++){
		x = bitsread(p*w);

		if((d + tableC[x].excessMin) < m)m = d + tableC[x].excessMin;
		d += tableC[x].excess;
	}
	
	for(p=lb*w;p<=j;p++){
		d += (bv[p] == 1)? 1 : -1;
		if(d < m)m=d;
	}
	
	return m;
}

long long int RMMTree_Bin::maxBlock(long long int i,long long int j,int &d){
	long long int p,x, eM=-w;
	long long int fb = ceil((double)(i+1)/w);
	long long int lb = j/w;
	long long int lim = min(j,(fb*w)-1);
	
	for(p=i; p <= lim;p++){
		d += (bv[p] == 1)? 1 : -1;
		if(d > eM)eM=d;
	}
	
	if(j<=lim)return eM;
	
	for(p=fb+1; p <=lb;p++){
		x = bitsread(p*w);

		if((d + tableC[x].excessMax) > eM)eM = d + tableC[x].excessMax;
		d += tableC[x].excess;
	}
	
	for(p=lb*w;p<=j;p++){
		d += (bv[p] == 1)? 1 : -1;
		if(d > eM)eM=d;
	}

	return eM;
}

long long int RMMTree_Bin::minCountBlock(long long int i,long long int j,long long int m,int &d){
	long long int fb = ceil((double)(i+1)/w);
	long long int lb = j/w;
	long long int lim = min(j,(fb*w)-1);
	long long int n=0,x;

	for(long long int p=i; p<=lim; p++){
		d += (bv[p] ==1)? 1 : -1;
		if(d==m)n++;
	}

	if(j<=lim) return n;

	for(long long int p=fb+1; p<=lb; p++){
		x = bitsread(p*w);
		if((d + tableC[x].excessMin) == m)n++;
		d += tableC[x].excess;
	}

	for(long long int p=lb*w;p<=j;p++){
		d += (bv[p] == 1)? 1 : -1;
		if(d == m)n++;
	}
	return n;
}

long long int RMMTree_Bin::minSelectBlock(long long int i,long long int j,long long int m, long long int &t, int &d){	
	long long int fb = ceil((double)(i+1)/w);
	long long int lb = j/w;
	long long int lim = min(j,(fb*w)-1);
	long long int p,x;
	
	for(p=i; p<=lim; p++){
		d += (bv[p] ==1)? 1 : -1;
		if(d==m){
			t-=1;
			if(t==0)return p;
		}
	}
	
	if(j<=lim) return p;
	
	for(p=fb; p<=lb; p++){
		x = bitsread(p*w);
		
		if((d + tableC[x].excessMin) <=  m)break;
		
		d += tableC[x].excess;
	}

	if(p>lb)return lb*w;
	
	for(p=lb*w;p<=j;p++){
		d += (bv[p] == 1)? 1 : -1;
		if(d == m){
			t-=1;
			if(t==0)return p;
		}
	}
	
	return p;
}

long long int RMMTree_Bin::fwdSearch(long long int i,int d){
	assert((i+1)>=0 && (i+1)< size);

	long long int j,k,v;
	int dr=0;

	j= fwdBlock(i,d,dr);
    if(dr == d) return j;

	k = (i+1)/sizeBlock;
	v = leafInTree(k);
	
	while( ((v+1)&(v+2))!=0 && !( (dr+tree[v+1].excessMin <= d) && (d<=dr+tree[v+1].excessMax) ) ){
		if((v&1)==1)dr += tree[v+1].excess;
		v =(v-1)/2;
	}
	
	if(((v+1)&(v+2)) ==0)return size;
	
	v++;
	
	while(v < numberLeaves-1){
		if((dr + tree[(2*v)+1].excessMin <= d) && (dr + tree[(2*v)+1].excessMax >= d)){
			v = (2*v)+1;
		}
		else{
			dr += tree[(2*v)+1].excess;
			v = (2*v)+2;
		}
	}

	k = numLeaf(v);
	j = fwdBlock((k*sizeBlock)-1,d,dr);
	
	return (dr == d)? j : size;
}

long long int RMMTree_Bin::bwdSearch(long long int i,int d){
	assert(i>=0 && i < size);
	long long int j,k,v;
	int dr=0;

	if(i==0)return -1;

	j= bwdBlock(i,d,dr);
    if(dr == d) return j;

	k = i/sizeBlock;
	v = leafInTree(k);
	
	
	while(((v+1)&v) !=0 &&  !((dr - tree[v-1].excess + tree[v-1].excessMin <= d) && (d <= dr - tree[v-1].excess + tree[v-1].excessMax)) ){
		
		if((v&1)== 0){
			dr-=tree[v-1].excess;
		}
		v = (v-1)/2;
	}
	
	if( ((v+1)&v) ==0 )return -1;

	v--;
	
	while(v < numberLeaves-1){
		if( (dr - tree[(2*v)+2].excess +tree[(2*v)+2].excessMin <= d)&&(dr - tree[(2*v)+2].excess +tree[(2*v)+2].excessMax >= d)){
			v = (2*v) + 2;
		}else{
			dr -= tree[(2*v)+2].excess;
			v = (2*v) +1;
		}
	}

	k = numLeaf(v);
	if(dr == d)return ((k+1)*sizeBlock)-1;
	j=bwdBlock( ((k+1)*sizeBlock)-1,d,dr);

	return (dr==d)? j : -1;
}

long long int RMMTree_Bin::minExcess(long long int i,long long int j){
	assert(	(i<=j) && (i>=0 && j < size));

	int d=0;
	long long int k_i =i/sizeBlock;
	long long int k_j = (j+1)/sizeBlock;
	long long int m   = minBlock(i,min(((k_i+1)*sizeBlock)-1,j),d);
	
	if(j <= ((k_i+1)*sizeBlock)-1)return m;
	
	long long int v   = leafInTree(k_i);
	long long int v_j = leafInTree(k_j)+1;
	
	while(v+1 > v_j ||  (int)((v_j)/ (1<< (int)(fLog_2(v_j) - (int)fLog_2(v+2)) ) ) !=v+2){
		if((v&1)==1){
			if(d+tree[v+1].excessMin < m)m = d+tree[v+1].excessMin;
			d+= tree[v+1].excess; 
		}
		v = (v-1)/2;
		if(v==-1)break;
		
	}
	
	v++;
	
	while(v < numberLeaves-1){
		if(d+tree[v].excessMin >=m)return m;
		
		if( (int)((v_j)/ (1<<  (int)(fLog_2(v_j) - (int)fLog_2((2*v)+2)) ))  != (2*v)+2){
			if(d+tree[(2*v)+1].excessMin < m)m= d+tree[(2*v)+1].excessMin;
			d+= tree[(2*v)+1].excess;
			v = (2*v)+2;
		}
		else v = (2*v)+1;
	}
	if(d+tree[v].excessMin >=m)return m;
	int dr=0;
	long long int mr = minBlock(k_j*sizeBlock,j,dr);
	return (d + mr < m) ? (d + mr) : m;
}

long long int RMMTree_Bin::maxExcess(long long int i,long long int j){
	assert(	(i<=j) && (i>=0 && j < size));

	int d=0;
	long long int k_i =i/sizeBlock;
	long long int k_j = (j+1)/sizeBlock;
	long long int eM   = maxBlock(i,min(((k_i+1)*sizeBlock)-1,j),d);
	
	if(j <= (k_i+1)*sizeBlock-1)return eM;
	
	long long int v   = leafInTree(k_i);
	long long int v_j = leafInTree(k_j)+1;
	
	while(v+1 > v_j ||   (int)((v_j)/ (1<< (int)(fLog_2(v_j) - fLog_2(v+2)) ) ) !=v+2){
		if((v&1)==1){
			if(d+tree[v+1].excessMax > eM)eM = d+tree[v+1].excessMax;
			d+= tree[v+1].excess; 
		}
		v = (v-1)/2;
	}

	v++;
	
	while(v < numberLeaves-1){

		if(d+tree[v].excessMax <=eM)return eM;
		if( (int)((v_j)/ (1<<  (int)(fLog_2(v_j) - fLog_2((2*v)+2)) ))  != (2*v)+2 ){
			if(d+tree[(2*v)+1].excessMax > eM)eM= d+tree[(2*v)+1].excessMax;
			d+= tree[(2*v)+1].excess;
			v = (2*v)+2;
		}
		else v = (2*v)+1;
	}

	if(d+tree[v].excessMax <=eM)return eM;

	int dr=0;
	long long int mr = maxBlock(k_j*sizeBlock,j,dr);
	return (d + mr > eM) ? (d + mr) : eM;
}

long long int RMMTree_Bin::minCount(long long int i,long long int j){
	assert(	(i<=j) && (i>=0 && j < size));

	long long int m = minExcess(i,j);
	int d = 0;
	long long int k_i = i/sizeBlock;
	long long int k_j = (j+1)/sizeBlock;
	long long int n = minCountBlock(i,min(((k_i+1)*sizeBlock)-1,j),m,d);
	
	if(j<=((k_i+1)*sizeBlock)-1)return n;

	long long int v = leafInTree(k_i);
	long long int v_j = leafInTree(k_j)+1;

	while(v+1 > v_j ||  (int)((v_j)/ (1<< (int)(fLog_2(v_j) - fLog_2(v+2)) ) ) !=v+2){

		if((v&1)==1){
			if(d+tree[v+1].excessMin == m)n += tree[v+1].numberExcessMin;
			d+= tree[v+1].excess; 
		}
		v =(v-1)/2;
	}
	
	v++;
	
	
	while(v < numberLeaves-1){

		if(d+tree[v].excessMin > m)return n;
		if( (int)((v_j)/ (1<<  (int)(fLog_2(v_j) - fLog_2((2*v)+2)) ))  != (2*v)+2 ){
			if(d+tree[(2*v)+1].excessMin == m )n+= tree[(2*v)+1].numberExcessMin;
			d+= tree[(2*v)+1].excess;
			v = (2*v)+2;
		}
		else v = (2*v)+1;
	}
	
	if(d+tree[v].excessMin > m)return n;
	
	return n+minCountBlock(k_j*sizeBlock,j,m,d);
}

long long int RMMTree_Bin::minSelectExcess(long long int i,long long int j, long long int t){
	assert(	(i<=j) && (i>=0 && j < size));

	int m = minExcess(i,j);
	
	int d = 0;
	long long int k_i = i/sizeBlock;
	long long int k_j = j/sizeBlock;

	long long int p = minSelectBlock(i,min(((k_i+1)*sizeBlock)-1,j),m,t,d);
	
	if(j<=(k_i+1)*sizeBlock || t == 0)return p;

	long long int v = leafInTree(k_i);
	long long int v_j = leafInTree(k_j)+1;
	
	while(v+2 > v_j ||  (int)((v_j)/ (1<< (int)(fLog_2(v_j) - fLog_2(v+2)) ) ) !=v+2 ){

		if((v&1)==1){
			if(d+tree[v+1].excessMin == m){
				if(t<=tree[v+1].numberExcessMin)break;
				t -= tree[v+1].numberExcessMin;
			}
			d+= tree[v+1].excess; 
		}
		v =(v-1)/2;
	}

	v++;
	
	
	while(v < numberLeaves-1){
		if( (int)((v_j)/ (1<<  (int)(fLog_2(v_j) - fLog_2((2*v)+2))))== (2*v)+1 ){
			v= (2*v)+1;
		}
		else if(d + tree[(2*v)+1].excessMin > m){
			d += tree[(2*v)+1].excess;
			v = (2*v)+2;
		}
		else if(t<=tree[(2*v)+1].numberExcessMin){
			v = (2*v)+1;
		}
		else{
			t-=tree[(2*v)+1].numberExcessMin;
			d += tree[(2*v)+1].excess;
			v = (2*v)+2;
		}
	}
	
	long long int k = numLeaf(v);
	
	p = minSelectBlock(k*sizeBlock,min(((k+1)*sizeBlock)-1,j),m,t,d);
	return (t==0)? p : size;
}

long long int RMMTree_Bin::rmq(long long int i,long long int j){
	long long int m = minExcess(i,j);
	return fwdSearch(i-1,m);
}

long long int RMMTree_Bin::rMq(long long int i,long long int j){
	long long int m = maxExcess(i,j);
	
	if(i==0){
		m += (bv[i]==1)? -1 : 1;
		return fwdSearch(0,m);
	}
	return fwdSearch(i-1,m);
}

long long int RMMTree_Bin::findClose(long long int i){
	assert(i>=0 && i<size);

	if((i == 0) && bv[i]==1)return size -1;
	return (bv[i] == 0) ? i : fwdSearch(i,-1);
}

long long int RMMTree_Bin::findOpen(long long int i){
	assert(i>=0 && i < size);

	if((i == (int)(size-1)) && bv[i]==0)return 0;
	return (bv[i] == 1) ?  i : bwdSearch(i,0)+1;
}

long long int RMMTree_Bin::enclose(long long int i){
	assert(i>=0 && i < size);

	if(i==0)return size;
	if(bv[i]==0)return findOpen(i);
	return bwdSearch(i,-2) + 1;
}

bool RMMTree_Bin::isLeaf(long long int x){
	assert(x>=0 && x < size);

	return  (bv[x] == 1 && bv[x+1]==0);
}

bool RMMTree_Bin::isAncestor(long long int x, long long int y){
	assert(x >=0 && y < size );

	return (x <= y && y< findClose(x));
}

long long int RMMTree_Bin::depth(long long int x){
	assert(x>=0 && x<size);

	return (2*b_rank1(x))-x;
}

long long int RMMTree_Bin::parent(long long int x){
	return enclose(x);
}

long long int RMMTree_Bin::nextSibling(long long int x){
	long long int i = findClose(x)+1;
	return (i<size && bv[i]==1)? i : size;
}

long long int RMMTree_Bin::prevSibling(long long int x){
	return (x!=0 && bv[x-1]==0)? findOpen(x-1) : size;
}

long long int RMMTree_Bin::child(long long int x,long long int t){
	if(t==1)return firstChild(x);

	long long int j = findClose(x);
	
	long long int p = minSelectExcess(x+1,j-1,t-1)+1; 
	
	return (p >=j) ? size : p;
}

long long int RMMTree_Bin::lastChild(long long int x){
	return (!isLeaf(x))? findOpen(findClose(x)-1) : size;
}

long long int RMMTree_Bin::firstChild(long long int x){
	return (!isLeaf(x) )? x+1 : size;
}

long long int RMMTree_Bin::childRank(long long int x){
	if(x-1 ==0) return 1;
	return minCount(parent(x)+1,x);
}

long long int RMMTree_Bin::subtreeSize(long long int x){
	assert(x >=0 && x<size);

	return (bv[x]==1)? ceil((findClose(x)-x+1)/2) : size;
}

long long int RMMTree_Bin::levelAncestor(long long int x,int d){
	assert(x >=0 && x<size);

	return (bv[x]==1)? bwdSearch(x,-d-1)+1 : size;
}

long long int RMMTree_Bin::lca(long long int x,long long int y){
	if(isAncestor(x,y))return x;
	if(isAncestor(y,x))return y;

	return enclose(rmq(x,y)+1);
}

long long int RMMTree_Bin::levelNext(long long int x){
	long long int close = findClose(x);
	if(close < 0 || close >=size)return size;
	return fwdSearch(close,1);
}

long long int RMMTree_Bin::levelPrev(long long int x){
	long long int i = bwdSearch(x,0);
	if(i < 0 || i>=size)return -1;
	return findOpen(i+1);
}

long long int RMMTree_Bin::levelLeftMost(int d){
	return (d==1) ? 0 : fwdSearch(0,d-1);
}

long long int RMMTree_Bin::levelRightMost(int d){
	long long int i = bwdSearch(size-1,d);
	if(i+1<0 || i+1>=size)return -1;
	return (d==1) ? 0 : findOpen(i+1);
}

long long int RMMTree_Bin::deepestNode(long long int x){
	return rMq(x, findClose(x));
}

long long int  RMMTree_Bin::degree(long long int x){
	return minCount(x+1, findClose(x)-1);
}

long long int  RMMTree_Bin::leafRank(long long int x){
	assert(x>=0 && x <size);

	return b_rank10(x);
}

long long int  RMMTree_Bin::leafSelect(long long int t){
	return b_sel10(t)-1;
}

long long int RMMTree_Bin::leftMostLeaf(long long int x){
	return (!isLeaf(x))? leafSelect(leafRank(x-1)+1) : x;
}

long long int RMMTree_Bin::rightMostLeaf(long long int x){
	long long int i = findClose(x);
	if(i == size || i<0)return size;
	return (!isLeaf(x))? leafSelect(leafRank(i)) : x;
}

long long int RMMTree_Bin::preRank(long long int x){
	assert(x>=0 && x<size);
	return b_rank1(x);
}

long long int RMMTree_Bin::postRank(long long int x){
	long long int i = findClose(x);
	if(i==size || i <0)return size;
	return b_rank0(findClose(x));
}

long long int RMMTree_Bin::preSelect(long long int t){
	if(t==0 || t>size)return size;
	return b_sel1(t);
}

long long int RMMTree_Bin::postSelect(long long int t){
	if(t==0 || t>size)return size;
	return findOpen(b_sel0(t));
}
