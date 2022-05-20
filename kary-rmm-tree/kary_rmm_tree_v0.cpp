#include <iostream>
#include <assert.h>
#include <vector>
#include <math.h>
#include <sdsl/int_vector.hpp>
#include <sdsl/bp_support_sada.hpp>
#include "kary_rmm_tree_v0.h"

using namespace std;
using namespace sdsl;

long long int temp;
long long int bNotDivisible_k=0;

RMMTree_Kary::RMMTree_Kary(int_vector<1> &bv,  int sizeBlock,  int w, int order){
	assert(!(order&(order-1)));
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
    this->order = order;
	this->numberLeaves = ceil( (double)bv.size()/(sizeBlock*order) );
	this->numberNodes = ceil((double)(order*numberLeaves-1)/(order-1));
	this->height = cLog_m(this->numberLeaves,this->order);
	this->tree.resize(this->numberNodes);
	temp = ceil(sqrt(order));
}

uint16_t RMMTree_Kary::reverse_16(uint16_t x){
	uint16_t y;
	unsigned char *q = (unsigned char *)&y;
	unsigned char *p = (unsigned char *)&x;
	q[1] = BitReverseTable256[p[0]];
	q[0] = BitReverseTable256[p[1]];
	return y;
}

long long int RMMTree_Kary::bitsread(uint64_t idx){
	uint64_t word;
	if(int(idx+16)>=size){
		return bNotDivisible_k;
	}
 	word = bv.data()[idx>>6];
	auto x = reverse_16( (word >> (idx & 0x3f)) & bits::lo_set[16]);
	return x;
}

unsigned long long RMMTree_Kary::fLog_2(unsigned long long  n){
	return  (8*sizeof (unsigned long long) - __builtin_clzll(n) - 1);
}

unsigned long long RMMTree_Kary::cLog_2(unsigned long long  n){
	return  fLog_2(2*n -1);
}

long long int RMMTree_Kary::getNumberLeaves(){
	return numberLeaves;
}

unsigned long long RMMTree_Kary::cLog_m(unsigned long long  n,unsigned long long  m){
	return ceil((double)log2(n)/log2(m));
}

unsigned long long RMMTree_Kary::fLog_m(unsigned long long  n,unsigned long long  m){
	return log2(n)/log2(m);
}

long long int RMMTree_Kary::min(long long int a , long long int b){
	return (a < b )? a:b;
}

long long int RMMTree_Kary::leafInTree(long long int k){
	if(height==0)return 0;
	
	long long int nNodesPrevLevel = 1 << (temp*(height-1));
	long long int totalNodesExcludingLastLevel = ceil((double)(order*nNodesPrevLevel -1)/(order-1));
	if(k < numberNodes -totalNodesExcludingLastLevel){
		return totalNodesExcludingLastLevel+k; 
	}
	return totalNodesExcludingLastLevel  - numberLeaves + k;
}

long long int RMMTree_Kary::numLeaf(long long int v){
	if(height == 0) return 0;
	
	long long int nNodesPrevLevel = 1 << (temp*(height-1));
	long long int totalNodesExcludingLastLevel = (order*nNodesPrevLevel -1)/(order-1);

	if(v >= totalNodesExcludingLastLevel){
		return v - totalNodesExcludingLastLevel;
	}
	return v - totalNodesExcludingLastLevel + numberLeaves;
}


int RMMTree_Kary::numKey(int k,long long int i){
	long int startBlock= k*sizeBlock*order;
	return (i-startBlock)/sizeBlock;
}

void RMMTree_Kary::buildingTree(){
    buildingTableC();
	buildingLeaves();
	buildingInternalNodesRoot();
}

void RMMTree_Kary::buildingTableC(){
	Key key;
	key.excess = 0;
	key.excessMax = 0 - w;
	key.excessMin = w;
	key.numberExcessMin = 0;

	for(long long int i=0; i < (1<<w); i++){
		tableC.push_back(key);

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
		bNotDivisible_k = tableC.size();
		tableC.push_back(key);
		
		tableC[bNotDivisible_k].excess += (bv[i])? 1 : -1;
		tableC[bNotDivisible_k].excessMax= tableC[bNotDivisible_k].excess;
		tableC[bNotDivisible_k].excessMin = tableC[bNotDivisible_k].excess;
		tableC[bNotDivisible_k].numberExcessMin = 1;
		
		for(i=i+1; i < size; i++){
			tableC[bNotDivisible_k].excess += (bv[i] ==1)? 1 : -1;
			if(tableC[bNotDivisible_k].excessMax < tableC[bNotDivisible_k].excess)tableC[bNotDivisible_k].excessMax = tableC[bNotDivisible_k].excess;
			if(tableC[bNotDivisible_k].excessMin > tableC[bNotDivisible_k].excess){
				tableC[bNotDivisible_k].excessMin = tableC[bNotDivisible_k].excess;
				tableC[bNotDivisible_k].numberExcessMin =1;
			}
			else if(tableC[bNotDivisible_k].excessMin == tableC[bNotDivisible_k].excess){
				tableC[bNotDivisible_k].numberExcessMin++;
			}

		}
	}
}

void RMMTree_Kary::buildingLeaves(){
	long long int x,v;
	int numKey=0;
	Key key;
	
	for(long long int k=0; k < numberLeaves;k++){
		v = leafInTree(k);
		while (tree[v].nKeys < order && numKey < ceil((double)size/sizeBlock)){
			key.excess = 0;
			key.excessMax = 0 -w;
			key.excessMin = 0 +w;
			key.numberExcessMin = 0;
			
			for(long long int p = (numKey*(sizeBlock/w))+1;p<=((numKey+1)*sizeBlock)/w && ((p-1)*w)<size;p++){
				x = bitsread((p-1)*w);
				if(key.excess + tableC[x].excessMax > key.excessMax){
					key.excessMax = key.excess + tableC[x].excessMax;
				}
				if(key.excess + tableC[x].excessMin < key.excessMin){
					key.excessMin = key.excess + tableC[x].excessMin;
					key.numberExcessMin = tableC[x].numberExcessMin;
				}else if(key.excess + tableC[x].excessMin == key.excessMin){
					key.numberExcessMin = key.numberExcessMin + tableC[x].numberExcessMin;
				}
				key.excess += tableC[x].excess;
			}
			tree[v].keys.push_back(key);
			tree[v].nKeys++;
			numKey++;
		}
	}
}

void RMMTree_Kary::buildingInternalNodesRoot(){
	Key key;
	for(long long int v= numberNodes - numberLeaves -1; v>=0;v--){
		for(long long int child=1; child <= order && (v*order+child)<numberNodes;child++){
			key.excess = 0;
			key.excessMax = 0 -w;
			key.excessMin = 0+w;
			key.numberExcessMin = 0;
			
			long long int targetSon = (order*v)+child;
			
			for(long long int i=0; i < tree[targetSon].nKeys;i++){
				if(key.excessMax < key.excess + tree[targetSon].keys[i].excessMax){
					key.excessMax = key.excess + tree[targetSon].keys[i].excessMax;
				}
				if(key.excessMin > key.excess + tree[targetSon].keys[i].excessMin){
					key.excessMin = key.excess + tree[targetSon].keys[i].excessMin;
					key.numberExcessMin = tree[targetSon].keys[i].numberExcessMin;
				}else if(key.excessMin == key.excess + tree[targetSon].keys[i].excessMin){
					key.numberExcessMin+= tree[targetSon].keys[i].numberExcessMin;
				}
				key.excess += tree[targetSon].keys[i].excess;
			}
			tree[v].keys.push_back(key);
			tree[v].nKeys++;
		}
	}
}

long long int RMMTree_Kary::fwdKey(long long int i,long long int v,int key,long long int k,int d,int &dr){
	long long int j;
	for(;key < tree[v].nKeys;key++){
		if((dr + tree[v].keys[key].excessMin <= d) && (d<= dr +tree[v].keys[key].excessMax)){
			j= fwdBlock(i,d,dr);
			if(dr == d) return j;
		}
		dr +=  tree[v].keys[key].excess;
		i = (order*k+key+1)*sizeBlock -1;

		if(dr==d)return i;
	}
	return size;
}

long long int RMMTree_Kary::fwdBlock(long long int i,int d,int &dr){
	long long int p;
	long long int fb = ceil((double)(i+1)/w);
	long long int lb = ceil((double)(i+2)/sizeBlock)* (sizeBlock/w);
	
	for(long long int j=i+1;j<=(fb*w)-1 && j<size;j++){
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
	
	for(long long int j= (p-1)*w; j <= (p*w)-1 && j<size;j++){
		dr += (bv[j] ==1)? 1:-1;
		if(dr == d)return j;
	}
    return size;
}

long long int RMMTree_Kary::fwdVerifySibling(long long int &v, int &dr, int d){
	long long int parent = (v-1)/order;
	long long int child = (parent*order+tree[parent].nKeys) -v;
	
	if(child==0)return size;
	v++;
	
	for( ;child >0 && v<numberNodes;child--){
		for(long long int key=0;key<tree[v].nKeys;key++){
			if((dr + tree[v].keys[key].excessMin <= d) && (d<= dr +tree[v].keys[key].excessMax)){
				if(v<numberNodes - numberLeaves){
					v = (v*order)+1+key;
					key =0;
				}
				return key;
			}
			dr +=  tree[v].keys[key].excess;
		}
		if(child-1 > 0)v++;
	}
	
	return size;
}

long long int RMMTree_Kary::fwdSearch(long long int i, int d){
	assert((i+1)>=0 && (i+1)< size);
	
	int dr=0;
	long long int k = (i+1)/(sizeBlock*order);
	long long int v = leafInTree(k);
	int key = numKey(k,i+1); 
	long long int j = fwdBlock(i,d,dr);
	if(dr == d)return j;
	
	key++;
	if(key < tree[v].nKeys){
		j= fwdKey((order*k+key)*sizeBlock -1,v,key, k,d,dr);
		if(dr == d)return j;
	}

	while( v!=0 && (key=fwdVerifySibling(v,dr,d))==size ){
		v = floor((v-1)/order);	
	} 
	
	if(v==0 && key==size)return size;

	while(v < numberNodes - numberLeaves){
		
		for(key=0;key<tree[v].nKeys;key++){
			if((dr+tree[v].keys[key].excessMin <=d)&& (dr + tree[v].keys[key].excessMax >=d) ){
				v = (v*order)+1+key;
				key=0;
				break;
			}
			else{
				dr += tree[v].keys[key].excess;
			}
		}
	}
	
	k = numLeaf(v);
	j = fwdKey((order*k+key)*sizeBlock -1,v,key, k,d,dr);
	return (dr == d)? j : size;
}

long long int RMMTree_Kary::bwdKey(long long int i,long long int v,int key,long long int k,int d, int &dr){
	long long int j;
	for(; key>=0;key--){
		if((dr - tree[v].keys[key].excess + tree[v].keys[key].excessMin <= d) && (d <= dr - tree[v].keys[key].excess + tree[v].keys[key].excessMax)){
			j = bwdBlock(i,d,dr);
			if(dr==d)return j;
		}
		dr-=tree[v].keys[key].excess;
		i =  (order*k+key)*sizeBlock-1;
		if(dr==d)return i;
	}
	return -1;
}

long long int RMMTree_Kary::bwdBlock(long long int i,int d, int &dr){
	long long int p,x,j;
	long long int fb = i/w;
	long long int lb = (i/sizeBlock) * (sizeBlock/w);
	for(j=min(i,size);j>=fb*w;j--){
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

long long int RMMTree_Kary::bwdVerifySibling(long long int &v, int &dr, int d){
	
	long long int parent = (v-1)/order;
	long long int child = v - (parent*order+1); 
	if(child==0)return -1;
	v--;
	for(; child >0 && v >0;child--){
		for(long long int key =tree[v].nKeys-1; key>=0;key--){
			if((dr - tree[v].keys[key].excess + tree[v].keys[key].excessMin <= d) && (d <= dr - tree[v].keys[key].excess + tree[v].keys[key].excessMax)){
				if(v<numberNodes - numberLeaves){
					v = (v*order)+1+key;
					key = tree[v].nKeys -1;
				}
				return key;
			}

			dr-=tree[v].keys[key].excess;
		}
		if(child-1 >0)v--;
	}
	return -1;
}

long long int RMMTree_Kary::bwdSearch(long long int i,int d){
	assert(i>=0 && i < size);
	if(i==0)return -1;

	int dr=0;
	long long int k = i/(sizeBlock*order);
	long long int v = leafInTree(k);
	int key = numKey(k,i);
	
	long long int j = bwdBlock(i,d,dr);
	if(dr==d) return j;
	
	key--;
	if(key>=0){
		j=bwdKey((order*k+key+1)*sizeBlock-1,v,key,k,d,dr);
		if(dr==d) return j;
	}
	
	while(v!=0 && (key=bwdVerifySibling(v,dr,d))==-1){
		v = (v-1)/order;
	}
	
	if(v==0 && key==-1)return -1;
	
	while(v < numberNodes - numberLeaves){
		for(key=tree[v].nKeys-1;key>=0;key--){
			if( (dr - tree[v].keys[key].excess +tree[v].keys[key].excessMin <= d)&&(dr - tree[v].keys[key].excess +tree[v].keys[key].excessMax >= d)){
				v = (v*order)+1+key;
				key = tree[v].nKeys -1;
				break;
			}
			else{
				dr -= tree[v].keys[key].excess;
			}
		}
	}
	
	k = numLeaf(v);
	if(dr == d)return (k*order*sizeBlock) + ((key+1)*sizeBlock) -1;
	j = bwdKey(k*(order*sizeBlock) + ((key+1)*sizeBlock) -1,v,key,k,d,dr);
	return (dr==d)? j : -1;
}

long long int RMMTree_Kary::minExcess(long long int i, long long int j){
	return 0;
}

long long int RMMTree_Kary::findClose(long long int i){
	assert(i>=0 && i<size);
	if((i == 0) && bv[i]==1)return size -1;
	return (bv[i] == 0) ? i : fwdSearch(i,-1);
}

long long int RMMTree_Kary::findOpen(long long int i){
	assert(i>=0 && i < size);

	if((i == (int)(size-1)) && bv[i]==0)return 0;
	return (bv[i] == 1) ?  i : bwdSearch(i,0)+1;
}

long long int RMMTree_Kary::rmq(long long int i,long long int j){
	return 0;
}


long long int RMMTree_Kary::enclose(long long int i){
	if(i==0)return size;

	if(bv[i]==0)return findOpen(i);

	return bwdSearch(i,-2) + 1;
}

bool RMMTree_Kary::isLeaf(long long int x){
	assert(x>=0 && x < size);

	return  (bv[x] == 1 && bv[x+1]==0);
}

bool RMMTree_Kary::isAncestor(long long int x, long long int y){
	assert(x >=0 && y < size );

	return (x <= y && y< findClose(x));
}

long long int RMMTree_Kary::depth(long long int x){
	assert(x>=0 && x<size);

	return (2*b_rank1(x))-x;
}

long long int RMMTree_Kary::parent(long long int x){
	return enclose(x);
}

long long int RMMTree_Kary::nextSibling(long long int x){
	long long int i = findClose(x)+1;
	return (i<size && bv[i]==1)? i : size;
}

long long int RMMTree_Kary::prevSibling(long long int x){
	return (x!=0 && bv[x-1]==0)? findOpen(x-1) : size;
}

long long int RMMTree_Kary::lastChild(long long int x){
	return (!isLeaf(x))? findOpen(findClose(x)-1) : size;
}

long long int RMMTree_Kary::firstChild(long long int x){
	return (!isLeaf(x) )? x+1 : size;
}

long long int RMMTree_Kary::subtreeSize(long long int x){
	assert(x >=0 && x<size);

	return (bv[x]==1)? ceil((findClose(x)-x+1)/2) : size;
}

long long int RMMTree_Kary::levelAncestor(long long int x,int d){
	assert(x >=0 && x<size);

	return (bv[x]==1)? bwdSearch(x,-d-1)+1 : size;
}

long long int RMMTree_Kary::levelNext(long long int x){
	long long int close = findClose(x);
	if(close <0 || close ==size)return size;
	return fwdSearch(close,1);
}

long long int RMMTree_Kary::levelPrev(long long int x){
	long long int i = bwdSearch(x,0);
	if(i <0 || i==size-1)return -1;
	return findOpen(i+1);
}

long long int RMMTree_Kary::levelLeftMost(int d){
	return (d==1) ? 0 : fwdSearch(0,d-1);
}

long long int RMMTree_Kary::levelRightMost(int d){
	long long int i = bwdSearch(size-1,d);
	if(i + 1 < 0 || i+1 >=size) return -1;
	return (d==1) ? 0 : findOpen(i+1);
}

long long int  RMMTree_Kary::leafRank(long long int x){
	assert(x>=0 && x <size);

	return b_rank10(x);
}

long long int  RMMTree_Kary::leafSelect(long long int t){
	return b_sel10(t)-1;
}

long long int RMMTree_Kary::leftMostLeaf(long long int x){
	return (!isLeaf(x))? leafSelect(leafRank(x-1)+1) : x;
}

long long int RMMTree_Kary::rightMostLeaf(long long int x){
	long long int i = findClose(x);
	if(i<0 || i == size)return size;
	return (!isLeaf(x))? leafSelect(leafRank(i)) : x;
}

long long int RMMTree_Kary::preRank(long long int x){
	assert(x>=0 && x<size);
	return b_rank1(x);
}

long long int RMMTree_Kary::postRank(long long int x){
	long long int i = findClose(x);
	if(i<0 || i == size)return size;
	return b_rank0(i);
}

long long int RMMTree_Kary::preSelect(long long int t){
	if(t==0 || t>size)return size;
	return b_sel1(t);
}

long long int RMMTree_Kary::postSelect(long long int t){
	if(t==0 || t>size)return size;
	return findOpen(b_sel0(t));
}

void RMMTree_Kary::printNode(vector<Key> vector, long long int i){
	cout << "v[" << i<< "].e = " << vector[i].excess << "; "
		 << "v[" << i<< "].m = " << vector[i].excessMin << "; "
		 << "v[" << i<< "].M = " << vector[i].excessMax << "; "
	     << "v[" << i << "].n = " << vector[i].numberExcessMin << ".\n"
		 << endl;
}

void RMMTree_Kary::printTableC(){
	cout << "---- Tabela C, para acelerar a construção da RMM-tree ----" << "\n\n";
	for(long long int i=0;i<(1 << w);i++)printNode(tableC,i);
}

void RMMTree_Kary::printTree(){
	long long int v,leaf;
	cout << " ----- Root and Internal nodes ----- \n";
	for(v=0;v<numberNodes-numberLeaves;v++){
		cout << " Nó " << v << "\n";
		for(long long int k=0;k<tree[v].nKeys; k++){
			cout << "Chave " << k << "\n";
			printNode(tree[v].keys , k);
		}
	}
	
	cout << " ----- Folhas -----" << endl;
	for(long long int k=0;v<numberNodes;v++,k++){
		leaf=numLeaf(v);
		cout << leaf << "-th folha " << " - nó " << v << ": área de cobertura: B[" << leaf*sizeBlock*order << "," << min((leaf+1)*sizeBlock*tree[v].nKeys,size-1)<< "]\n";
		for(long long int k=0;k<tree[v].nKeys; k++){
			cout << "Chave " << k << "\n";
			printNode(tree[v].keys , k);
		}
	}
}

void RMMTree_Kary::printInfoTree(){
    cout << "Tamanho de bloco: " << sizeBlock << '\n'
		 << "Ordem da árvore " << order << '\n'
		 << "Quantidade de folhas: " << numberLeaves << '\n'
		 << "Quantidade de nós: " << numberNodes << '\n'
		 << "Altura da árvore: " << height << '\n'
		 << endl;
}

