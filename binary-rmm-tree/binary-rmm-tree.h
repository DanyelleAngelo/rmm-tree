#ifndef RMMTREE_BIN_H
#define RMMTREE_BIN_H

#include <sdsl/int_vector.hpp>
#include <sdsl/bp_support_sada.hpp>
#include <vector>


using namespace sdsl;
using std::vector;

static const unsigned char BitReverseTable256_bin[] ={
        0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
        0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
        0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
        0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
        0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
        0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
        0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
        0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
        0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
        0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
        0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
        0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
        0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
        0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
        0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
        0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF
};

typedef struct Node_bin{
	long int excess;
	long int excessMax;
	long int excessMin;
	long int numberExcessMin;
}Node_bin;

class RMMTree_Bin{
    public:
		bit_vector bv;				// Vetor de bits que armazena a sequ??ncia de par??nteses balanceados
		long long int size;			// Tamanho da sequ??ncia de par??nteses balanceados
		vector<Node_bin> tree;		// Vetor do tipo Node, usado para armazenar a Range-min-max tree

        RMMTree_Bin(int_vector<1> &bv, int sizeBlock, int w);
		
		/*!
		*	@brief ler os bits de idx at?? idx+8, e contabiliza o inteiro correspondente a esta sequ??ncia bin??ria
		*	@param idx: ponto de partida da leitura
		*	@return inteiro correspondente aos bits lidos.
		*/
		long long int bitsread(uint64_t idx);

		/*!
		*	@param n: n??mero a ser retirado
		*	@return ch??o do logar??tmo na base 2 de n
		*/
		unsigned long long fLog_2(unsigned long long  n);

		/*!
		*	@param n: n??mero a ser retirado
		*	@return teto do logar??tmo na base 2 de n
		*/
		unsigned long long cLog_2(unsigned long long  n);
		long long int getNumberLeaves();
		/*!
		*	@brief verifica se a k-th folha est?? no ??ltimo ou pen??ltimo n??vel da ??rvore e calcula sua posi????o.
		*	A ordem das folhas aqui vai de 0 ?? numberLeaves-1.
		*	@param k = k-th folha
		*	@return ??ndice da k-th folha na rmM-tree
		*/
		long long int leafInTree(long long int k);

		/*!
		*	@brief Dado um ??ndice na rmM-tree,verifica o n??vel em que a folha est??, e calcula a sua ordem.
		*	A ordem das folhas vai de 0 ?? numberLeaves-1
		*	@param v : ??ndice da folha na rmM-tree
		*	@return ordem da folha
		*/
		long long int numLeaf(long long int v);

		/*!
		*	@brief constr??i a estrutura da rmM-tree, chamando primeiro a fun????o para pr??-computar
		*	uma tabela que ir?? acelerar o processo, e chamando depois a fun????o para construir as folhas
		*	da rmM-tree e por ??ltimo seus n??s internos.
		*/
		void buildingTree();

		/*!
		*	@brief Imprime a ra??z, depois n??s internos e por ??ltimo folhas da rmM-tree atrav??s da 
		*	fun????o printNode()
		*/
		void printTree();

		/*!
		*	@brief Imprime os valores da tabela de acelera????o atrav??s da fun????o printNode()
		*/
		void printTableC();

		/*!
		*	@brief Imprime as informa????es b??sicas da rmM-tree, n??mero de n??s, n??mero de folhas, altura
		*	e tamanho de cobertura de bloco.
		*/
		void printInfoTree();

		/*!
		*	@brief Realiza um percusor de subida e depois desciada na rmM-tree,  nos n??s a direita
		* 	da folha (e na pr??pria folha)  que cont??m a posi????o i+1 em busca do excesso d.
		*	@param i: Posi????o a partir da qual quero realizar a busca (ser?? incrementado de 1)
		*	@param d: excesso de "1" buscado (profundidade)
		*	@return primeira posi????o j > i em bv onde ocorre o excesso d.
		*/
		long long int fwdSearch(long long int i,int d);

		/*!
		*	@brief Realiza um percusor de subida e depois desciada na rmM-tree atrav??s dos nos n??s 
		* 	a esquerda da folha (e na pr??pria folha)  que cont??m a posi????o i em busca do excesso d.
		*	@param i: Posi????o a partir da qual quero realizar a busca.
		*	@param d: excesso de "1" buscado (profundidade)
		*	@return primeira posi????o j < i em bv onde ocorre o excesso d.
		*/
		long long int bwdSearch(long long int i,int d);

		/*!
		*	@brief Realiza um percurso de subida e descida na rmM-tree em busca do execesso m??nimo
		*	presente no intervalo i,j. O percurso de subida termina ao encontrar o n?? ancestral da
		*	folha que cobre a posi????o j.
		*	@param i: ??ndice a partir do qual deve-se realizar a busca
		*	@param j: ??ndice onde devemos terminar a busca
		*	@return inteiro representando o excesso m??nimo no intervalo i,j 
		*/
		long long int minExcess(long long int i,long long int j);

		/*!
		*	@brief Realiza um percurso de subida e descida na rmM-tree em busca do execesso m??ximo
		*	presente no intervalo i,j. O percurso de subida termina ao encontrar o n?? ancestral da
		*	folha que cobre a posi????o j.
		*	@param i: ??ndice a partir do qual deve-se realizar a busca
		*	@param j: ??ndice onde devemos terminar a busca
		*	@return inteiro representando o excesso m??ximo no intervalo i,j 
		*/
		long long int maxExcess(long long int i,long long int j);

		/*!
		*	@brief Chama minExcess() para contabilizar o excesso m??nimo no intervalo i,j
		*	e depois realiza um percurso de subida e descida na ??rvore, afim de contabilizar o
		*	n??mero de vezes que o excesso aparece na ??rea
		*	@param i: ??ndice a partir do qual deve-se realizar a busca
		*	@param j: ??ndice onde devemos terminar a busca
		*	@return n??mero de vezes que o excesso m??nimo aparece em bv[i,j]
		*/
		long long int minCount(long long int i,long long int j);

		/*!
		*	@brief Chama minExcess() para contabilizar o excesso m??nimo no intervalo i,j
		*	e depois realiza um percurso de subida e descida na ??rvore, afim de obter a posi????o exata
		*	em que o t-th execesso m??nimo ocorre.
		*	@param i: ??ndice a partir do qual deve-se realizar a busca
		*	@param j: ??ndice onde devemos terminar a busca
		*	@param t: ordem em que o excesso m??nimo deve ocorrer
		*	@return ??ndice, p, com i <= p <= j, onde ocorre o t-th excesso m??nimo 
		*/
		long long int minSelectExcess(long long int i,long long int j, long long int t);

		/*!
		*	@brief Chama minExcess() para contabilizar o excesso m??nimo no intervalo i,j
		*	e depois chama a fun????o fwdSearch(), passando como par??metro o excesso m??nimo computado,
		*	e a posi????o i-1, a fim de encontar a posi????o exata em que esse excesso ocorre.
		*	@param i: ??ndice a partir do qual deve-se realizar a busca
		*	@param j: ??ndice onde devemos terminar a busca
		*	@return ??ndice, p, com i <= p <= j, onde ocorre pela primeira vez o excesso m??nimo
		*/
		long long int rmq(long long int i,long long int j);

		/*!
		*	@brief Chama maxExcess() para contabilizar o excesso m??ximo no intervalo i,j
		*	e depois chama a fun????o fwdSearch(), passando como par??metro o excesso m??ximo computado,
		*	e a posi????o i-1, a fim de encontar a posi????o exata em que esse excesso ocorre.
		*	@param i: ??ndice a partir do qual deve-se realizar a busca
		*	@param j: ??ndice onde devemos terminar a busca
		*	@return ??ndice, p, com i <= p <= j, onde ocorre pela primeira vez o excesso m??ximo
		*/
		long long int rMq(long long int i,long long int j);

		/*!
		*	@brief busca o par??nteses de fechamento, corresponde ao par??nteses  de abertura i.
		*	para isso, busca atrav??s da fun????o fwdSearch() a primeira posi????o j >i onde ocorre
		*	o excesso -1.
		*	@param i: ??ndice de um par??nteses de abertura em bv, que cofica um n??
		*	@return ??ndice j do par??nteses de fechamento corresponde ao par??nteses de abertura i, ou i, caso
		*	i codifique um par??nteses de fechamento.
		*/
		long long int findClose(long long int i);

		/*!
		*	@brief busca o par??nteses de abertura, corresponde ao par??nteses  de fechamento i.
		*	para isso, busca atrav??s da fun????o bwdSearch() a primeira posi????o j < i onde ocorre
		*	o excesso 0.
		*	@param i: ??ndice de um par??nteses de fechamento em bv
		*	@return ??ndice j do par??nteses de abertura corresponde ao par??nteses de fechamento i, ou i, caso
		*	i codifique um par??nteses de abertura.
		*/		
		long long int findOpen(long long int i);

		/*!
		*	@brief busca a posi????o j < i , tal que BV[j,findclose(j)] contenha a posi????o i.
		*	@param i: ??ndice de um par??nteses de abertura em bv, que cofica um n??
		*	@return bv.size() caso i=0, findOpen(i) caso bv[i]=1, ou n?? pai de i (computado atrav??s da fun????o bwdSearch(i,-2)+1)
		*/
		long long int enclose(long long int i);

		/*!
		*	@brief  verifica se bv[x] ?? um bit 1 e se o elemento que o sucede ?? um bit 0 para decidir se este ?? um n?? folha.
		*	@param x: ??ndice do par??ntes que codifica o n?? analisado.
		*	@return: true se x codificar um n?? folha, false caso contr??rio.
		*/
		bool isLeaf(long long int x);

		/*!
		*	@brief  verifica se o n?? x ?? ancestral do n?? y, checando se o primeiro cont??m o segundo.
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x
		*	@param y: ??ndice do vetor de par??ntese balanceados que codifica o n??
		*	@return: true se x ?? ancestral de y, false caso contr??rio.
		*/
		bool isAncestor(long long int x, long long int y);

		/*!
		*	@brief  contabiliza o excesso de 1 no intervalo bv[0,x] para saber a profundidade do n?? x.
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x 
		*	@return: profundidade do n?? x
		*/
		long long int depth(long long int x);

		/*!
		*	@brief  Busca o n?? que cont??m o x, mais ?? esquerda de x.
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x 
		*	@return: ??ndice do n?? pai de x
		*/
		long long int parent(long long int x);

		/*!
		*	@brief  busca o irm??o ?? direita de x
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x 
		*	@return: retorna o ??ndice do irm??o ?? direita de x
		*/
		long long int nextSibling(long long int x);

		/*!
		*	@brief  busca o irm??o ?? esquerda de x
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x 
		*	@return: retorna o ??ndice do irm??o ?? esquerda de x
		*/
		long long int prevSibling(long long int x);

		/*!
		* 	@brief usa minSelect para computar o t-th filho do n?? x (se houver)
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x
		*	@param t: ordem do filho de x a ser buscado
		*	@return : ??ndice do par??nteses de abertura que codifica o t-th filho de x, se houver, e  bv,size() caso contr??rio
		*/
		long long int  child(long long int x,long long int t);

		/*!
		*	@brief  verifica se o elemento anterior a close de i, ?? um n?? (se i n??o ?? folha)
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x 
		*	@return: retorna o ??ndice do i-th filho de x, se houver, e  bv,size() caso contr??rio
		*/
		long long int lastChild(long long int x);

		/*!
		*	@brief verifica se o n?? ?? uma folha, se n??o for retorna o primeiro elemento contido por ele
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x 
		*	@return: Primeiro filho do n?? x, se houver, e  bv,size() caso contr??rio
		*/
		long long int firstChild(long long int x);

		/*! 
		*	@brief contabiliza o n??mero de irm??os a esquerda do n?? x
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x 
		*	@return a quantidade de n??s a esquerda de x
		*/
		long long int childRank(long long int x);

		/*!
		*	@brief contabiliza o n??mero de bits 1 no intervalo b[x,close(x)]
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x
		*	@return tamanho da sub??rvore enra??zada em x
		*/
		long long int subtreeSize(long long int x);

		/*!
		*	@brief  busca o ancestral de x que est?? d n??veis acima dele
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x
		* 	@param d: quantidade de n??veis acima de x onde teremos a resposta
		*	@return: ??ndice j do ancestral de x tal que depth(j)=depth(x)-d
		*/
		long long int levelAncestor(long long int x,int d);

		/*!
		*	@brief  busca o ancestral comum (mais pr??ximo) de x e y
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x
		* 	@param y: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? y
		*	@return: ??ndice do ancestral comum de x e y.
		*/
		long long int lca(long long int x,long long int y);

		/*!
		*	@brief usa fwdSearch para encontrar o pr??ximo n?? a direita de x (n??o necessariamente o irm??o de x,
		*	a ideia ?? pegar o pr??ximo elemento do n??vel)
		* 	com a mesma profundidade de x
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x
		*	@return : ??ndice do primeiro n?? a direita de x, com a mesma profundidade de x; ou size se a resposta n??o for encotrada
		*/
		long long int levelNext(long long int x);

		/*!
		*	@brief usa bwdSearch para encontrar o elemnto a esquerda de x (n??o necessariamente o irm??o de x,
		*	a ideia ?? pegar o  elemento anterior do memso n??vel n??vel)
		* 	com a mesma profundidade de x
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x
		*	@return : ??ndice do primeiro n?? a esquerda de x, com a mesma profundidade de x; ou menos 1 se a resposta n??o for encontrada.
		*/
		long long int levelPrev(long long int x);

		/*!
		* 	@brief usa fwdSearch para encontrar o n?? mais a esquerda com profundidade d
		*	@param d: profundidade desejada
		*	@return ??ndice do par??nteses de abertura do n?? mais a esquerda da ??rvore com profundidade d
		*/
		long long int levelLeftMost(int d);

		/*!
		* 	@brief usa bwdSearch (a partir de bv.size()) para encontrar o n?? mais a direita com profundidade d
		*	@param d: profundidade desejada
		*	@return ??ndice do par??nteses de abertura do n?? mais a direita da ??rvore com profundidade d
		*/
		long long int levelRightMost(int d);

		/*!
		* 	@brief usa rMq para buscar o filho de x que tenha a maior profundidade (mais a esquerda). 
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x
		*	@return : ??ndice do filho de x com maior profundidade
		*/
		long long int  deepestNode(long long int x);

		/*!
		* 	@brief usa minCount para computar o n??mero de vezes que o excesso m??nimo (ocorre quando x engloba um filho)
		*	aparece em x. 
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x
		*	@return : quantidade de filhos de x
		*/
		long long int  degree(long long int x);

		/*!
		*	@brief contabiliza a quantida de ocorr??ncias do bit 1 seguido do bit 0, no intervalo B[0,x].
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x
		*	@return: quantidade de folhas a esquerda de x
		*/
		long long int leafRank(long long int x);

		/*!
		*	@brief busca o ??ndice da t-th ocorr??ncia do par de bits 10 (ou seja da t-th folha).
		*	@param t: ordem da folha buscada
		*	@return: ??ndice da t-th folha.
		*/
		long long int leafSelect(long long int t);

		/*!
		* 	@brief usa leafSelect para encontrar a folha mais a esquerda de x.
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x
		*	@return ??ndice da folha mais ?? esquerda de x. Se x for uma folha retorna seu irm??o mais a esquerda,
		*	se ele n??o existir, retorna o pr??prio x. Se x for um n??? pai, retorna o seu filho com maior profundidade,mais a esqueda.
		*/
		long long int leftMostLeaf(long long int x);

		/*!
		* 	@brief usa leafSelect para encontrar a folha mais a direita de x.
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x
		*	@return ??ndice da folha mais ?? direita de x. Se x for uma folha retorna seu irm??o mais a direita,
		*	se ele n??o existir, retorna o pr??prio x. Se x for um n??? pai, retorna o seu filho com maior profundidade,mais a direita.
		*/
		long long int rightMostLeaf(long long int x);

		/*!
		*	@brief visita x usando um percurso pr??-ordem
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x
		*	@return rank1 no intervalo B[0,x]
		*/
		long long int preRank(long long int x);

		/*!
		*	@brief calcula o rank de x, a partir de um percurso pos-order
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x
		*	@return rank0 no intervalo B[0,close(x)]
		*/
		long long int postRank(long long int x);

		/*!
		*	@brief calcula a posi????o do t-th par??nteses de abertura
		*	@param x: par??nteses de abertura no vetor de par??nteses balanceados que codifica o n?? x
		*	@return t-th par??nteses de abertura
		*/
		long long int preSelect(long long int t);

		long long int postSelect(long long int t);

	private:
		rank_support_v<1> b_rank1;			// Fornece suporte a opera????op  rank, tendo como alvo bit 1
		rank_support_v<0> b_rank0;			// Fornece suporte a opera????op  rank, tendo como alvo bit 0
		rank_support_v<10,2> b_rank10;		// Fornece suporte a opera????op  rank, tendo como alvo a ocorr??ncia do bit 1,seguido do bit 0
		select_support_mcl<1> b_sel1;		// Fornece suporte a opera????op  select, tendo como alvo bit 1
		select_support_mcl<0> b_sel0;		// Fornece suporte a opera????op  select, tendo como alvo bit 0
		select_support_mcl<10,2> b_sel10;	// Fornece suporte a opera????op  select, tendo como alvo a ocorr??ncia do bit 1,seguido do bit 0
		int sizeBlock;						// Tamanho do intervalocoberto por um n?? folha
		int w;								// Divisor de sizeBlock. usado para pecorrer os bits de bv, de w em w, e assim acelerar o processo.
		long int numberLeaves;				// Quantidade de folhas na rmM-tree
		long int numberNodes;				// N??mero de n??s da rmM-tree
		int height;							// Altura da rmM-tree
		vector<Node_bin> tableC;			// Tabela de bits, com valores de excesso pr??-computados,usados para acelar a constru????o da rmM-tree
		
		/*!
		*	@brief reverte os bits lidos por bitsread
		*/
		uint16_t reverse_16(uint16_t x);

		/*!
		*	@return retorna a se a < b e b se b <=a. 
		*/
		long long int min(long long int a , long long int b);

		/*!
		*	@brief Pr??-computa uma tabela de excessos C, para agilizar a constru????o e as opera????es da RMM-tree.
		*/
		void buildingTableC();

		/*!
		*	@brief Constro?? as folhas de cada n??vel da RmM-tree.
		*/
		void buildingLeaves();

		/*!
		*	@brief Constro?? os n??s internos e a ra??z da RmM-tree
		*/
		void buildingInternalNodesRoot();

		/*!
		*	@brief Imprime as informa????es de excesso de um n??
		*	@param vector: estrutura (??rvore ou tabela) que ter?? seu n?? impresso
		*	@param i: ??ndice do elemento da estrutura.
		*/
		void printNode(vector<Node_bin> vector, long long int i);
		
		/*!
		*	@brief Pecorre para frente cada subbloco de tamanho "w" do bloco pertencente ?? "i".
		*	@param i: Posi????o a partir da qual devo buscar o excesso (i ?? adicionado de 1)
		*	@param d: Excesso buscado
		*	@param dr: Excesso relativo (atualizado a cada posi????o que avan??amos no bloco)
		*	@return a posi????o em que ocorre o excesso d ou bv.size() caso o excesso n??o se encontre neste bloco.
		*/
		long long int fwdBlock(long long int i,int d,int &dr);

		/*!
		*	@brief Pecorre para tr??s cada subbloco de tamanho "w" do bloco pertencente ?? "i".
		*	@param i: Posi????o a partir da qual devo iniciar a busca para tr??s do excesso
		*	@param d: Excesso buscado
		*	@param dr: Excesso relativo (atualizado a cada posi????o que avan??amos no bloco)
		*	@return a posi????o em que ocorre o excesso d ou ou bv.size() caso o excesso n??o se encontre neste bloco.
		*/
		long long int bwdBlock(long long int i,int d, int &dr);

		/*!
		*	@brief Pecorre para frente cada subbloco de tamanho "w" do bloco pertencente ?? "i", em busca do menor excesso na ??rea.
		*	@param i: Posi????o a partir da qual devo iniciar a busca.
		*	@param j: Intervalo superior da busca. Varremos at?? "j" ou at?? chegarmos ao limite do bloco de i, o que vier primeiro.
		*	@param d: Excesso relativo.
		*	@return o excesso m??nimo no intervalo definido.
		*/
		long long int minBlock(long long int i,long long int j,int &d);

		/*!
		*	@brief Pecorre para frente cada subbloco de tamanho "w" do bloco pertencente ?? "i", em busca do maior excesso na ??rea.
		*	@param i: Posi????o a partir da qual devo iniciar a busca.
		*	@param j: Intervalo superior da busca. Varremos at?? "j" ou at?? chegarmos ao limite do bloco de i, o que vier primeiro.
		*	@param d: Excesso relativo.
		*	@return o excesso m??ximo no intervalo definido.
		*/
		long long int maxBlock(long long int i,long long int j,int &d);

		/*!
		*	@brief Contabiliza a quantidade de vezes que o excesso m??nimo ocorre no intervalo i e j
		*	@param i: Posi????o a partir da qual devo iniciar a busca.
		*	@param j: Intervalo superior da busca (limite superior da folha de i).
		*	@param m : excesso m??nimo a ser computado.
		*	@param d : excesso relativo (contabilizado a cada bit varrido), par??metro passado por refer??ncia
		*	@return o n??mero de vezes que o excesso m??nimo aparece no intervalo definido.
		*/
		long long int minCountBlock(long long int i,long long int j,long long int m,int &d);

		/*!
		*	@brief Procura a posi????o em que o excesso m??nimo ocorre pela t-th no intervalo i e j
		*	@param i: Posi????o a partir da qual devo iniciar a busca.
		*	@param j: Intervalo superior da busca (limite superior da folha de i).
		*	@param m : excesso m??nimo a ser computado.
		*	@param t: t-th excesso m??nimo que quero encontrar
		*	@param d : excesso relativo (contabilizado a cada bit varrido), par??metro passado por refer??ncia
		*	@return a posi????o p em que o corre o t-th excesso m??nimo m.
		*/
		long long int minSelectBlock(long long int i,long long int j,long long int m,long long int &t, int &d);

};

#endif
