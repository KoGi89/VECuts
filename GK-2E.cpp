//usage: ./CutEdgesHigh <input_graph> <output_file>

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

int n, m;
int* edges;
int* G; int* firstOut;

void create_adj();
void find_cut_edges();
int* cut_edges;
int n_cut_edges;
int n_cut_pairs;

int* dfs, * idfs, * p;
int* high;
void compute_high();
int* ufparent, * ufrank, * representative;
void ufInit();
int find(int);
void unite(int, int, int);

int main(int n_args, char** args){
	FILE* fp = fopen(args[1], "r");
	fscanf(fp, "%d %d", &n, &m);
	edges = (int*)malloc(sizeof(int)*2*m);
	for(int i=0;i<m;i++){
		int x, y;
		fscanf(fp, "%d %d", &x, &y);
		edges[2*i]=x-1; edges[2*i+1]=y-1;
	}
	fclose(fp);

	create_adj();

	struct timeval begin, end;
	gettimeofday(&begin, 0);
	find_cut_edges();

	gettimeofday(&end, 0);
	long seconds = end.tv_sec - begin.tv_sec;
	long microseconds = end.tv_usec - begin.tv_usec;
	double elapsed = seconds + microseconds*1e-6;
	printf("Total time= %g\n", elapsed);

	fp = fopen(args[2], "w");
	fprintf(fp, "%d\n", n_cut_edges);
	for(int i=0;i<n_cut_edges;i++){
		fprintf(fp, "%d %d\n", cut_edges[2*i]+1, cut_edges[2*i+1]+1);
	}

	fprintf(fp, "%d\n", n_cut_pairs);
	fprintf(fp, "%lf\n", elapsed);
	fclose(fp);

	return 0;
}

void find_cut_edges(){
	cut_edges = (int*)malloc(sizeof(int)*(4*n-4));
	n_cut_edges = 0;

	dfs = (int*)malloc(sizeof(int)*n);
	idfs = (int*)malloc(sizeof(int)*n);
	p = (int*)malloc(sizeof(int)*n);
	int* l = (int*)malloc(sizeof(int)*n);
	int* low = (int*)malloc(sizeof(int)*n);
	int* currentChild = (int*)malloc(sizeof(int)*n);
	int* nextSibling = (int*)malloc(sizeof(int)*n);
	int* prevSibling = (int*)malloc(sizeof(int)*n);
	int* L = (int*)malloc(sizeof(int)*n);
	int* R = (int*)malloc(sizeof(int)*n);
	int* M = (int*)malloc(sizeof(int)*n);
	int* nextM = (int*)malloc(sizeof(int)*n);
	int* bCount = (int*)malloc(sizeof(int)*n);

	int* stack = (int*)malloc(sizeof(int)*n);
	int* temp_out = (int*)malloc(sizeof(int)*n);

	for(int i=0;i<n;i++){
		dfs[i]=-1; l[i]=i; low[i]=i; currentChild[i]=-1; bCount[i]=0; nextSibling[i]=-1;prevSibling[i]=-1;
	}

	//perform DFS
	int Nr=0;
	dfs[0]=Nr; idfs[Nr++]=0;
	stack[0]=0; temp_out[0]=firstOut[0];
	int SP=0;
	while(SP!=-1){
		int v=stack[SP];
		char down=0;
		for(int i=temp_out[SP];i<firstOut[v+1];i++){
			int u=G[i];
			if(dfs[u]==-1){
				dfs[u]=Nr; idfs[Nr++]=u; p[u]=v;
				if(currentChild[v]!=-1){
					nextSibling[currentChild[v]]=u;
					prevSibling[u]=currentChild[v];
				}
				currentChild[v]=u;
				stack[SP+1]=u; temp_out[SP+1]=firstOut[u]; temp_out[SP]=i;
				down=1; break;
			}
			if(dfs[u]<dfs[v] && u!=p[v]){
				if(dfs[u]<dfs[l[v]]){ l[v]=u; }
				if(dfs[u]<dfs[low[v]]){ low[v]=u; }
				bCount[v]++; bCount[u]--;
			} else if(v==p[u]){
				if(dfs[low[u]]<dfs[low[v]]){ low[v]=low[u]; }
				bCount[v]+=bCount[u];
			}
		}
		if(down){ SP++;continue; }
		SP--;
	}

	//compute the high points
	compute_high();

	//initialize all L and R
	for(int v=0;v<n;v++){ L[v]=-1; }
	for(int i=1;i<n;i++){
		int u=idfs[i];
		if(dfs[low[u]]<dfs[p[u]]){
			if(L[p[u]]==-1){
				L[p[u]]=u; R[p[u]]=u;
			} else{
				R[p[u]]=u;
			}
		}
	}

	//calculate M and nextM
	for(int v=0;v<n;v++){ nextM[v]=-1; }
	for(int i=n-1;i>0;i--){
		int v=idfs[i];
		if(l[v]!=v){ M[v]=v;continue; }
		if(L[v]!=R[v]){ M[v]=v;continue; }
		int c=L[v];
		int m=M[c];
		while(1){
			if(dfs[l[m]]<i){ M[v]=m;break; }
			while(dfs[low[L[m]]]>=i){ L[m]=nextSibling[L[m]]; }
			while(dfs[low[R[m]]]>=i){ R[m]=prevSibling[R[m]]; }
			if(L[m]!=R[m]){ M[v]=m;break; }
			c=L[m];
			m=M[c];
		}
		nextM[c]=v;
	}

	char* tree_edge_is_cut_edge = (char*)malloc(sizeof(char)*n);
	char* back_edge_is_cut_edge = (char*)malloc(sizeof(char)*n);

	for(int i=0;i<n;i++){
		tree_edge_is_cut_edge[i]=0; back_edge_is_cut_edge[i]=0;
	}

	//find cut-edges
	for(int v=1;v<n;v++){
		if(bCount[v]==1){
			back_edge_is_cut_edge[M[v]]=1;
			tree_edge_is_cut_edge[v]=1;
		}
		int next=nextM[v];
		if(next!=-1 && high[v]==high[next]){
			tree_edge_is_cut_edge[v]=1;
			tree_edge_is_cut_edge[next]=1;
		}
	}

	//count cut-pairs
	n_cut_pairs=0;
	for(int v=1;v<n;v++){
		if(bCount[v]==1){
			n_cut_pairs++;
		}
	}
	for(int v=1;v<n;v++){
		if(v!=M[v]){ continue; }
		int u=v;
		while(u!=-1){
			int next=nextM[u];
			if(next==-1){ break; }
			int n_edges=0;
			int h=high[u];
			while(h==high[next]){
				n_edges++;
				next=nextM[next];
				if(next==-1){ break; }
			}
			n_cut_pairs+=(n_edges*(n_edges+1))/2;
			u=next;
		}
	}

	for(int v=1;v<n;v++){
		if(tree_edge_is_cut_edge[v]){
			cut_edges[2*n_cut_edges]=v;
			cut_edges[2*n_cut_edges+1]=p[v];
			n_cut_edges++;
		}
		if(back_edge_is_cut_edge[v]){
			cut_edges[2*n_cut_edges]=v;
			cut_edges[2*n_cut_edges+1]=low[v];
			n_cut_edges++;
		}
	}
}

void compute_high(){
	high = (int*)malloc(sizeof(int)*n);
	ufInit();
	for(int i=n-1;i>=0;i--){
		int v=idfs[i];
		for(int j=firstOut[v];j<firstOut[v+1];j++){
			int u=G[j];
			if(dfs[v]<dfs[u] && v!=p[u]){
				int x = representative[find(u)];
				while(x!=v){
					high[x] = v;
					int next = representative[find(p[x])];
					unite(x, p[x], next);
					x = next;
				}
			}
		}
	}
}

void ufInit(){
	ufparent = (int*)malloc(sizeof(int)*n);
	ufrank = (int*)malloc(sizeof(int)*n);
	representative = (int*)malloc(sizeof(int)*n);
	for(int i=0;i<n;i++){ ufparent[i]=i; ufrank[i]=0; representative[i]=i; }
}

int find(int x){
	int r=x;
	while(ufparent[r]!=r){ r=ufparent[r]; }
	while(ufparent[x]!=x){ int next=ufparent[x]; ufparent[x]=r; x=next; }
	return r;
}

void unite(int x, int y, int w){
	int r1=find(x);
	int r2=find(y);
	if(r1==r2){ representative[r1]=w;return; }
	int cmp=ufrank[r1]-ufrank[r2];
	if(cmp<0){
		ufparent[r1]=r2;
		representative[r2]=w;
	} else if(cmp>0){
		ufparent[r2]=r1;
		representative[r1]=w;
	} else{
		ufparent[r1]=r2;
		ufrank[r2]++;
		representative[r2]=w;
	}
}

void create_adj(){
	G = (int*)malloc(sizeof(int)*4*m);
	firstOut = (int*)malloc(sizeof(int)*(n+1));
	for(int i=0;i<=n;i++){ firstOut[i]=0; }
	for(int i=0;i<m;i++){ firstOut[edges[2*i]+1]++; firstOut[edges[2*i+1]+1]++; }
	int* nextOut = (int*)malloc(sizeof(int)*(n+1));
	nextOut[0]=0;
	for(int i=1;i<=n;i++){ firstOut[i]+=firstOut[i-1]; nextOut[i]=firstOut[i]; }
	for(int i=0;i<m;i++){
		int x=edges[2*i]; int y=edges[2*i+1];
		G[nextOut[x]++]=y;
		G[nextOut[y]++]=x;
	}
}










