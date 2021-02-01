//usage: ./LinearVE <input_graph> <output_file>

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

int n, m;
int* edges;
int* G; int* firstOut;

void create_adj();
void compute_count();
int* count;

int* dfs, * idfs, * parent, * T, * firstChild, * tempChild;
int* l, * low, * highp, * bcount, * up, * L, * R, * M, * Mp;
int* invHighp, * firstInvHighp;
int* invM, * invMfirst, * invMp, * invMpfirst;


void construct_tree();
void compute_M();
void computeInvM();
void compute_highp();
void sort_ChildrenHighp();

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
	compute_count();
	gettimeofday(&end, 0);
	long seconds = end.tv_sec - begin.tv_sec;
	long microseconds = end.tv_usec - begin.tv_usec;
	double elapsed = seconds + microseconds*1e-6;
	printf("Total time= %g\n", elapsed);

	fp = fopen(args[2], "w");
	for(int i=0;i<n;i++){
		fprintf(fp, "%d\n", count[i]);
	}
	fprintf(fp, "%lf\n", elapsed);
	fclose(fp);

	return 0;
}

void compute_count(){
	count = (int*)malloc(sizeof(int)*n);
	for(int i=0;i<n;i++){ count[i]=0; }

	dfs = (int*)malloc(sizeof(int)*n);
	idfs = (int*)malloc(sizeof(int)*n);
	parent = (int*)malloc(sizeof(int)*n);
	tempChild = (int*)malloc(sizeof(int)*n);
	l = (int*)malloc(sizeof(int)*n);
	low = (int*)malloc(sizeof(int)*n);
	bcount = (int*)malloc(sizeof(int)*n);
	up = (int*)malloc(sizeof(int)*n); //up[v]=#{back-edges (u,p(v)), where u is a descendant of v}

	//perform DFS
	for(int i=0;i<n;i++){ dfs[i]=-1; l[i]=i; low[i]=i; bcount[i]=0; up[i]=0; }
	int* temp_vertex = (int*)malloc(sizeof(int)*n);
	int* temp_out = (int*)malloc(sizeof(int)*n);
	int Nr=0;
	dfs[0]=Nr; idfs[Nr++]=0; parent[0]=-1;
	temp_vertex[0]=0; temp_out[0]=firstOut[0];
	int SP=0;
	while(SP!=-1){
		int v=temp_vertex[SP];
		char descend=0;
		for(int i=temp_out[SP];i<firstOut[v+1];i++){
			int u=G[i];
			if(dfs[u]==-1){
				dfs[u]=Nr; idfs[Nr++]=u; parent[u]=v; tempChild[v]=u;
				temp_vertex[SP+1]=u; temp_out[SP+1]=firstOut[u]; temp_out[SP]=i;
				descend=1; break;
			}
			if(dfs[u]<dfs[v] && u!=parent[v]){
				if(dfs[u]<dfs[l[v]]){
					l[v]=u;
					if(dfs[u]<dfs[low[v]]){
						low[v]=u;
					}
				}
				bcount[v]++; up[tempChild[u]]++;
			} else if(v==parent[u]){
				if(dfs[low[u]]<dfs[low[v]]){
					low[v]=low[u];
				}
				bcount[v]+=bcount[u];
			}
		}
		if(descend){ SP++;continue; }
		bcount[v]-=up[v];
		SP--;
	}

	//gather parameters
	construct_tree();
	compute_M();
	computeInvM();
	compute_highp();

	//e is a back-edge
	for(int i=2;i<n;i++){
		int v=idfs[i];
		count[parent[v]]+=bcount[v]==1;
	}

	//high(u)=v
	for(int v=0;v<n;v++){
		int u=firstInvHighp[v];
		int c=firstChild[v];
		int endH=firstInvHighp[v+1];
		int endT=firstChild[v+1];
		while(u!=endH){
			if(up[invHighp[u]]>0){ u++;continue; } //check if high(u)=highp(u)
			while(!(c==endT-1 || dfs[T[c+1]]>dfs[invHighp[u]])){
				c++;
			}
			if(low[invHighp[u]]==v || dfs[invHighp[u]]<=dfs[Mp[T[c]]]){
				count[v]++;
			}
			u++;
		}
	}

	//high(u)<v
	for(int x=0;x<n;x++){
		int u=invMfirst[x];
		int c=invMpfirst[x];
		int endM=invMfirst[x+1];
		int endMp=invMpfirst[x+1];
		while(u!=endM && c!=endMp){
			while(c!=endMp && dfs[invMp[c]]>=dfs[invM[u]]){ c++; }
			if(c==endMp){ break; }
			if(up[invM[u]]==0 && dfs[highp[invM[u]]]<dfs[parent[invMp[c]]]){
				int n_edges=0;
				int h=dfs[highp[invM[u]]];
				while(c!=endMp && h<dfs[parent[invMp[c]]]){
					while(u!=endM && dfs[invMp[c]]<dfs[invM[u]]){
						n_edges++;
						u++;
					}
					count[parent[invMp[c]]]+=n_edges;
					c++;
				}
			} else{
				u++;
			}
		}
	}

	//M(u)=v
	sort_ChildrenHighp(); //sort the lists of children in decreasing order w.r.t. the highp of their elements

	for(int v=1;v<n;v++){
		if(invMfirst[v]==invMfirst[v+1]){ continue; }
		int u=invMfirst[v]+1;
		int c=firstChild[v];
		int endM=invMfirst[v+1];
		int endT=firstChild[v+1];
		int min=v;
		while(u!=endM && c!=endT){
			min=highp[T[c]];
			while(u!=endM && dfs[invM[u]]>dfs[min]){
				count[v]++;
				u++;
			}
			min=low[T[c]];
			c++;
			while(c!=endT && dfs[highp[T[c]]]>=dfs[min]){
				if(dfs[low[T[c]]]<dfs[min]){
					min=low[T[c]];
				}
				c++;
			}
			while(u!=endM && dfs[invM[u]]>dfs[min]){ u++; }
		}
		while(u!=endM){
			if(dfs[invM[u]]<=dfs[min]){
				count[v]++;
			}
			u++;
		}
	}

	//M(u)>v
	for(int x=0;x<n;x++){
		int c = invMpfirst[x];
		int u = invMfirst[x];
		int endMp = invMpfirst[x+1];
		int endM = invMfirst[x+1];
		while(c!=endMp && u!=endM){
			while(u!=endM && dfs[invM[u]]>=dfs[parent[invMp[c]]]){ u++; }
			if(u==endM){ break; }
			if(dfs[highp[invMp[c]]]<dfs[invM[u]]){
				int n_edges=0;
				int first = u;
				while(u!=endM && dfs[highp[invMp[c]]]<dfs[invM[u]]){
					n_edges++;
					u++;
				}
				int last = u-1;
				count[parent[invMp[c]]]+=n_edges;
				c++;
				while(c!=endMp && dfs[invM[last]]<dfs[parent[invMp[c]]]){
					while(dfs[invM[first]]>=dfs[parent[invMp[c]]]){
						n_edges--;
						first++;
					}
					count[parent[invMp[c]]]+=n_edges;
					c++;
				}
			} else{
				c++;
			}
		}
	}
}

void construct_tree(){
	//sort the list of the children of every vertex in increasing order
	T = (int*)malloc(sizeof(int)*n);
	firstChild = (int*)malloc(sizeof(int)*(n+1));
	for(int i=0;i<=n;i++){ firstChild[i]=0; }
	for(int i=1;i<n;i++){ firstChild[parent[idfs[i]]+1]++; }
	for(int i=0;i<n;i++){ tempChild[i]=0; }
	int tc_r=0;
	tempChild[0]=firstChild[1];
	for(int i=1;i<n;i++){ firstChild[i+1]+=firstChild[i]; tempChild[i]=firstChild[i+1]; }
	for(int i=1;i<n;i++){
		int v=idfs[i];
		if(parent[v]==0){ T[tc_r++]=v; } else{ T[tempChild[parent[v]-1]++]=v; }
	}
}

void compute_M(){
	L = (int*)malloc(sizeof(int)*n);
	R = (int*)malloc(sizeof(int)*n);
	M = (int*)malloc(sizeof(int)*n);
	Mp = (int*)malloc(sizeof(int)*n);

	//calculate all M and Mp
	for(int i=n-1;i>0;i--){
		int v=idfs[i];

		//initialize L and R to calculate M[v]
		L[v]=firstChild[v];
		R[v]=firstChild[v+1]-1;

		//compute M
		if(dfs[l[v]]<dfs[v]){ M[v]=v; } else if(L[v]!=R[v]){ M[v]=v; } else{
			int c=T[L[v]];
			int m=M[c];
			while(1){
				if(dfs[l[m]]<dfs[v]){ M[v]=m;break; }
				while(dfs[low[T[L[m]]]]>=dfs[v]){ L[m]++; }
				while(dfs[low[T[R[m]]]]>=dfs[v]){ R[m]--; }
				if(L[m]!=R[m]){ M[v]=m;break; }
				c=T[L[m]];
				m=M[c];
			}
		}

		if(i==1){ break; } //Mp is undefined for the child of the root

		//compute Mp
		if(dfs[l[v]]<dfs[parent[v]]){ Mp[v]=v;continue; }

		//update L[v] and R[v] to calculate Mp[v]
		while(dfs[low[T[L[v]]]]>=dfs[parent[v]]){ L[v]++; }
		while(dfs[low[T[R[v]]]]>=dfs[parent[v]]){ R[v]--; }

		if(L[v]!=R[v]){ Mp[v]=v; } else{
			int c=T[L[v]];
			int m=Mp[c];
			while(1){
				if(dfs[l[m]]<dfs[parent[v]]){ Mp[v]=m;break; }
				while(dfs[low[T[L[m]]]]>=dfs[parent[v]]){ L[m]++; }
				while(dfs[low[T[R[m]]]]>=dfs[parent[v]]){ R[m]--; }
				if(L[m]!=R[m]){ Mp[v]=m;break; }
				c=T[L[m]];
				m=Mp[c];
			}
		}
	}
}

void computeInvM(){
	//calculate the inverse lists, and have their elements sorted in decreasing order
	invM = new int[n];
	invMfirst = new int[n+1];
	for(int i=0;i<=n;i++){ invMfirst[i]=0; }
	for(int i=1;i<n;i++){
		int v=idfs[i];
		invMfirst[M[v]+1]++;
	}
	for(int i=1;i<=n;i++){ invMfirst[i]+=invMfirst[i-1]; }
	int* invMnext = new int[n+1];
	for(int i=0;i<=n;i++){ invMnext[i]=invMfirst[i]; }
	for(int i=n-1;i>0;i--){
		int v=idfs[i];
		invM[invMnext[M[v]]++]=v;
	}

	invMp = new int[n];
	invMpfirst = new int[n+1];
	for(int i=0;i<=n;i++){ invMpfirst[i]=0; }
	for(int i=2;i<n;i++){
		int v=idfs[i];
		invMpfirst[Mp[v]+1]++;
	}
	for(int i=1;i<=n;i++){ invMpfirst[i]+=invMpfirst[i-1]; }
	for(int i=0;i<=n;i++){ invMnext[i]=invMpfirst[i]; }
	for(int i=n-1;i>1;i--){
		int v=idfs[i];
		invMp[invMnext[Mp[v]]++]=v;
	}
}

void compute_highp(){
	highp = (int*)malloc(sizeof(int)*n);
	ufInit();
	for(int i=n-1;i>=0;i--){
		int v=idfs[i];
		for(int j=firstOut[v];j<firstOut[v+1];j++){
			int u=G[j];
			if(dfs[v]<dfs[u] && v!=parent[u]){
				int x = representative[find(u)];
				while(parent[x]!=v){
					highp[x] = v;
					int next = representative[find(parent[x])];
					unite(x, parent[x], next);
					x = next;
				}
			}
		}
	}

	//sort the elements of the inverse lists highp^-1 in increasing order
	invHighp = (int*)malloc(sizeof(int)*n);
	firstInvHighp = (int*)malloc(sizeof(int)*(n+1));
	for(int i=0;i<=n;i++){ firstInvHighp[i]=0; }
	for(int i=2;i<n;i++){ firstInvHighp[highp[idfs[i]]+1]++; }
	int temp0=0;
	tempChild[0]=firstInvHighp[1];
	for(int i=1;i<n;i++){ firstInvHighp[i+1]+=firstInvHighp[i]; tempChild[i]=firstInvHighp[i+1]; }
	for(int i=2;i<n;i++){
		int v=idfs[i];
		if(highp[v]==0){ invHighp[temp0++]=v; } else{ invHighp[tempChild[highp[v]-1]++]=v; }
	}
}

void sort_ChildrenHighp(){
	for(int i=0;i<n;i++){ tempChild[i]=firstChild[i]; }
	for(int i=n-1;i>=0;i--){
		int v=idfs[i];
		for(int i=firstInvHighp[v];i<firstInvHighp[v+1];i++){
			int c=invHighp[i];
			T[tempChild[parent[c]]++]=c;
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










