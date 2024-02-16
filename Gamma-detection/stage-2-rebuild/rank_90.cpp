#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<algorithm>
#include<cmath>
#include<queue>
#include<vector>
#include<ctime>
#include<map>
#include<bitset>
#include<set>
#include<assert.h>
#include<chrono>
#include<random>
#include<iostream>
#define LL long long
#define mp(x,y) make_pair(x,y)
#define pll pair<long long,long long>
#define pii pair<int,int>
#define SZ(x) ((int)x.size())
#define VI vector<int>
#define ull unsigned long long
using namespace std;
mt19937 rnd(chrono::steady_clock::now().time_since_epoch().count());
inline LL read()
{
	LL f=1,x=0;char ch=getchar();
	while(ch<'0'||ch>'9'){if(ch=='-')f=-1;ch=getchar();}
	while(ch>='0'&&ch<='9'){x=x*10+ch-'0';ch=getchar();}
	return x*f;
}
int stack[20];
template<typename T>inline void write(T x)
{
	if(x<0){putchar('-');x=-x;}
    if(!x){putchar('0');return;}
    int top=0;
    while(x)stack[++top]=x%10,x/=10;
    while(top)putchar(stack[top--]+'0');
}
template<typename T>inline void pr1(T x){write(x);putchar(' ');}
template<typename T>inline void pr2(T x){write(x);putchar('\n');}
template<typename T>inline void chkmin(T &x,T y){x=x<y?x:y;}
template<typename T>inline void chkmax(T &x,T y){x=x>y?x:y;}
const double PI = acos(-1.0);
const int MAXN = 40;
struct vec
{
	double x,y,z;
	vec(){}
	vec(double _x,double _y,double _z){x=_x;y=_y;z=_z;}
};
double X[MAXN],Y[MAXN],Z[MAXN],E_test[MAXN],T_test[MAXN],id[MAXN];
double vec_dot(vec vec1, vec vec2){return vec1.x*vec2.x+vec1.y*vec2.y+vec1.z*vec2.z;}
double calculate_theta(int pre, int now, int suf) // ������Լнǵ�cos_thetaֵ 
{
	vec vec1 = vec(X[now]-X[pre], Y[now]-Y[pre], Z[now]-Z[pre]);
	vec vec2 = vec(X[suf]-X[now], Y[suf]-Y[now], Z[suf]-Z[now]);
	double cos_theta = vec_dot(vec1, vec2) / ( sqrt(vec_dot(vec1, vec1)) * sqrt(vec_dot(vec2, vec2)) );
	return cos_theta;
}

double calculate_prob(double x_std, double sig_std, double x_test) // ��̬�ֲ� N(x_std, sig_std) ȡx_test�ĸ��� 
{
	double up1 = - (x_test - x_std) * (x_test - x_std);
	up1 /= (2 * sig_std * sig_std);
	up1 = exp(up1);
	up1 /= (sqrt(2 * PI) * sig_std);
	return up1;
}
const double A_energy_sigma = 0.01873269;
const double B_energy_sigma = 0.77153127;
const double C_energy_sigma = 0.01036436;
double calculate_energy_sigma(double E_std)
{
	return E_std * (A_energy_sigma * exp(- B_energy_sigma * E_std) + C_energy_sigma);
}
int n;

int P[MAXN], vis[MAXN], P_ans[MAXN];

double ans;

int exit_limit;

void dfs(int k, double fval, double tval, double eval)
{
	if(fval < exit_limit) return ;
	if(k == n)
	{
		double E_final_sigma = calculate_energy_sigma(eval);// ���������ЧӦ�ĸ���
		
		double prob = calculate_prob(eval, E_final_sigma, E_test[P[n]]);
		fval = fval + log10(prob);
		
		if(fval > ans)
		{
			ans = fval;
			for(int i=1;i<=n;i++)P_ans[i] = P[i];
		}
		return ;
	}
	for(int i=1;i<=n;i++)if(!vis[i])
	{
		vec vec_relative = vec(X[i]-X[P[k]], Y[i]-Y[P[k]], Z[i]-Z[P[k]]);
		double vec_len = sqrt(vec_dot(vec_relative, vec_relative)); // ��λ�� mm 
		double t_cost = vec_len / 299.792458; // �˴��Ѿ��ѵ�λͳһ���� ns ( C * 1e-9 * 1e3 ) 
		
		double t_next = tval + t_cost;
		double prob_t = log10(calculate_prob(t_next, 1.5, T_test[i]));
		
		if(prob_t < exit_limit) continue;
		double f_next = fval + prob_t, E_next = eval;
		
		if(k!=0) // ���� P[k-1] P[k] �� i �Ŀ��ն�ɢ����� 
		{
			double photon_scatter_cos_theta = calculate_theta(P[k-1], P[k], i); // �������ɢ��� 
			E_next = eval / (1 + eval / 0.511 * (1 - photon_scatter_cos_theta) );
			double E_scatter_elec_exact = eval - E_next;
			
			double E_scatter_elec_sigma = calculate_energy_sigma(E_scatter_elec_exact);
			
			double prob_e = log10(calculate_prob(E_scatter_elec_exact, E_scatter_elec_sigma, E_test[P[k]]));
			if(prob_e < exit_limit) continue;
			f_next += prob_e;
		}
		
		P[k+1] = i; vis[i] = 1;
		dfs(k+1, f_next, t_next, E_next);
		P[k+1] = 0; vis[i] = 0;
	}
}

int main()
{
	#ifdef Rose
		double BeginJudgeTime=clock();
	#endif
	
	int liny;
	n=read();
	for(int i=1;i<=n;i++)
		scanf("%d%lf%lf%lf%lf%lf",&id[i],&X[i],&Y[i],&Z[i],&E_test[i],&T_test[i]);
	
	// ��֦ 
	exit_limit = -90;
	while(1)
	{ 
		ans = -999999999;
		
		dfs(0, 0, 0, 0.662);
	
		if(ans < - 2000) exit_limit -=100;
		else if(exit_limit < - 500)
		{
			cerr<<"����Ҳ�����"<<endl;
			break;
		}
		else break;
	}
	cerr<<ans<<endl;
	for(int i=1;i<=n;i++)
	{
		for(int j=1;j<=n;j++)if(P_ans[j]==i)pr2(j-1);
	}
	
	
	#ifdef Rose
		double EndJudgeTime=clock();
		cerr<<"JudgeTime is"<<" ";
		cerr<<(EndJudgeTime-BeginJudgeTime)/CLOCKS_PER_SEC<<endl;
	#endif
	return 0;
}

