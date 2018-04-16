#include "api.h"

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZ_p.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>

using namespace std;
using namespace NTL;

#if N == 1024
	const long CDT_length = 11;
	const unsigned long CDT[22] = {0, 4402564254475628998UL, 11764982215938697676UL, 16069007141612633800UL, 17828139942285597337UL, 18330820414912403750UL, 18431248444640379412UL, 18445276175012283164UL, 18446646076749434596UL, 18446739608909152776UL, 0, 0, 14371828625661980934UL, 4718557133137876784UL, 11277549837058342664UL, 11589062832341743150UL, 11926477327048184144UL, 3913014129104688847UL, 2742149467003257727UL, 4527561242905591802UL, 14730590305512100561UL, 0};
	const long CDT_inv_min[256] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 5};
	const long CDT_inv_max[256] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6};
#else
	const long CDT_length = 7;
	const unsigned long CDT[14] = {0, 7094901854892740860UL, 16010399019983881949UL, 18221726696665494136UL, 18438247182874980419UL, 18446616369313222554UL, 0, 0, 1132165336062584669UL, 8115838531820800606UL, 8387848400301917718UL, 11300671729619518355UL, 17590659254897005002UL, 0};
	const long CDT_inv_min[256] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3};
	const long CDT_inv_max[256] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4};
#endif


static unsigned char filter[8] = {1, 2, 4, 8, 16, 32, 64, 128};
static int bias[17] = {0, 445, 888, 1333, 1776, 2221, 2666, 3109, 3554, 3997, 4442, 4885, 5330, 5775, 6218, 6663, 7106};
unsigned long step = sizeof(unsigned long)*8;


void Sample(NTL::vec_ZZ_p &ret, NTL::RR sigma, int flag);
//int randombytes(unsigned char *x, unsigned long long xlen);
void Sample(NTL::vec_ZZ_p &ret, NTL::RR sigma, int flag)
{
// Sampling
// flag = 1, return 1 time of sampled value; flag = 2, return 2 times of sampled value
      for(int zz = 0; zz < N; zz++)
	{
		long r0, min, max;
		unsigned char ch;

		r0 = NTL::RandomBits_ulong(8);
	

		min = CDT_inv_min[r0];
		max = CDT_inv_max[r0];

		if (max-min < 2)
		{
			long val = (NTL::RandomBits_ulong(1)) ? min : -min;
			conv(ret[zz], flag*val);
			continue;
		}

		unsigned long r1;
		unsigned long cur;
		unsigned long mask_index;
		unsigned long r2;


		mask_index = step - 8;
		r1 = ((unsigned long) r0) << mask_index;
		r2 = (0xFF << mask_index);
		cur = (min+max)/2;

		while(true)
		{
			if(r1 > CDT[cur])
				min = cur;
			else if(r1 < (CDT[cur] & r2))
				max = cur;
			else
			{
				if(!mask_index)
					break;
				mask_index-= 8;

				r2 |= 0xFF << mask_index;
				r1 |= NTL::RandomBits_ulong(8) << mask_index;
			}
			if(max-min < 2)
			{
				long val = (NTL::RandomBits_ulong(1)) ? min : -min;
				conv(ret[zz], flag*val);
				continue;
			}
			cur = (min+max)/2;
		}

		r2 = NTL::RandomBits_ulong(step);
		while(true)
		{
			if (r1 < CDT[cur] || ((r1 == CDT[cur]) && (r2 < CDT[cur+CDT_length])))
				max = cur;
			else
				max = cur;
			cur = (min+max)/2;
			if(max-min < 2)
			{
				long val = (NTL::RandomBits_ulong(1)) ? min : -min;
				conv(ret[zz], val);
				continue;
			}
		}
	}
}


void Compute(Vec <ZZ_p>&b_ZZp,ZZ_pX& a_px,Vec <ZZ_p>&s_ZZp,Vec <ZZ_p>&e_ZZp);
void GenerateSet (Vec <ZZ_p>& I0,Vec <ZZ_p>& I1,Vec <ZZ_p>& sum_0,Vec <ZZ_p>& sum_1,Vec <ZZ>& E,ZZ& q);
void CrossRounding (Vec <ZZ>& cross_ZZ,Vec <ZZ>& R,long& q);
void ModularRounding (Vec <ZZ>& round_ZZ,Vec <ZZ>& R,long& q);
void Rec(Vec <ZZ>& coeff_R1_e,Vec <ZZ>& ZZ_vec_cross,Vec <ZZ>&ZZ_vec_rec,Vec <ZZ_p>& sum_0,Vec <ZZ_p>& sum_1);
int main()
{
	ZZ q;
	cout<<"Please enter q= ";
  	cin >> q;
	long q_long;
	conv(q_long,q);
	ZZ_p::init(q); 
 	RR sigma;
        sigma = 2.6;

  	Vec <ZZ> coeff_R1,coeff_R1_e,ZZ_vec_cross,ZZ_vec_modular,ZZ_vec_rec,E;
  	Vec<ZZ_p> I0,I1,sum_0,sum_1;

	GenerateSet(I0,I1,sum_0,sum_1,E,q);//generate I0,I1,E,I0+E,I1+E

	long n;
	n = N;
	NTL::vec_ZZ_p s_ZZp(NTL::INIT_SIZE, n), e_ZZp(NTL::INIT_SIZE, n);
	NTL::vec_ZZ_p s1_ZZp(NTL::INIT_SIZE, n), e1_ZZp(NTL::INIT_SIZE, n);
	NTL::vec_ZZ_p e2_ZZp(NTL::INIT_SIZE, n),e3_ZZp(NTL::INIT_SIZE, n);
	NTL::div(sigma, sigma, sqrt(2*PI));	// sigma /= sqrt(2*pi)
	//sample s,e
	Sample(s_ZZp, sigma, 1);
	Sample(e_ZZp, sigma, 1);

	//sample s1,e1
	Sample(s1_ZZp, sigma, 1);
	Sample(e1_ZZp, sigma, 1);

	//sample e2
	Sample(e2_ZZp, sigma, 1);

	NTL::vec_ZZ_p b_ZZp(NTL::INIT_SIZE, n);
	NTL::vec_ZZ_p b1_ZZp(NTL::INIT_SIZE, n);
	NTL::vec_ZZ_p v_ZZp(NTL::INIT_SIZE, n);
	NTL::vec_ZZ_p w_ZZp(NTL::INIT_SIZE, n);

	NTL::vec_ZZ v_ZZ(NTL::INIT_SIZE, n);
	NTL::vec_ZZ w_ZZ(NTL::INIT_SIZE, n);
	NTL::vec_ZZ cross_v_ZZ(NTL::INIT_SIZE, n);
	NTL::vec_ZZ round_v_ZZ(NTL::INIT_SIZE, n);
	NTL::vec_ZZ rec_v_ZZ(NTL::INIT_SIZE, n);

	NTL::ZZ_pX a_px,b_px,b1_px;
	random(a_px, n);//generate a random and uniformly
  
	Compute(b_ZZp,a_px,s_ZZp,e_ZZp);//compute b=as+e
	Compute(b1_ZZp,a_px,s1_ZZp,e1_ZZp);//compute b1=as1+e1

	conv(b_px,b_ZZp);
	Compute(v_ZZp,b_px,s1_ZZp,e2_ZZp);//compute v=bs1+e2

	conv(v_ZZ,v_ZZp);
	CrossRounding (cross_v_ZZ,v_ZZ,q_long);//compute cross_v_ZZ=<v>
	ModularRounding (round_v_ZZ,v_ZZ,q_long);//compute round_v_ZZ=round(v)
	//cout<<"v_ZZ"<<v_ZZ<<endl;
	//cout<<"cross(v)="<<cross_v_ZZ<<endl;
	cout<<"round(v)="<<round_v_ZZ<<endl;
	conv(b1_px,b1_ZZp);
	Compute(w_ZZp,b1_px,s_ZZp,e3_ZZp);//compute w=b1s+e3,e3=0
	conv(w_ZZ,w_ZZp);
	Rec(w_ZZ,cross_v_ZZ,rec_v_ZZ,sum_0,sum_1);//compute rec_v_ZZ=rec(w,<v>)
	//cout<<"w_ZZ"<<w_ZZ<<endl;
	cout<<"rec(w,b)="<<rec_v_ZZ<<endl;
	/*cout<<"s="<<s_ZZp<<endl;
	cout<<"e="<<e_ZZp<<endl;
	cout<<"a="<<a_px<<endl;
	cout<<"s="<<s_ZZp<<endl;
	cout<<"e="<<e_ZZp<<endl;
	cout<<"b=as+e="<<b_ZZp<<endl;
	cout<<"s1="<<s1_ZZp<<endl;
	cout<<"e1="<<e1_ZZp<<endl;
	cout<<"b1=as1+e1="<<b1_ZZp<<endl;
	cout<<"v=bs1+e2="<<v_ZZp<<endl;
	cout<<"w=b1s="<<w_ZZp<<endl;*/
}

/*input :a,s,e
  output:b=as+e*/
void Compute(Vec <ZZ_p>&b_ZZp,ZZ_pX &a_px,Vec <ZZ_p>&s_ZZp,Vec <ZZ_p>&e_ZZp)
{
	long n;
	n = N;
	NTL::ZZ_pX s_px, temp_px;
	s_px = to_ZZ_pX(s_ZZp);

	NTL::ZZ_pX u;
	u.SetMaxLength(n);
	SetCoeff(u, 0, 1);
	SetCoeff(u, n, 1);//u=1+x^N
	NTL::ZZ_pXModulus F(u);

	MulMod(temp_px, a_px, s_px, F);//temp_pX = (a_px*s_px)%(1+x^n)
	//cout<<"temp_px="<<temp_px<<endl;

	NTL::vec_ZZ_p pk_ZZp(NTL::INIT_SIZE, n);

	conv(b_ZZp, VectorCopy(temp_px, n));
	add(b_ZZp, b_ZZp, e_ZZp);
}

/*
  input:q
  output:I0,I1,E,sum_0=I0+E,sum_1=I1+E
*/
void GenerateSet (Vec <ZZ_p>& I0,Vec <ZZ_p>& I1,Vec <ZZ_p>& sum_0,Vec <ZZ_p>& sum_1,Vec <ZZ>& E,ZZ& q)
{

  Vec<ZZ> I0_ZZ,I1_ZZ,sum_0_ZZ,sum_1_ZZ;
  double s1,s2;
  long q1,q2,q3,q4,q5;
  long i,j;

  conv (s1,q);
  s1 = s1 / 4;

  conv (s2,q);
  s2 = s2 / 8;
  q3 = ceil(-s2);

  if (q%8==0)
      s2 = s2 - 1;
  q4 = floor(s2);

  q5 = q4 - q3 +1;
  E.SetLength(q5);
  
  for(i=q3,j=0;i<=q4;i++,j++)//compute E=[-q/4,q/4)
  {
     E[j] = i;
  }

  q1 = round(s1);
  q2 = floor(s1);

  I0.SetLength(q1);
  I1.SetLength(q2);
  I0_ZZ.SetLength(q1);
  I1_ZZ.SetLength(q2);

  for(i=0;i<q1;i++)
  {
    I0[i] = i;
  }
  for(i=(-q2),j=0;i<0;i++,j++)
  {
    I1[j]=i;
  }


  conv(I0_ZZ,I0);
  conv(I1_ZZ,I1);

  ZZ min_0,max_0,min_1,max_1;
  min_0 = E[0];//I0_ZZ[0]=0;
  max_0 = I0_ZZ[q1-1] + E[q5-1];
  min_1 = I1_ZZ[0] + E[0];
  max_1 = I1_ZZ[q2-1] + E[q5-1];
 
  ZZ t_0,t_1;
  t_0 = max_0 - min_0 + 1;
  t_1 = max_1 - min_1 + 1;
  long t0,t1;
  conv(t0,t_0);
  conv(t1,t_1);
  sum_0_ZZ.SetLength(t0);
  sum_1_ZZ.SetLength(t1);
  

  for(i=0,t_0=min_0;t_0<=max_0;i++,t_0++)
      sum_0_ZZ[i] = t_0;
  conv(sum_0,sum_0_ZZ);//calculte I0+E
  
  for(i=0,t_1=min_1;t_1<=max_1;i++,t_1++)
      sum_1_ZZ[i] = t_1;
  conv(sum_1,sum_1_ZZ);//calculte I1+E


  //cout << "I0= "<< I0 << endl;
  //cout << "I1= "<< I1 << endl;
  //cout << "E= "<< E << endl;
  //cout << "I0+E= "<< sum_0 << endl;
  //cout << "I1+E= "<< sum_1 << endl;
}


/*function:<R>
  input:R,q
  output:cross_ZZ
*/
void CrossRounding (Vec <ZZ>& cross_ZZ,Vec <ZZ>& R,long& q)
{
  long n;
  n = R.length();
  long i;
  for (i=0;i<n;i++)
  {
    div(cross_ZZ[i],R[i]*4,q);
    cross_ZZ[i] = cross_ZZ[i] % 2;
  }
}
/*function:<R>
  input:R,q
  output:round_ZZ
*/
void ModularRounding (Vec <ZZ>& round_ZZ,Vec <ZZ>& R,long& q)
{
  long n;
  n = R.length();

  long i;
  double s;
  for (i = 0;i<n;i++)
  {
    conv (s,R[i]);
    s = (2*s) / q;
    round_ZZ[i] = round(s);
    round_ZZ[i] = round_ZZ[i] % 2;
  }
}
/*function:rec(w,b)
  input:w,cross_zz,
  output:rec_zz
*/
void Rec(Vec <ZZ>& coeff_R1_e,Vec <ZZ>& ZZ_vec_cross,Vec <ZZ>&ZZ_vec_rec,Vec <ZZ_p>& sum_0,Vec <ZZ_p>& sum_1)
{
  long i,j,flag;
  Vec<ZZ> sum_0_ZZ,sum_1_ZZ;
  conv(sum_0_ZZ,sum_0);
  conv(sum_1_ZZ,sum_1);

  for(i=0;i<ZZ_vec_rec.length();i++)
  {
    flag=0;
    if (ZZ_vec_cross[i]==0)
      {
         for(j=0;j<sum_0.length();j++)
         {
            if(sum_0_ZZ[j]==coeff_R1_e[i])
               { 
                 flag=1;
                 break;
               }    
         }
         if(flag==1)
           { ZZ_vec_rec[i]=0;
             continue;
           }
         else
           { ZZ_vec_rec[i]=1;
             continue;
           }         
       }  


    else
      {
         for(j=0;j<sum_1.length();j++)
         {
            if(sum_1_ZZ[j]==coeff_R1_e[i])
               { 
                 flag=1;
                 break;
               }    
         }
         if(flag==1)
           { ZZ_vec_rec[i]=0;
             continue;
           }
         else
           { ZZ_vec_rec[i]=1;
             continue;
           }         
       }  
  }

}
