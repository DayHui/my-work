#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZ_p.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include "api.h"

using namespace std;
using namespace NTL;

unsigned long step = sizeof(unsigned long)*8;
int main()
{

		RR sigma;
                cin>>sigma;

		long CDT_length;
		unsigned long* CDT;
		long CDT_inv_min[256];
		long CDT_inv_max[256];
		long tau=6;


		RR f, ff;

		mul(f, sigma, sigma);
		f = 2.0*f;

		conv(CDT_length, tau * sigma + 1.0);
		CDT = new unsigned long[CDT_length *2];


		// Compute CDT
		NTL::RR t, z;
		t=0;
		for(long i=1; i<CDT_length; i++)
		{
			conv(z, i-1);
			mul(z, z, -z);
			div(z, z, f);
			exp(z, z);
			if(i==1)
				z/=2.0;
			add(t, t, z);
		}

		NTL::RR y, tt;

		power2(tt, step);
		y=0;

		NTL::ZZ temp;

		for(long i=1; i<CDT_length; i++)
		{
			conv(z, i-1);
			mul(z, z, -z);
			div(z, z, f);
			exp(z, z);
			if(i==1)
				z/=2.0;

			div(z, z, t);
			add(y, y, z);

			z = y;

			for(long j=0; j<2; j++) 
			{
				mul(z, z, tt);
				FloorToZZ(temp, z);
				conv(CDT[i+j*CDT_length], z);
				sub(z, z, CDT[i+j*CDT_length]);
			}
		}

		for(long j=0; j<2; j++)
			CDT[j*CDT_length] = 0;


		long min=0, max = 0;
		unsigned long val;
		unsigned long mask = 0xFF << (step-8);


		for(long i=0;i<256;i++)
		{
			val = ((unsigned long) i) << (step-8);

			while(CDT[min+1]<val)
				min++;

			while((max+1 < CDT_length) && ((CDT[max] & mask) <= val))
				max++;

			CDT_inv_min[i]=min;
			CDT_inv_max[i]=max;
		}
		cout<<"CDT=";
		for(long i=0;i<CDT_length *2;i++)
			cout<<CDT[i]<<" ";
                cout<<endl;
		cout<<"CDT_inv_min=";
		for(long i=0;i<256;i++)
			cout<<CDT_inv_min[i]<<" ";
                cout<<endl;
		cout<<"CDT_inv_max=";
		for(long i=0;i<256;i++)
			cout<<CDT_inv_max[i]<<" ";
		cout<<endl;
}

