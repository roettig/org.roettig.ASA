package org.roettig.ASC.kernels;

import java.util.Arrays;
import java.util.List;

import org.roettig.MLToolbox.base.instance.DualInstance;
import org.roettig.MLToolbox.kernels.KernelFunction;
import org.roettig.maths.statistics.Statistics;


public class RauschStringKernel extends KernelFunction<DualInstance<String>>
{

	@Override
	public double compute(DualInstance<String> x, DualInstance<String> y) throws Exception
	{
		String sx = x.getPayload();
		String sy = y.getPayload();
		
		if(sx.length()!=sy.length())
        	throw new Exception("error evaluating kernel function "+this.getClass().getCanonicalName()+" : supplied strings have different lengths");
        
		
		double s = 0.0;
		for(int i=0;i<sx.length();i++)
		{
			char xi = sx.charAt(i);
			char yi = sy.charAt(i);
			
			for(int z=1;z<=12;z++)
			{
				s += getDescriptor(z, xi) * getDescriptor(z, yi);  
			}
			
		}
		return s;
	}

	private static Double D1[] = new Double[26];
	private static Double D2[] = new Double[26];
	private static Double D3[] = new Double[26];
	private static Double D4[] = new Double[26];
	private static Double D5[] = new Double[26];
	private static Double D6[] = new Double[26];
	private static Double D7[] = new Double[26];
	private static Double D8[] = new Double[26];
	private static Double D9[] = new Double[26];
	private static Double D10[] = new Double[26];
	private static Double D11[] = new Double[26];
	private static Double D12[] = new Double[26];

	static
	{
		//1 aa-alpha-helix.aaindex
		D1[0]=1.420;D1[1]=0.000;D1[2]=0.700;D1[3]=1.010;D1[4]=1.510;D1[5]=1.130;D1[6]=0.570;D1[7]=1.000;D1[8]=1.080;D1[9]=0.000;D1[10]=1.160;D1[11]=1.210;D1[12]=1.450;D1[13]=0.670;D1[14]=0.000;D1[15]=0.570;D1[16]=1.110;D1[17]=0.980;D1[18]=0.770;D1[19]=0.830;D1[20]=0.000;D1[21]=1.060;D1[22]=1.080;D1[23]=0.000;D1[24]=0.690;D1[25]=0.000;
		//2 aa-beta-sheet.aaindex
		D2[0]=0.830;D2[1]=0.000;D2[2]=1.190;D2[3]=0.540;D2[4]=0.370;D2[5]=1.380;D2[6]=0.750;D2[7]=0.870;D2[8]=1.600;D2[9]=0.000;D2[10]=0.740;D2[11]=1.300;D2[12]=1.050;D2[13]=0.890;D2[14]=0.000;D2[15]=0.550;D2[16]=1.100;D2[17]=0.930;D2[18]=0.750;D2[19]=1.190;D2[20]=0.000;D2[21]=1.700;D2[22]=1.370;D2[23]=0.000;D2[24]=1.470;D2[25]=0.000;
		//3 aa-beta-turn.aaindex
		D3[0]=0.740;D3[1]=0.000;D3[2]=0.960;D3[3]=1.520;D3[4]=0.950;D3[5]=0.660;D3[6]=1.560;D3[7]=0.950;D3[8]=0.470;D3[9]=0.000;D3[10]=1.190;D3[11]=0.500;D3[12]=0.600;D3[13]=1.460;D3[14]=0.000;D3[15]=1.560;D3[16]=0.960;D3[17]=1.010;D3[18]=1.430;D3[19]=0.980;D3[20]=0.000;D3[21]=0.590;D3[22]=0.600;D3[23]=0.000;D3[24]=1.140;D3[25]=0.000;
		//4 aa-hydrogenbond.aaindex
		D4[0]=0.000;D4[1]=0.000;D4[2]=0.000;D4[3]=1.000;D4[4]=1.000;D4[5]=0.000;D4[6]=0.000;D4[7]=1.000;D4[8]=0.000;D4[9]=0.000;D4[10]=2.000;D4[11]=0.000;D4[12]=0.000;D4[13]=2.000;D4[14]=0.000;D4[15]=0.000;D4[16]=2.000;D4[17]=4.000;D4[18]=1.000;D4[19]=1.000;D4[20]=0.000;D4[21]=0.000;D4[22]=1.000;D4[23]=0.000;D4[24]=1.000;D4[25]=0.000;
		//5 aa-hydrophobicity-neu1.aaindex
		D5[0]=0.060;D5[1]=0.000;D5[2]=-0.560;D5[3]=0.970;D5[4]=0.850;D5[5]=-0.990;D5[6]=0.320;D5[7]=0.150;D5[8]=-1.000;D5[9]=0.000;D5[10]=1.000;D5[11]=-0.830;D5[12]=-0.680;D5[13]=0.700;D5[14]=0.000;D5[15]=0.450;D5[16]=0.710;D5[17]=0.800;D5[18]=0.480;D5[19]=0.380;D5[20]=0.000;D5[21]=-0.750;D5[22]=-0.570;D5[23]=0.000;D5[24]=-0.350;D5[25]=0.000;
		//6 aa-hydrophobicity-neu2.aaindex
		D6[0]=-0.250;D6[1]=0.000;D6[2]=-0.400;D6[3]=-0.080;D6[4]=-0.100;D6[5]=0.180;D6[6]=-0.320;D6[7]=-0.030;D6[8]=-0.030;D6[9]=0.000;D6[10]=0.320;D6[11]=0.050;D6[12]=-0.010;D6[13]=-0.060;D6[14]=0.000;D6[15]=0.230;D6[16]=-0.020;D6[17]=0.190;D6[18]=-0.150;D6[19]=-0.100;D6[20]=0.000;D6[21]=-0.190;D6[22]=0.310;D6[23]=0.000;D6[24]=0.400;D6[25]=0.000;
		//7 aa-hydrophobicity-neu3.aaindex
		D7[0]=0.250;D7[1]=0.000;D7[2]=-0.140;D7[3]=0.080;D7[4]=-0.050;D7[5]=0.150;D7[6]=0.280;D7[7]=-0.100;D7[8]=0.100;D7[9]=0.000;D7[10]=0.110;D7[11]=0.010;D7[12]=0.040;D7[13]=0.170;D7[14]=0.000;D7[15]=0.410;D7[16]=0.120;D7[17]=-0.410;D7[18]=0.230;D7[19]=0.290;D7[20]=0.000;D7[21]=0.030;D7[22]=0.340;D7[23]=0.000;D7[24]=-0.020;D7[25]=0.000;
		//8 aa-isoelectric.aaindex
		D8[0]=6.000;D8[1]=0.000;D8[2]=5.050;D8[3]=2.770;D8[4]=3.220;D8[5]=5.480;D8[6]=5.970;D8[7]=7.590;D8[8]=6.020;D8[9]=0.000;D8[10]=9.740;D8[11]=5.980;D8[12]=5.740;D8[13]=5.410;D8[14]=0.000;D8[15]=6.300;D8[16]=5.650;D8[17]=10.760;D8[18]=5.680;D8[19]=5.660;D8[20]=0.000;D8[21]=5.960;D8[22]=5.890;D8[23]=0.000;D8[24]=5.660;D8[25]=0.000;
		//9 aa-polar-grantham.aaindex
		D9[0]=8.100;D9[1]=0.000;D9[2]=5.500;D9[3]=13.000;D9[4]=12.300;D9[5]=5.200;D9[6]=9.000;D9[7]=10.400;D9[8]=5.200;D9[9]=0.000;D9[10]=11.300;D9[11]=4.900;D9[12]=5.700;D9[13]=11.600;D9[14]=0.000;D9[15]=8.000;D9[16]=10.500;D9[17]=10.500;D9[18]=9.200;D9[19]=8.600;D9[20]=0.000;D9[21]=5.900;D9[22]=5.400;D9[23]=0.000;D9[24]=6.200;D9[25]=0.000;
		//10 aa-polar-radzicka.aaindex
		D10[0]=-0.060;D10[1]=0.000;D10[2]=1.360;D10[3]=-0.800;D10[4]=-0.770;D10[5]=1.270;D10[6]=-0.410;D10[7]=0.490;D10[8]=1.310;D10[9]=0.000;D10[10]=-1.180;D10[11]=1.210;D10[12]=1.270;D10[13]=-0.480;D10[14]=0.000;D10[15]=1.100;D10[16]=-0.730;D10[17]=-0.840;D10[18]=-0.500;D10[19]=-0.270;D10[20]=0.000;D10[21]=1.090;D10[22]=0.880;D10[23]=0.000;D10[24]=0.330;D10[25]=0.000;
		//11 aa-polar-zimmerman.aaindex
		D11[0]=0.000;D11[1]=0.000;D11[2]=1.480;D11[3]=49.700;D11[4]=49.900;D11[5]=0.350;D11[6]=0.000;D11[7]=51.600;D11[8]=0.130;D11[9]=0.000;D11[10]=49.500;D11[11]=0.130;D11[12]=1.430;D11[13]=3.380;D11[14]=0.000;D11[15]=1.580;D11[16]=3.530;D11[17]=52.000;D11[18]=1.670;D11[19]=1.660;D11[20]=0.000;D11[21]=0.130;D11[22]=2.100;D11[23]=0.000;D11[24]=1.610;D11[25]=0.000;
		//12 aa-volume.aaindex
		D12[0]=90.000;D12[1]=0.000;D12[2]=103.300;D12[3]=117.300;D12[4]=142.200;D12[5]=191.900;D12[6]=64.900;D12[7]=160.000;D12[8]=163.900;D12[9]=0.000;D12[10]=167.300;D12[11]=164.000;D12[12]=167.000;D12[13]=124.700;D12[14]=0.000;D12[15]=122.900;D12[16]=149.400;D12[17]=194.000;D12[18]=95.400;D12[19]=121.500;D12[20]=0.000;D12[21]=139.000;D12[22]=228.200;D12[23]=0.000;D12[24]=197.000;D12[25]=0.000;
	}
	
	public RauschStringKernel()
	{
		standardize(D1);
		standardize(D2);
		standardize(D3);
		standardize(D4);
		standardize(D5);
		standardize(D6);
		standardize(D7);
		standardize(D8);
		standardize(D9);
		standardize(D10);
		standardize(D11);
		standardize(D12);
	}

	public static double getDescriptor(int idx, char x)
	{
		if(x=='-'||x=='X')
			return 0.0;
		switch(idx)
		{
		case 1:
			return D1[ ((char) x)-65 ];
		case 2:
			return D2[ ((char) x)-65 ];
		case 3:
			return D3[ ((char) x)-65 ];
		case 4:
			return D4[ ((char) x)-65 ];
		case 5:
			return D5[ ((char) x)-65 ];
		case 6:
			return D6[ ((char) x)-65 ];
		case 7:
			return D7[ ((char) x)-65 ];
		case 8:
			return D8[ ((char) x)-65 ];
		case 9:
			return D9[ ((char) x)-65 ];
		case 10:
			return D10[ ((char) x)-65 ];
		case 11:
			return D11[ ((char) x)-65 ];
		case 12:
			return D12[ ((char) x)-65 ];   	     
		}
		return 0;
	}
	
	public static void standardize(Double D[])
	{
		double m=0.0; double s =1.0;
		List<Double> l = Arrays.asList(D);
		m = Statistics.mean(l);
		s = Statistics.std(l);
		for(int i=0;i<D.length;i++)
		{
			D[i]-=m;
			if(s>0)
				D[i]/=s;
		}
	}
}