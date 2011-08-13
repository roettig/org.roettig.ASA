package org.roettig.ASC.kernels;

import org.roettig.MLToolbox.base.instance.DualInstance;
import org.roettig.MLToolbox.kernels.KernelFunction;


public class AtchleyStringKernel extends KernelFunction<DualInstance<String>> 
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
			
			s += getC1(xi)*getC1(yi);
			s += getC2(xi)*getC2(yi);
			s += getC3(xi)*getC3(yi);
			s += getC4(xi)*getC4(yi);
			s += getC5(xi)*getC5(yi);
		}
		return s;
	}
	
	private static double[] CC1 = new double[26];
	private static double[] CC2 = new double[26];
	private static double[] CC3 = new double[26];
	private static double[] CC4 = new double[26];
	private static double[] CC5 = new double[26];

	static
	{
		CC1[0]=-0.591;CC1[17]=1.538;CC1[13]=0.945;CC1[3]=1.05;CC1[2]=-1.343;CC1[16]=0.931;CC1[4]=1.357;CC1[6]=-0.384;CC1[7]=0.336;CC1[8]=-1.239;CC1[11]=-1.019;CC1[10]=1.831;CC1[12]=-0.663;CC1[5]=-1.006;CC1[15]=0.189;CC1[18]=-0.228;CC1[19]=-0.032;CC1[22]=-0.595;CC1[24]=0.26;CC1[21]=-1.337;
		CC2[0]=-1.302;CC2[17]=-0.055;CC2[13]=0.828;CC2[3]=0.302;CC2[2]=0.465;CC2[16]=-0.179;CC2[4]=-1.453;CC2[6]=1.652;CC2[7]=-0.417;CC2[8]=-0.547;CC2[11]=-0.987;CC2[10]=-0.561;CC2[12]=-1.524;CC2[5]=-0.59;CC2[15]=2.081;CC2[18]=1.399;CC2[19]=0.326;CC2[22]=0.009;CC2[24]=0.83;CC2[21]=-0.279;
		CC3[0]=-0.733;CC3[17]=1.502;CC3[13]=1.299;CC3[3]=-3.656;CC3[2]=-0.862;CC3[16]=-3.005;CC3[4]=1.477;CC3[6]=1.33;CC3[7]=-1.673;CC3[8]=2.131;CC3[11]=-1.505;CC3[10]=0.533;CC3[12]=2.219;CC3[5]=1.891;CC3[15]=-1.628;CC3[18]=-4.76;CC3[19]=2.213;CC3[22]=0.672;CC3[24]=3.097;CC3[21]=-0.544;
		CC4[0]=1.57;CC4[17]=0.44;CC4[13]=-0.169;CC4[3]=-0.259;CC4[2]=-1.02;CC4[16]=-0.503;CC4[4]=0.113;CC4[6]=1.045;CC4[7]=-1.474;CC4[8]=0.393;CC4[11]=1.266;CC4[10]=-0.277;CC4[12]=-1.005;CC4[5]=-0.397;CC4[15]=0.421;CC4[18]=0.67;CC4[19]=0.908;CC4[22]=-2.128;CC4[24]=-0.838;CC4[21]=1.242;
		CC5[0]=-0.146;CC5[17]=2.897;CC5[13]=0.933;CC5[3]=-3.242;CC5[2]=-0.255;CC5[16]=-1.853;CC5[4]=-0.837;CC5[6]=2.064;CC5[7]=-0.078;CC5[8]=0.816;CC5[11]=-0.912;CC5[10]=1.648;CC5[12]=1.212;CC5[5]=0.412;CC5[15]=-1.392;CC5[18]=-2.647;CC5[19]=1.313;CC5[22]=-0.184;CC5[24]=1.512;CC5[21]=-1.262;
	}
	
	public static double getC1(char x)
	{
		if(x=='-'||x=='X')
			return 0.0;
		return CC1[ ((char) x)-65 ];
	}

	public static double getC2(char x)
	{
		if(x=='-'||x=='X')
			return 0.0;
		return CC2[ ((char) x)-65 ];
	}

	public static double getC3(char x)
	{
		if(x=='-'||x=='X')
			return 0.0;
		return CC3[ ((char) x)-65 ];
	}

	public static double getC4(char x)
	{
		if(x=='-'||x=='X')
			return 0.0;
		return CC4[ ((char) x)-65 ];
	}

	public static double getC5(char x)
	{
		if(x=='-'||x=='X')
			return 0.0;
		return CC5[ ((char) x)-65 ];
	}

}