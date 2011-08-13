package org.roettig.ASC.kernels;

import org.roettig.MLToolbox.base.PrimalEncoder;

public class WoldEncoder implements PrimalEncoder<String>
{
	@Override
	public double[] encode(String t)
	{
		int L = t.length();
		int K = 3;
		double[] fts = new double[L*K];
		
		int z=0;
		for(int l=0;l<L;l++)
		{
			fts[z++] = WoldStringKernel.getZ1(t.charAt(l));
			fts[z++] = WoldStringKernel.getZ2(t.charAt(l));
			fts[z++] = WoldStringKernel.getZ3(t.charAt(l));
		}
		return fts;
	}
	
}
