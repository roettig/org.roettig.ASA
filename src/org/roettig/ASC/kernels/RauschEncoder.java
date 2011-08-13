package org.roettig.ASC.kernels;

import org.roettig.MLToolbox.base.PrimalEncoder;

public class RauschEncoder implements PrimalEncoder<String>
{
	@Override
	public double[] encode(String t)
	{
		int L = t.length();
		int K = 12;
		double[] fts = new double[L*K];
		
		int z=0;
		for(int l=0;l<L;l++)
		{
			for(int f=1;f<=12;f++)
				fts[z++] = RauschStringKernel.getDescriptor(f, t.charAt(l));
		}
		return fts;
	}
	
}
