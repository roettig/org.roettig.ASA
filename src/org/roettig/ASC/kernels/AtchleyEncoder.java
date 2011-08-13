package org.roettig.ASC.kernels;

import org.roettig.MLToolbox.base.PrimalEncoder;

public class AtchleyEncoder implements PrimalEncoder<String>
{
	@Override
	public double[] encode(String t)
	{
		int L = t.length();
		int K = 5;
		double[] fts = new double[L*K];
		
		int z=0;
		for(int l=0;l<L;l++)
		{
			fts[z++] = AtchleyStringKernel.getC1(t.charAt(l));
			fts[z++] = AtchleyStringKernel.getC2(t.charAt(l));
			fts[z++] = AtchleyStringKernel.getC3(t.charAt(l));
			fts[z++] = AtchleyStringKernel.getC4(t.charAt(l));
			fts[z++] = AtchleyStringKernel.getC5(t.charAt(l));
		}
		return fts;
	}
	
}