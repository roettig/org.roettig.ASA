package org.roettig.ASC.kernels;

import org.roettig.MLToolbox.base.instance.DualInstance;
import org.roettig.MLToolbox.kernels.KernelFunction;
import org.roettig.SequenceTools.SeqTools;


public class BLOSUMStringKernel extends KernelFunction<DualInstance<String>>
{
	@Override
	public double compute(DualInstance<String> x, DualInstance<String> y)
			throws Exception
	{
		String sx = x.getPayload();
		String sy = y.getPayload();
		
		if(sx.length()!=sy.length())
        	throw new Exception("error evaluating kernel function "+this.getClass().getCanonicalName()+" : supplied strings have different lengths");
        
		
		double sim = 0.0;
		
		for(int i=0;i<sx.length();i++)
		{
 			sim += SeqTools.getNormalizedBLOSUM62Score(sx.charAt(i),sy.charAt(i));
		}
		
		return sim; 
	}
}
