/**
 * 
 */
package org.roettig.ASC.tools;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.seq.Sequence;
import org.roettig.ASC.kernels.WoldStringKernel;
import org.roettig.MLToolbox.base.instance.PrimalInstance;


/**
 * @author roettig
 *
 */
public class ASCTools
{
	public static PrimalInstance encodePrimalWold(Sequence s)
	{

		int W = s.seqString().length();
		double x[] = new double[W*3];
		int idx = 0;
		for(int i=0;i<W;i++)
		{
			x[idx] = WoldStringKernel.getZ1(s.seqString().charAt(i));
			idx++;
			x[idx] = WoldStringKernel.getZ2(s.seqString().charAt(i));
			idx++;
			x[idx] = WoldStringKernel.getZ3(s.seqString().charAt(i));
			idx++;
		}
		PrimalInstance pi = new PrimalInstance(null,x);

		return pi;
	}
}
