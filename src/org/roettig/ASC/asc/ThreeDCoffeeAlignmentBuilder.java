/**
 * 
 */
package org.roettig.ASC.asc;

import org.roettig.SequenceTools.MSA;
import org.roettig.SequenceTools.SequenceSet;
import org.roettig.SequenceTools.ThreeDCoffeeAlignment;

/**
 * @author roettig
 *
 */
public class ThreeDCoffeeAlignmentBuilder implements AlignmentBuilder
{
	private String jobid = null;

	ThreeDCoffeeAlignmentBuilder()
	{
	}

	ThreeDCoffeeAlignmentBuilder(String _jobid)
	{
		jobid = _jobid;
	}

	@Override
	public MSA build(SequenceSet seqs) throws Exception
	{
		ThreeDCoffeeAlignment ret =  new ThreeDCoffeeAlignment(seqs);
		if(jobid==null)
			return ret.align();
		else
			return ret.fetchMSA(jobid);
	}

}
