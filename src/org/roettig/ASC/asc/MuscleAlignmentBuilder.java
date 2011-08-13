/**
 * 
 */
package org.roettig.ASC.asc;

import org.roettig.SequenceTools.MSA;
import org.roettig.SequenceTools.SequenceSet;

/**
 * @author roettig
 *
 */
public class MuscleAlignmentBuilder implements AlignmentBuilder
{
	@Override
	public MSA build(SequenceSet seqs) throws Exception
	{
		MSA ret = MSA.createMuscleMSA(seqs);
		return ret;
	}

}
