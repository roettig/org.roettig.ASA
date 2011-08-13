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
public class MockAlignmentBuilder implements AlignmentBuilder
{
	private MSA msa;

	public MockAlignmentBuilder(MSA _msa)
	{
		msa = _msa;
	}

	@Override
	public MSA build(SequenceSet seqs) throws Exception
	{
		return MSA.align(msa, seqs);
	}

}
