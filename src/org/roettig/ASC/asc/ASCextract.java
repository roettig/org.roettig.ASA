/**
 * 
 */
package org.roettig.ASC.asc;

import java.util.Map;
import java.util.TreeMap;

import org.roettig.PDBTools.SequenceAnchoredResidue;
import org.roettig.SequenceTools.MSA;
import org.roettig.SequenceTools.SeqTools;
import org.roettig.SequenceTools.SequenceSet;

/**
 * @author roettig
 *
 */
public class ASCextract
{
	public static void main(String[] args) throws Exception
	{
		ASAJob   jd     = ASAJob.load("ASCresult");
		ASAModel ascmdl = ASAModel.load("ASCmodel");

		SequenceSet seqs = SequenceSet.readFromFile(args[0]);
		seqs.add(SeqTools.makeProteinSequence("pdb",jd.getTemplateSequence()));
		MSA mmsa = MSA.align(ascmdl.getMSA(), seqs);
		mmsa.store("test.afa");

		Map<Integer,SequenceAnchoredResidue> ASCidx = new TreeMap<Integer,SequenceAnchoredResidue>();

		for(SequenceAnchoredResidue aa: jd.getASCresidues())
		{
			ASCidx.put(aa.seqIdx, aa);
		}

		SequenceSet xsigs = mmsa.getSubSequences(ASCidx.keySet(), "pdb");
		xsigs.store("xsigs.afa");
	}
}
