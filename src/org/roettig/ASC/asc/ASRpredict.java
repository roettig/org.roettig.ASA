/**
 * 
 */
package org.roettig.ASC.asc;

import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;
import org.biojava.bio.seq.Sequence;
import org.roettig.ASC.kernels.WoldEncoder;
import org.roettig.MLToolbox.base.Prediction;
import org.roettig.MLToolbox.base.PrimalEncoder;
import org.roettig.MLToolbox.base.impl.DefaultInstanceContainer;
import org.roettig.MLToolbox.base.instance.InstanceContainer;
import org.roettig.MLToolbox.base.instance.PrimalInstance;
import org.roettig.MLToolbox.model.Model;
import org.roettig.PDBTools.SequenceAnchoredResidue;
import org.roettig.SequenceTools.MSA;
import org.roettig.SequenceTools.SeqTools;
import org.roettig.SequenceTools.SequenceSet;

/**
 * @author roettig
 *
 */
public class ASRpredict
{

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception
	{
		ASAJob   jd     = ASAJob.load("ASAdata");
		ASAModel ascmdl = ASAModel.load("ASRmodel");

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

		Model<PrimalInstance> m = ascmdl.getModel();
		SequenceSet out = new SequenceSet();
		PrimalEncoder<String> encoder = new WoldEncoder();
		
		for(Sequence seq: xsigs)
		{
			double[] fts = encoder.encode(seq.seqString());
			PrimalInstance pi = new PrimalInstance(null,fts);
			InstanceContainer<PrimalInstance> cnt = new DefaultInstanceContainer<PrimalInstance>();
			cnt.add(pi);
			
			List<Prediction> pred = m.predict(cnt);
			System.out.println(pred+" "+seq.getName());
			Sequence oseq = SeqTools.makeProteinSequence(String.format("%s#%s",pred,seq.getName()),seq.seqString());
			out.add(oseq);
		}

		out.store("xsigs.afa");
	}

}
