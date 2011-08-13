package org.roettig.ASC.asc;

import org.biojava.bio.seq.Sequence;
import org.roettig.ML2.base.DualInstance;
import org.roettig.ML2.base.Instance;
import org.roettig.ML2.base.Label;
import org.roettig.ML2.kernels.SimilarityMeasure;
import org.roettig.SequenceTools.AlignedSequenceIdentity;
import org.roettig.SequenceTools.GlobalSequenceIdentity;
import org.roettig.SequenceTools.SeqTools;


public class SequenceSimilarity implements SimilarityMeasure<DualInstance<String>>
{

	private static final long	serialVersionUID	= 2782899364793354337L;
	
	public static int GLOBAL  = 0;
	public static int ALIGNED = 1;
	private int type = GLOBAL;

	public SequenceSimilarity(int _type)
	{
		type = _type; 
	}

	@Override
	public double getSimilarity(DualInstance<String> x, DualInstance<String> y) throws Exception
	{
		Label lab = x.getLabel();
		Sequence seq1 = (Sequence) SeqTools.makeProteinSequence("a",x.getPayload());
		Sequence seq2 = (Sequence) SeqTools.makeProteinSequence("b",y.getPayload());
		double pid = 0.0;
		if(type==GLOBAL)         
			pid = GlobalSequenceIdentity.getInstance().calculate(seq1.seqString(), seq2.seqString());
		else
			pid = AlignedSequenceIdentity.getInstance().calculate(seq1.seqString(), seq2.seqString());
		if(lab.toDouble()==3.0)
			pid-=0.4;
		return pid;
	}

	@Override
	public double getDistance(DualInstance<String> x, DualInstance<String> y) throws Exception
	{
		Sequence seq1 = (Sequence) SeqTools.makeProteinSequence("a",x.getPayload());
		Sequence seq2 = (Sequence) SeqTools.makeProteinSequence("b",y.getPayload());
		double pid = 0.0;
		if(type==GLOBAL)         
			pid = GlobalSequenceIdentity.getInstance().calculate(seq1.seqString(), seq2.seqString());
		else
			pid = AlignedSequenceIdentity.getInstance().calculate(seq1.seqString(), seq2.seqString());
		return 1.0-pid;
	}

}
