package org.roettig.ASC.asc.interceptors;

import org.biojava.bio.seq.Sequence;
import org.roettig.ASC.asc.ASAFlow;
import org.roettig.ASC.asc.ASAFlow.ASAFlowContext;
import org.roettig.SequenceTools.SequenceSet;

public class SignatureFilterInterceptor implements ASAFlowInterceptor
{

	@Override
	public void intercept(ASAFlow flow, ASAFlowContext ctx)
	{
		SequenceSet keep = new SequenceSet();
		for(Sequence seq : ctx.sigseqs)
		{
			String s = seq.seqString();
			int gaps = 0;
			int N    = s.length();
			for(int i=0;i<N;i++)
			{
				if(s.charAt(i)=='-')
					gaps++;
			}
			keep.add(seq);
			/*
			if(gaps<=N*0.15)
				keep.add(seq);
			else
			{
				flow.log("signatureFilter","dismissing signature sequence "+ctx.intern_seq_name_2_orig_seqname.get(seq.getName()));
			}
			*/
		}
		ctx.sigseqs = keep;
	}

}
