package org.roettig.ASC.asc.interceptors;

import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.roettig.ASC.asc.ASAFlow;
import org.roettig.ASC.asc.ASAFlow.ASAFlowContext;
import org.roettig.MLToolbox.base.label.Label;
import org.roettig.SequenceTools.AlignedSequenceIdentity;
import org.roettig.SequenceTools.PairwiseAlignment;
import org.roettig.SequenceTools.SequenceSet;

public class DefaultSequenceFilterInterceptor implements ASAFlowInterceptor
{
	@Override
	public void intercept(ASAFlow flow, ASAFlowContext ctx)
	{
		SequenceSet keep = new SequenceSet();
		Map<String,Integer> occupancies = new HashMap<String,Integer>();
		Map<String,Double>  seqids = new HashMap<String,Double>();
		Map<String,String>  sequenceclasses = new HashMap<String,String>();

		for(Sequence seq : ctx.allseqs)
		{
			double pid = align(ctx.template,seq);
			ctx.seqdists.put(seq,pid);

			String orig_name = ctx.intern_seq_name_2_orig_seqname.get(seq.getName());

			if(pid>=ctx.jd.getPid_threshold())
			{
				flow.log("filter","Keeping sequence "+orig_name+"  "+pid);
				keep.add(seq);
				Label lab = ctx.intern_seq_name_2_class.get(seq.getName());
				String clsname = ctx.label_2_classname.get(lab);
				if(!occupancies.containsKey(clsname))
					occupancies.put(clsname,1);
				occupancies.put(clsname,occupancies.get(clsname)+1);
				seqids.put(orig_name, pid);
				sequenceclasses.put(orig_name,clsname);
			}
			else
				flow.log("filter","Filtered out sequence "+orig_name+"  "+pid);	    
		}
		ctx.jd.setOccupancies(occupancies);
		ctx.jd.setSeqIds(seqids);
		ctx.jd.setSequenceclasses(sequenceclasses);
		ctx.allseqs = keep;

		if(keep.size()==0)
		{
			flow.exit("critical error - No sequence passed the filtering step");
		}
	}
	
	protected PairwiseAlignment pwa = new PairwiseAlignment();
	
	public double align(Sequence x, Sequence y)
	{
		double pid = 0.0;
		try
		{
			pid = pwa.align(x, y, AlignedSequenceIdentity.getInstance());
		}
		catch( IllegalSymbolException e)
		{
			e.printStackTrace();
		}
		return pid;
	} 

}
