package org.roettig.ASC.asc.interceptors;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.seq.Sequence;
import org.roettig.ASC.asc.ASAFlow;
import org.roettig.ASC.asc.AlignmentBuilder;
import org.roettig.ASC.asc.ASAFlow.ASAFlowContext;
import org.roettig.SequenceTools.HMM;

public class AlignmentInterceptor implements ASAFlowInterceptor 
{
	AlignmentBuilder alignmentbuilder;
	
	public AlignmentInterceptor(AlignmentBuilder _builder)
	{
		alignmentbuilder = _builder;
	}
	
	@Override
	public void intercept(ASAFlow flow, ASAFlowContext ctx)
	{
		Map<Double,Integer> sids = new HashMap<Double,Integer>();

		for(Sequence seq : ctx.allseqs)
		{
			double pid   = ctx.seqdists.get(seq);
			Double range = seqRange(pid);
			Integer seq_class = ctx.seq2class.get(seq);
			if(!sids.containsKey(range))
			{
				sids.put(range, 1);
				ctx.coreseqs.add(seq);
				flow.log("buildAlignment","adding sequence "+seq.getName()+" into core set "+seq_class+" dist="+pid);
			}
			else
			{
				if(sids.get(range)>=2)
				{
					continue;
				}
				else
				{
					sids.put(range, sids.get(range)+1);
					ctx.coreseqs.add(seq);
					flow.log("buildAlignment","adding sequence "+seq.getName()+" into core set "+seq_class+" dist="+pid);
				}
			}

		}
		
		ctx.allseqs.add(ctx.template);
		ctx.coreseqs.add(ctx.template);

		ctx.coreseqs.store(ctx.outputdir+"/core.fa");

		try
		{
			ctx.coremsa = alignmentbuilder.build(ctx.coreseqs);
		} 
		catch (Exception e)
		{
			e.printStackTrace();
			flow.exit("critical error: building multiple sequence alignment failed - MSA step");
		}

		flow.log("buildAlignment","built core MSA");
		ctx.coremsa.store(ctx.outputdir+"/core.afa");

		try
		{
			ctx.corehmm = new HMM(ctx.coremsa);
		} 
		catch (Exception e)
		{
			e.printStackTrace();
			flow.exit("critical error: building multiple sequence alignment failed - HMM build step");
		}
		flow.log("buildAlignment","built core HMM");

		try
		{
			ctx.corehmm.save(ctx.outputdir+"/raus.hmm");
		} 
		catch (IOException e)
		{    
			e.printStackTrace();
		}

		try
		{
			ctx.allmsa = ctx.corehmm.align(ctx.allseqs);
		} 
		catch (Exception e)
		{
			e.printStackTrace();
			flow.exit("critical error: building multiple sequence alignment failed - HMM align step");
		}

		ctx.allmsa.store(ctx.outputdir+"/all.afa");
	}
	
	public static Double seqRange(double pid)
	{
		if(pid>0.9)
			return 0.9;
		if(pid>0.8)
			return 0.8;
		if(pid>0.7)
			return 0.7;
		if(pid>0.6)
			return 0.6;
		if(pid>0.5)
			return 0.5;
		if(pid>0.4)
			return 0.4;
		if(pid>0.3)
			return 0.3;
		if(pid>0.2)
			return 0.2;
		return 0.1;
	}
}
