/**
 * 
 */
package org.roettig.ASC.asc;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import org.biojava.bio.seq.Sequence;
import org.roettig.MLToolbox.base.instance.PrimalInstance;
import org.roettig.MLToolbox.base.label.Label;
import org.roettig.MLToolbox.base.label.NumericLabel;
import org.roettig.MLToolbox.kernels.RBFKernel;
import org.roettig.MLToolbox.model.NuSVRModel;
import org.roettig.MLToolbox.model.SelectedModel;
import org.roettig.MLToolbox.validation.ModelValidation;
import org.roettig.PDBTools.ResidueLocator;
import org.roettig.SequenceTools.SeqTools;
import org.roettig.SequenceTools.SequenceSet;
import org.roettig.SequenceTools.exception.FileParseErrorException;

/**
 * @author roettig
 *
 */
public class ASR extends ASAFlow
{
	public ASR(ASAJob _jd)
	{
		super(_jd);
	}

	@Override
	protected void readSequences()
	{
		
		String filename = ctx.jd.getFilenames().get(0);
		int idx = 1;
		SequenceSet seqs = null;
		try
		{
			seqs = SequenceSet.readFromFile(filename);
		} 
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
			exit("critical error - could not find sequence file "+filename);
		} 
		catch (FileParseErrorException e)
		{
			e.printStackTrace();
			exit("critical error - could not parse sequence file "+filename);
		}

		for(Sequence seq: seqs)
		{
			String toks[] = seq.getName().split("#");
			double y = Double.parseDouble(toks[toks.length-1]);
			Label lab = new NumericLabel(y);
			String sid = String.format(Locale.ENGLISH,"%.6f|%d",y,idx);
			Sequence tmpseq = SeqTools.makeProteinSequence(sid,seq.seqString());
			ctx.intern_seq_name_2_orig_seqname.put(sid,seq.getName());
			ctx.intern_seq_name_2_class.put(sid,lab);
			ctx.allseqs.add(tmpseq);
			idx++;
		}
	}
	
	

	@Override
	protected void analyzeImportance()
	{
	}

	@Override
	protected void validateModel()
	{
		RBFKernel k = new RBFKernel();
		double g = RBFKernel.estimateGamma(ctx.samples);
		k.setGamma(g/8.0,g/4.0,g/2.0,g,2*g,4*g,8*g);
		
		NuSVRModel<PrimalInstance> m = new NuSVRModel<PrimalInstance>(k);
		m.setC(0.1,1.0,10.0,100.0,100.0);
		m.setNU(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9);
		
		
		SelectedModel<PrimalInstance> bestModel = null;
		log("validateModel","### Doing Nested CV ####");

		int nInner = ctx.jd.getNInnerFolds();
		int nOuter = ctx.jd.getNOuterFolds();
		

		try
		{
			bestModel = ModelValidation.SimpleNestedCV(2,2,ctx.samples,m);
		}
		catch(Exception e)
		{
			e.printStackTrace();
			exit("critical error - an error occured during Nested CV model checking");
		}
		ctx.model = bestModel.model;
		log("validateModel","qual="+bestModel.qual);			
		try
		{
			bestModel.model.store(ctx.outputdir+"/ASRmodel");
		} 
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args)
	{
		ASAJob jd = new ASAJob();
		jd.setPdbfile("/tmp/1AMU.pdb");
		jd.setOutputDir("/tmp");
		ResidueLocator resloc = new ResidueLocator("A","PHE",566);
		List<ResidueLocator> reslocs = new ArrayList<ResidueLocator>();
		reslocs.add(resloc);
		jd.setReslocs(reslocs);
		jd.setDistance(8.0);
		List<String> filenames = new ArrayList<String>();
		//filenames.add("/tmp/3.fas");
		filenames.add("/tmp/reg.fa");
		//filenames.add("/tmp/tyr.fa");
		//filenames.add("/tmp/asp.fa");
		jd.setFilenames(filenames);
		List<String> classnames = new ArrayList<String>();
		//classnames.add("ala");
		//classnames.add("tyr");
		//classnames.add("asp");
		jd.setClassnames(classnames);
		ASAFlow flow = new ASR(jd);
		flow.run();		
	}
	
}
