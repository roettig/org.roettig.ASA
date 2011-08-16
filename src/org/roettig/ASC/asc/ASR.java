package org.roettig.ASC.asc;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.seq.Sequence;
import org.roettig.ASC.kernels.BLOSUMStringKernel;
import org.roettig.ASC.kernels.IndexedStringKernel;
import org.roettig.MLToolbox.base.Prediction;
import org.roettig.MLToolbox.base.impl.DefaultInstanceContainer;
import org.roettig.MLToolbox.base.instance.DualInstance;
import org.roettig.MLToolbox.base.instance.InstanceContainer;
import org.roettig.MLToolbox.base.instance.PrimalInstance;
import org.roettig.MLToolbox.base.label.FactorLabel;
import org.roettig.MLToolbox.base.label.Label;
import org.roettig.MLToolbox.base.label.NumericLabel;
import org.roettig.MLToolbox.kernels.RBFKernel;
import org.roettig.MLToolbox.model.CSVCModel;
import org.roettig.MLToolbox.model.NuSVRModel;
import org.roettig.MLToolbox.model.SelectedModel;
import org.roettig.MLToolbox.validation.GMeasure;
import org.roettig.MLToolbox.validation.ModelValidation;
import org.roettig.PDBTools.ResidueLocator;
import org.roettig.SequenceTools.SeqTools;
import org.roettig.SequenceTools.SequenceSet;
import org.roettig.SequenceTools.exception.FileParseErrorException;

/**
 * The ASR class is used to conduct active site regression analysis.
 * 
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
			String toks[] = seq.getName().split("\\|");
			
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
		final class PositionImportance
		{
			public PositionImportance(int _pos, double _imp)
			{
				pos = _pos; imp = _imp;  
			}

			public boolean equals(Object o)
			{
				if (!(o instanceof PositionImportance))
					return false;

				PositionImportance p = (PositionImportance) o;
				return pos==p.pos;
			}

			public int hashCode()
			{
				return pos;
			}


			public int pos;
			public double imp;
		}

		

		Logger ml2logger = Logger.getLogger("org.roettig.MLToolbox.validation.ModelValidation");
		Level old_lvl = ml2logger.getLevel();
		ml2logger.setLevel(Level.OFF);


		int K = ctx.ASRidx.size();

		double score[] = new double[K];
		
		Vector<PositionImportance> imps = new Vector<PositionImportance>();
		
		for(int k=0;k<K;k++)
		{
			IndexedStringKernel f = new IndexedStringKernel( new BLOSUMStringKernel() );
			f.addFixedIndex(k);


			NuSVRModel<DualInstance<String>> m = new NuSVRModel<DualInstance<String>>(f);
			
			m.setNU(0.3);
			m.setC(100.0);

			double qual = 0.0;
			try
			{
				qual = ModelValidation.CV(2, ctx.stringsamples, m);
			} 
			catch (Exception e)
			{
				e.printStackTrace();
			}
			
			PositionImportance imp = new PositionImportance((int)k+1,qual);
			imps.add(imp);
		}

		Collections.sort(imps, new Comparator< PositionImportance >() {
			public int compare(PositionImportance e1, PositionImportance e2) 
			{
				Double accu1 = e1.imp;Double accu2 = e2.imp; 
				return accu1.compareTo(accu2);
			}});

		log("analyzeImportance","calculating ASC residue importance");
		for(PositionImportance imp : imps)
		{
			log("analyzeImportance",String.format(Locale.ENGLISH,"pos=%2d accu=%.2f",imp.pos,imp.imp));
			score[imp.pos-1] = imp.imp;
		}
		ctx.jd.setScore(0, 1, score);
		
		ml2logger.setLevel(old_lvl);
		
	}

	@Override
	protected void validateModel()
	{
		RBFKernel k = new RBFKernel();
		double g = RBFKernel.estimateGamma(ctx.samples);
		k.setGamma(g/8.0,g/4.0,g/2.0,g,2*g,4*g,8*g);
		
		NuSVRModel<PrimalInstance> m = new NuSVRModel<PrimalInstance>(k);
		m.setC(0.1,1.0,10.0,100.0,1000.0);
		m.setNU(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9);
		
		SelectedModel<PrimalInstance> bestModel = null;
		log("validateModel","### Doing Nested CV ####");

		int nInner = ctx.jd.getNInnerFolds();
		int nOuter = ctx.jd.getNOuterFolds();
		

		try
		{
			bestModel = ModelValidation.SimpleNestedCV( nInner, nOuter, ctx.samples, m);
		}
		catch(Exception e)
		{
			e.printStackTrace();
			exit("critical error - an error occured during Nested CV model checking");
		}
		
		ctx.model = bestModel.model;
		ctx.jd.setQuality(bestModel.qual);
		ctx.jd.setQualityName(bestModel.model.getQualityMeasure().getClass().getCanonicalName());
		
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

	@Override
	protected void buildFinalModel()
	{	
	}
	
}
