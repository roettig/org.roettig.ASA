package org.roettig.ASC.asc;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.seq.Sequence;

import org.roettig.MLToolbox.base.Prediction;
import org.roettig.MLToolbox.base.instance.DualInstance;
import org.roettig.MLToolbox.base.instance.PrimalInstance;
import org.roettig.MLToolbox.base.label.FactorLabel;
import org.roettig.MLToolbox.base.label.Label;
import org.roettig.MLToolbox.kernels.KernelFunction;
import org.roettig.MLToolbox.kernels.RBFKernel;
import org.roettig.MLToolbox.model.CSVCModel;
import org.roettig.MLToolbox.model.SelectedModel;
import org.roettig.MLToolbox.validation.ModelValidation;
import org.roettig.PDBTools.ResidueLocator;
import org.roettig.SequenceTools.SeqTools;
import org.roettig.SequenceTools.SequenceSet;
import org.roettig.SequenceTools.exception.FileParseErrorException;
import org.roettig.maths.statistics.Statistics;

public class ASC extends ASAFlow
{
	public ASC(ASAJob _jd)
	{
		super(_jd);
	}
	
	@Override
	protected void readSequences()
	{
		int idx = 1;
		int C = ctx.jd.getFilenames().size();
		for(int c=1;c<=C;c++)
		{
			String filename = ctx.jd.getFilenames().get(c-1);
			
			Label lab = null;
			
			if(ctx.jd.getClassnames().get(c-1).toLowerCase().equals("unlabeled"))
				lab =  new FactorLabel(ctx.jd.getClassnames().get(c-1),0.0);
			else
				lab =  new FactorLabel(ctx.jd.getClassnames().get(c-1));
			
			ctx.labels.add(lab);
			
			log("readSequences"," "+filename+": mapping "+ctx.jd.getClassnames().get(c-1)+" --> "+lab);
			
			ctx.label_2_classname.put(lab,ctx.jd.getClassnames().get(c-1));

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
				String sid = String.format("%d|%d",c,idx);
				
				Sequence tmpseq = SeqTools.makeProteinSequence(sid,seq.seqString());

				ctx.intern_seq_name_2_orig_seqname.put(sid,seq.getName());
				ctx.intern_seq_name_2_class.put(sid,lab);
				ctx.seq2class.put(tmpseq,c);
				ctx.allseqs.add(tmpseq);
				idx++;
			}
		}

	}

	@Override
	protected void analyzeImportance()
	{
		/*
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

		int nC = ctx.labels.size();

		Logger ml2logger = Logger.getLogger("org.roettig.ML2.base.evaluation.ModelValidation");
		Level old_lvl = ml2logger.getLevel();
		ml2logger.setLevel(Level.OFF);
		
		
		for(int i=0;i<nC;i++)
		{
			for(int j=(i+1);j<nC;j++)
			{
				Vector<PositionImportance> imps = new Vector<PositionImportance>();

				Label pos = new FactorLabel("pos",1.0);
				Label neg = new FactorLabel("neg",2.0);
				List<DualInstance<String>> subset = ModelValidation.binarySubset(ctx.stringsamples, ctx.labels.get(i), ctx.labels.get(j), pos, neg);
				int K = ctx.ASRidx.size();
				double score[] = new double[K];
				for(int k=0;k<K;k++)
				{
					KernelFunction<DualInstance<String>> f = new IndexedSimilarityKernel( new StringHammingSimilarityDelegate() );
					Parameter<Integer> kIdx = new Parameter<Integer>(IndexedSimilarityKernel.IDX, new Integer[]{k});
					f.addParameter(kIdx);
					Model<DualInstance<String>> m = new CSVCmodel<DualInstance<String>>(f);
					Parameter<Double> C = new Parameter<Double>(CSVCmodel.C, new Double[]{100.0});
					m.addParameter(C);
					
					double qual = 0.0;
					try
					{
						qual = ModelValidation.CV(2, subset, m);
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

				log("analyzeImportance","calculating ASC residue importance for : "+ctx.labels.get(i)+" vs. "+ ctx.labels.get(j));
				for(PositionImportance imp : imps)
				{
					log("analyzeImportance",String.format(Locale.ENGLISH,"pos=%2d accu=%.2f",imp.pos,imp.imp));
					score[imp.pos-1] = imp.imp;
				}
				ctx.jd.setScore(i,j, score);
			}
		}
		ml2logger.setLevel(old_lvl);
		*/
	}

	@Override
	protected void validateModel()
	{
		if(ctx.jd.getNumberOfClasses()==1)
		{
			log("validateModel","no model validation since extract-only mode");
			return;
		}

		
		RBFKernel k = new RBFKernel();
		k.setGamma(0.001,0.01,0.1,1.0,10.0);
		
		
		CSVCModel<PrimalInstance> m = new CSVCModel<PrimalInstance>(k);
		m.setC(0.1,1.0,10.0,100.0,100.0);
		
		
		
		SelectedModel<PrimalInstance> bestModel = null;
		
		log("validateModel","### Doing Nested CV ####");

		int nInner = ctx.jd.getNInnerFolds();
		int nOuter = ctx.jd.getNOuterFolds();
		
		
		try
		{
			bestModel = ModelValidation.SimpleNestedCV( nOuter, nInner, ctx.samples, m);
		}
		catch(Exception e)
		{
			e.printStackTrace();
			exit("critical error - an error occured during Nested CV model checking");
		}
		
		ctx.model = bestModel.model;
		ctx.jd.setFmeasure(bestModel.qual);
		
		List<Label> yp = new ArrayList<Label>(); List<Label> yt = new ArrayList<Label>();
		Prediction.split(bestModel.predictions, yt, yp);
		
		Map<String,Double> precs = new HashMap<String,Double>();
		Map<String,Double> recs  = new HashMap<String,Double>();
		
		for(Label lab: ctx.labels)
		{
			List<Double> prec_rec = Statistics.calcPrecRec(lab, yp, yt);
			log("###",String.format(Locale.ENGLISH,"prec(%s)=%.5f rec(%s)=%.5f",lab.toString(),prec_rec.get(0),lab.toString(),prec_rec.get(1)));
			String classname = ctx.label_2_classname.get(lab);
			precs.put(classname,prec_rec.get(0));
			recs.put(classname,prec_rec.get(1));
		}

		ctx.jd.setPrecs(precs);
		ctx.jd.setRecs(recs);

		ctx.jd.setQualityName(bestModel.model.getQualityMeasure().getClass().getCanonicalName());
		ctx.jd.setQuality(bestModel.qual);
				
		log("validateModel","qual="+bestModel.qual);
		
		try
		{
			bestModel.model.store(ctx.outputdir+"/ASCmodel");
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
