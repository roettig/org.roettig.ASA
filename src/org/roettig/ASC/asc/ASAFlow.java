package org.roettig.ASC.asc;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.PDBFileReader;
import org.roettig.ASC.asc.interceptors.ASAFlowInterceptor;
import org.roettig.ASC.asc.interceptors.AlignmentInterceptor;
import org.roettig.ASC.asc.interceptors.SequenceFilterInterceptor;
import org.roettig.ASC.asc.interceptors.SignatureFilterInterceptor;
import org.roettig.ASC.kernels.WoldEncoder;
import org.roettig.MLToolbox.base.Prediction;
import org.roettig.MLToolbox.base.PrimalEncoder;
import org.roettig.MLToolbox.base.impl.DefaultInstanceContainer;
import org.roettig.MLToolbox.base.instance.DualInstance;
import org.roettig.MLToolbox.base.instance.InstanceContainer;
import org.roettig.MLToolbox.base.instance.PrimalInstance;
import org.roettig.MLToolbox.base.label.FactorLabel;
import org.roettig.MLToolbox.base.label.Label;
import org.roettig.MLToolbox.kernels.RBFKernel;
import org.roettig.MLToolbox.model.CSVCModel;
import org.roettig.MLToolbox.model.Model;
import org.roettig.MLToolbox.model.SelectedModel;
import org.roettig.MLToolbox.validation.ModelValidation;
import org.roettig.PDBTools.PDBTools;
import org.roettig.PDBTools.ResidueLocator;
import org.roettig.PDBTools.SequenceAnchoredResidue;
import org.roettig.SequenceTools.ASCMSA;
import org.roettig.SequenceTools.HMM;
import org.roettig.SequenceTools.MSA;
import org.roettig.SequenceTools.SeqTools;
import org.roettig.SequenceTools.SequenceSet;
import org.roettig.SequenceTools.exception.FileParseErrorException;
import org.roettig.maths.statistics.Statistics;

public abstract class ASAFlow
{
	final public class ASAFlowContext
	{
		public String 								outputdir = "";	
		public ASAJob 								jd;
		public Structure 							struc;
		public Sequence 							template;
		public Map<Integer,SequenceAnchoredResidue> ASRidx;
		public List<SequenceAnchoredResidue> 		asrs;
		public List<Label>    						labels  = new Vector<Label>();
		public Map<Label,String> 					label_2_classname = new HashMap<Label,String>();
		public Map<String,String>  					intern_seq_name_2_orig_seqname = new HashMap<String,String>();
		public Map<String,Label> 					intern_seq_name_2_class = new HashMap<String,Label>();
		public SequenceSet      					coreseqs = new SequenceSet();
		public SequenceSet      					allseqs  = new SequenceSet();
		public Map<Sequence,Integer> 				seq2class = new HashMap<Sequence,Integer>();		
		public MSA              					allmsa   = null;
		public Map<Sequence,Double>  				seqdists = new HashMap<Sequence,Double>();
		public MSA              					coremsa;
		public HMM									corehmm;
		public SequenceSet      					sigseqs;
		public ASCMSA 								ascmsa;
		public InstanceContainer<PrimalInstance> 	samples = new DefaultInstanceContainer<PrimalInstance>();
		public InstanceContainer<DualInstance<String>> stringsamples = new DefaultInstanceContainer<DualInstance<String>>();
		public Model<?> 							model;
	}
	
	protected ASAFlowContext ctx = new ASAFlowContext();
	
	public ASAFlow(ASAJob _jd)
	{
		ctx.jd = _jd;
		ctx.outputdir = ctx.jd.getOutputDir();

		// disconnect standard log handlers
		clearLogHandlers();

		turnOnLogging();

		// connect ASCs own log handler
		Handler loghandler = new ASALogHandler(this);
		Logger.getLogger("org.roettig.SequenceTools.ThreeDCoffeeAlignment").addHandler(loghandler);
		Logger.getLogger("org.roettig.ML2.base.evaluation.ModelValidation").addHandler(loghandler);
		Logger.getLogger("org.roettig.PDBTools").addHandler(loghandler);


	}
	
	protected String lastStage;
	protected PrintStream logstream = System.out;
	
	public void log(String stage, String message)
	{		
		if(lastStage!=stage)
			logstream.println(String.format("%-20s",stage)+"|   "+message);
		else
			logstream.println(String.format("%-20s","")+"|   "+message);
		lastStage = stage;
	}

	public void turnOnLogging()
	{
		Logger.getLogger("").setLevel(Level.INFO);	
	}

	public void turnOffLogging()
	{
		Logger.getLogger("").setLevel(Level.OFF);
	}

	public void clearLogHandlers()
	{
		// clear previous log handlers
		Handler[] handlers = Logger.getLogger("").getHandlers();
		for(Handler h: handlers)
		{
			Logger.getLogger("").removeHandler(h);
		}	
	}
	
	public void exit(String message)
	{
		System.err.println(message);
		logstream.println(message);
		PrintWriter out = null;
		try
		{
			out =  new PrintWriter(ctx.outputdir+"/error");
		} 
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
		}
		out.println(message);
		out.close();
		System.exit(1);
	}
	
	public void run()
	{
		readInput();
		readStructure();
		extractASRs();
		readSequences();
		filterSequences();
		buildAlignment();
		extractSignatures();
		filterSignatures();
		prepareTrainingData();
		validateModel();
		analyzeImportance();
		buildFinalModel();
		store();
	}
	
	protected void readInput()
	{
		
	}

	protected void readStructure()
	{
		PDBFileReader pdbreader = new PDBFileReader();

		ctx.struc = null;
		try
		{
			ctx.struc = pdbreader.getStructure(ctx.jd.getPdbfile());
		}
		catch(Exception e)
		{
			exit("critical error: could not load PDB file "+ctx.jd.getPdbfile());
		}
		log("readStructure","loaded PDB structure "+ctx.struc.getPDBCode()+" from file "+ctx.jd.getPdbfile());   
		
		// extract template sequence from pdb structure
		ctx.template = PDBTools.getSequence(ctx.struc.getChain(0),null,null);
		
		ctx.jd.setTemplateSequence(ctx.template.seqString());
		log("readStructure","extracted sequence from first chain "+ctx.template.seqString());
	}

	protected void extractASRs()
	{
		List<ResidueLocator> locs = ctx.jd.getReslocs();
		ctx.asrs = PDBTools.getASCResidues(ctx.struc, locs, ctx.jd.getPositions(), ctx.jd.getDistance());

		List<Group> centers     = PDBTools.getGroups(ctx.struc, locs);

		List<Group> groups = PDBTools.getASCGroups(ctx.struc, locs, ctx.jd.getPositions(), ctx.jd.getDistance());

		PrintWriter out = null;
		try
		{
			out  = new PrintWriter(ctx.outputdir+"/activesite.pdb");
		} 
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
			e.printStackTrace(logstream);
			exit("critical error: could not store activesite.pdb");
		}

		PDBTools.writeAtomData(out, groups);  
		PDBTools.writeAtomData(out, centers);

		out.close();


		ctx.jd.setASCresidues(ctx.asrs);

		if(ctx.asrs.size()==0)
		{
			exit("critical error: no ASC residues found");
		}
		else
		{
			log("extractASCs","### Extracted ASCs ###");
			for(SequenceAnchoredResidue asc: ctx.asrs)
			{
				log("extractASCs",asc.toString());
			}
		}
		ctx.ASRidx = new TreeMap<Integer,SequenceAnchoredResidue>();
		for(SequenceAnchoredResidue aa: ctx.asrs)
		{
			ctx.ASRidx.put(aa.seqIdx, aa);
		}
		
	}
	
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
				lab =  new FactorLabel(ctx.jd.getClassnames().get(c-1),(double) c);
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
				//Sequence tmpseq = SeqTools.makeProteinSequence(c+"|"+seq.getName(),seq.seqString());
				Sequence tmpseq = SeqTools.makeProteinSequence(sid,seq.seqString());

				ctx.intern_seq_name_2_orig_seqname.put(sid,seq.getName());
				ctx.intern_seq_name_2_class.put(sid,lab);
				ctx.seq2class.put(tmpseq,c);
				ctx.allseqs.add(tmpseq);
				idx++;
			}
		}
	}

	protected ASAFlowInterceptor seqfilterInterceptor = new SequenceFilterInterceptor();
	
	protected void filterSequences()
	{
		seqfilterInterceptor.intercept(this, ctx);
	}

	protected ASAFlowInterceptor alignmentInterceptor = new AlignmentInterceptor(new MuscleAlignmentBuilder());
	
	public void setAlignmentInterceptor(AlignmentInterceptor ai)
	{
		alignmentInterceptor = ai;
	}
	
	protected void buildAlignment()
	{
		alignmentInterceptor.intercept(this, ctx);
	}

	protected void extractSignatures()
	{
		// retrieve ASC residues columns from MSA
		ctx.sigseqs = ctx.allmsa.getSubSequences(ctx.ASRidx.keySet(),"pdb"); 
		ctx.ascmsa = new ASCMSA(ctx.allmsa);
		ctx.ascmsa.setASCIdx(ctx.ASRidx.keySet());
		ctx.ascmsa.store(ctx.outputdir+"/all.asc");
		ctx.sigseqs.store(ctx.outputdir+"/sigs.afa");	
	}
	
	protected ASAFlowInterceptor sigfilterInterceptor = new SignatureFilterInterceptor(); 

	protected void filterSignatures()
	{
		sigfilterInterceptor.intercept(this, ctx);
	}

	
	protected PrimalEncoder<String> encoder = new WoldEncoder();
	
	protected void prepareTrainingData()
	{
		int n = 0;
		for(Sequence seq: ctx.sigseqs)
		{
			if(seq.getName().substring(0,3).equals("pdb"))
				continue;
			
			String sid = seq.getName();
			
			Label lab = ctx.intern_seq_name_2_class.get(sid);
			
			
			double[] fts = encoder.encode(seq.seqString());
			
			PrimalInstance pi = new PrimalInstance(lab,fts);
			pi.addProperty("sid",  seq.getName());
			pi.addProperty("osid", ctx.intern_seq_name_2_orig_seqname.get(sid));
			pi.addProperty("id", n);
			
			ctx.samples.add(pi);
			
			DualInstance<String> di = new DualInstance<String>(lab,seq.seqString());
			
			di.addProperty("sid",  seq.getName());
			di.addProperty("osid", ctx.intern_seq_name_2_orig_seqname.get(sid));
			di.addProperty("id", n);
			
			ctx.stringsamples.add(di);

			n++;
		}
	}

	
	
	
	protected void validateModel()
	{
		if(ctx.jd.getNumberOfClasses()==1)
		{
			log("validateModel","no model validation since extract-only mode");
			return;
		}

		
		RBFKernel k = new RBFKernel();
		double g_initial = RBFKernel.estimateGamma(ctx.samples);
		
		k.setGamma(g_initial/10.0,g_initial/4.0,g_initial/2.0,g_initial,g_initial*2,g_initial*4,g_initial*10);
		
		CSVCModel<PrimalInstance> m = new CSVCModel<PrimalInstance>(k);
		m.setC(0.1,1.0,10.0,100.0,1000.0);
		
		
		SelectedModel<PrimalInstance> bestModel = null;
		log("validateModel","### Doing Nested CV ####");

		//int nInner = ctx.jd.getNInnerFolds();
		//int nOuter = ctx.jd.getNOuterFolds();
		//

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
		ctx.jd.setFmeasure(bestModel.qual);
		
		List<Label> yp = null; List<Label> yt = null;
		Prediction.split(bestModel.predictions, yt, yp);
		
		for(Label lab: ctx.labels)
		{
			List<Double> prec_rec = Statistics.calcPrecRec(lab, yp, yt);
			log("###",String.format(Locale.ENGLISH,"prec(%s)=%.3f rec(%s)=%.3f",lab.toString(),prec_rec.get(0),lab.toString(),prec_rec.get(1)));
		}
		/*
		Map<Object, List<Double>> stats = Statistics.calc

		Map<String,Double> precs = new HashMap<String,Double>();
		Map<String,Double> recs  = new HashMap<String,Double>();

		for(Object o: stats.keySet())
		{
			log("###",String.format(Locale.ENGLISH,"prec(%s)=%.3f rec(%s)=%.3f",o,stats.get(o).get(0),o,stats.get(o).get(1)));
			String classname = ctx.label_2_classname.get((Label) o);
			precs.put(classname,stats.get(o).get(0));
			recs.put(classname,stats.get(o).get(1));
		}
		ctx.jd.setPrecs(precs);
		ctx.jd.setRecs(recs);
		*/
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

	protected abstract void analyzeImportance();


	protected void buildFinalModel()
	{		
	}

	protected void store()
	{
		SequenceSet raus = new SequenceSet();
		for(Sequence s: ctx.sigseqs)
		{
			if(!ctx.intern_seq_name_2_orig_seqname.containsKey(s.getName()))
			{
				raus.add(s);
			}
			else
			{
				Sequence st = SeqTools.makeProteinSequence(ctx.intern_seq_name_2_orig_seqname.get(s.getName()).toString(), s.seqString());
				raus.add(st);
			}
		}
		raus.store(ctx.outputdir+"/sig-orig.afa");

		SequenceSet raus2 = new SequenceSet();
		for(Sequence s: ctx.allmsa)
		{
			if(!ctx.intern_seq_name_2_orig_seqname.containsKey(s.getName()))
			{
				raus2.add(s);
			}
			else
			{
				Sequence st = SeqTools.makeProteinSequence(ctx.intern_seq_name_2_orig_seqname.get(s.getName()).toString(), s.seqString());
				raus2.add(st);
			}	    
		}
		raus2.store(ctx.outputdir+"/all-orig.afa");

		try
		{
			ctx.jd.save(ctx.outputdir+"/ASCresult");
		} 
		catch (Exception e)
		{
			e.printStackTrace();
		}	
	}
	 
	
	public Model<?> getModel()
	{
		return ctx.model;
	}
	
	public MSA getMSA()
	{
		return ctx.allmsa; 
	}
}
