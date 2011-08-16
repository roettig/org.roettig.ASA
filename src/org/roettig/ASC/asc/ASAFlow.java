package org.roettig.ASC.asc;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
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
import org.roettig.ASC.asc.interceptors.DefaultSequenceFilterInterceptor;
import org.roettig.ASC.asc.interceptors.SignatureFilterInterceptor;
import org.roettig.ASC.kernels.WoldEncoder;
import org.roettig.MLToolbox.base.PrimalEncoder;
import org.roettig.MLToolbox.base.impl.DefaultInstanceContainer;
import org.roettig.MLToolbox.base.instance.DualInstance;
import org.roettig.MLToolbox.base.instance.InstanceContainer;
import org.roettig.MLToolbox.base.instance.PrimalInstance;
import org.roettig.MLToolbox.base.label.FactorLabel;
import org.roettig.MLToolbox.base.label.Label;
import org.roettig.MLToolbox.model.Model;
import org.roettig.PDBTools.PDBTools;
import org.roettig.PDBTools.ResidueLocator;
import org.roettig.PDBTools.SequenceAnchoredResidue;
import org.roettig.SequenceTools.ASCMSA;
import org.roettig.SequenceTools.HMM;
import org.roettig.SequenceTools.MSA;
import org.roettig.SequenceTools.SeqTools;
import org.roettig.SequenceTools.SequenceSet;
import org.roettig.SequenceTools.exception.FileParseErrorException;


public abstract class ASAFlow
{
	/**
	 * The class ASAFlowContext stores contextual information about the flow which is communicated
	 * to any ASAFlowInterceptor and can be modified by these interceptors.
	 * 
	 * @author roettig
	 *
	 */
	final public class ASAFlowContext
	{
		public String 									outputdir="";	
		public ASAJob 									jd;
		public Structure 								struc;
		public Sequence 								template;
		public Map<Integer,SequenceAnchoredResidue> 	ASRidx;
		public List<SequenceAnchoredResidue> 			asrs;
		public List<Label>    							labels  = new Vector<Label>();
		public Map<Label,String> 						label_2_classname = new HashMap<Label,String>();
		public Map<String,String>  						intern_seq_name_2_orig_seqname = new HashMap<String,String>();
		public Map<String,Label> 						intern_seq_name_2_class = new HashMap<String,Label>();
		public SequenceSet      						coreseqs = new SequenceSet();
		public SequenceSet      						allseqs  = new SequenceSet();
		public Map<Sequence,Integer> 					seq2class = new HashMap<Sequence,Integer>();		
		public MSA              						allmsa   = null;
		public Map<Sequence,Double>  					seqdists = new HashMap<Sequence,Double>();
		public MSA              						coremsa;
		public HMM										corehmm;
		public SequenceSet      						sigseqs;
		public ASCMSA 									ascmsa;
		public InstanceContainer<PrimalInstance> 		samples       = new DefaultInstanceContainer<PrimalInstance>();
		public InstanceContainer<DualInstance<String>> 	stringsamples = new DefaultInstanceContainer<DualInstance<String>>();
		public Model<?> 								model;
	}
	
	protected ASAFlowContext ctx = new ASAFlowContext();
	protected String lastStage;
	protected PrintStream logstream = System.out;
	
	protected PrimalEncoder<String> encoder = new WoldEncoder();
	
	// interceptors
	protected ASAFlowInterceptor seqfilterInterceptor = new DefaultSequenceFilterInterceptor();
	protected ASAFlowInterceptor alignmentInterceptor = new AlignmentInterceptor(new MuscleAlignmentBuilder());
	protected ASAFlowInterceptor sigfilterInterceptor = new SignatureFilterInterceptor();
	
	/**
	 * ctor with supplied job data.
	 * 
	 * @param _jd job data
	 */
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
		Logger.getLogger("org.roettig.MLToolbox.validation.ModelValidation").addHandler(loghandler);
		Logger.getLogger("org.roettig.PDBTools").addHandler(loghandler);


	}
	
	/**
	 * sends out log message for current stage.
	 * 
	 * @param stage
	 * @param message
	 */
	public void log(String stage, String message)
	{		
		if(lastStage!=stage)
			logstream.println(String.format("%-20s",stage)+"|   "+message);
		else
			logstream.println(String.format("%-20s","")+"|   "+message);
		lastStage = stage;
	}

	/**
	 * turns on logging.
	 */
	public void turnOnLogging()
	{
		Logger.getLogger("").setLevel(Level.INFO);	
	}

	/**
	 * turns off any logging.
	 */
	public void turnOffLogging()
	{
		Logger.getLogger("").setLevel(Level.OFF);
	}

	/**
	 * detach any log handlers.
	 */
	public void clearLogHandlers()
	{
		// clear previous log handlers
		Handler[] handlers = Logger.getLogger("").getHandlers();
		for(Handler h: handlers)
		{
			Logger.getLogger("").removeHandler(h);
		}	
	}
	
	/**
	 * terminates the execution of this flow with the supplied error message.
	 * 
	 * @param message
	 */
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
	
	/**
	 * read any input data.
	 */
	protected void readInput()
	{	
	}

	/**
	 * reads in the structure from a PDB file.
	 */
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

	/**
	 * extract active site residues (ASRs) from template structure.
	 */
	protected void extractASRs()
	{
		// get residue locators of ASRs from context
		List<ResidueLocator> locs = ctx.jd.getReslocs();
		// and extract ASRs from structre
		ctx.asrs = PDBTools.getASCResidues(ctx.struc, locs, ctx.jd.getPositions(), ctx.jd.getDistance());

		List<Group> centers = PDBTools.getGroups(ctx.struc, locs);

		List<Group> groups  = PDBTools.getASCGroups(ctx.struc, locs, ctx.jd.getPositions(), ctx.jd.getDistance());

		// write out active site atoms to pdb file
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
	
	/**
	 * read training sequences from files.
	 */
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

				Sequence tmpseq = SeqTools.makeProteinSequence(sid,seq.seqString());

				ctx.intern_seq_name_2_orig_seqname.put(sid,seq.getName());
				ctx.intern_seq_name_2_class.put(sid,lab);
				ctx.seq2class.put(tmpseq,c);
				ctx.allseqs.add(tmpseq);
				idx++;
			}
		}
	}

	public void setSequenceFilterInterceptor(ASAFlowInterceptor ai)
	{
		seqfilterInterceptor = ai;
	}
	
	/**
	 * filter sequences.
	 */
	protected void filterSequences()
	{
		seqfilterInterceptor.intercept(this, ctx);
	}

	public void setAlignmentInterceptor(ASAFlowInterceptor ai)
	{
		alignmentInterceptor = ai;
	}
	
	/**
	 * build multiple sequence alignment (MSA) of all training sequences and template.
	 */
	protected void buildAlignment()
	{
		alignmentInterceptor.intercept(this, ctx);
	}

	/**
	 * extract signatures from all training sequences using msa.
	 */
	protected void extractSignatures()
	{
		// retrieve ASC residues columns from MSA
		ctx.sigseqs = ctx.allmsa.getSubSequences(ctx.ASRidx.keySet(),"pdb"); 
		ctx.ascmsa = new ASCMSA(ctx.allmsa);
		ctx.ascmsa.setASCIdx(ctx.ASRidx.keySet());
		ctx.ascmsa.store(ctx.outputdir+"/all.asc");
		ctx.sigseqs.store(ctx.outputdir+"/sigs.afa");	
	}
	
	 

	/**
	 * filter signature sequences.
	 */
	protected void filterSignatures()
	{
		sigfilterInterceptor.intercept(this, ctx);
	}

	
	
	
	/**
	 * prepare training data : encode training sequences into PrimalInstances, etc.
	 */
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

	protected abstract void validateModel();
	
	protected abstract void analyzeImportance();

	protected abstract void buildFinalModel();
	

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
