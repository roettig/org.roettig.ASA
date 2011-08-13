package org.roettig.ASC.asc;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.biojava.bio.structure.Group;

import org.roettig.PDBTools.ResidueLocator;
import org.roettig.PDBTools.SequenceAnchoredResidue;
import org.roettig.PDBTools.Vector3;
import org.roettig.maths.statistics.Statistics;
import org.roettig.maths.util.Pair;




/**
 * @author roettig
 *
 */
public class ASAJob implements Serializable
{
	private static final long serialVersionUID = 5931979496857071036L;

	private String modeltype      = null;
	private String kerneltype     = null;
	private String THREEDCOFFEEID = null;
	private String MSAmethod      = "muscle";
	private String templatesequence;
	private String jobid;
	private double fmeasure;
	private int numberOfClasses;
	private Map<String,Double> seqids  = new HashMap<String,Double>();
	private Map<String,Double> precs   = new HashMap<String,Double>();
	private Map<String,Double> recs    = new HashMap<String,Double>();
	private Map<String,Integer> occupancies  = new HashMap<String,Integer>();
	private List<String> classnames  = new Vector<String>();
	private Map<String,String> sequenceclasses = new HashMap<String,String>();

	private double distance;
	private String templatename;
	private List<SequenceAnchoredResidue> ascresidues = new Vector<SequenceAnchoredResidue>();
	private List<ResidueLocator> reslocs = new Vector<ResidueLocator>();

	private double pid_threshold      = 0.20;
	private String pdbfile            = null; 
	private String outputdir	      = "";


	private List<String> filenames  = null;
	private double mttsi              = 1.0;
	private double frac               = 0.5;
	private boolean usesigs           = true;
	private int nInnerFolds           = 3;
	private int nOuterFolds           = 3;
	private List<Vector3> positions = new Vector<Vector3>();

	private Map<Integer, Map<Integer,double[]> > scores = new HashMap<Integer, Map<Integer,double[]> >();

	Map<String, List<Double> > hyperparams;




	public String getMSAmethod()
	{
		return MSAmethod;
	}


	public void setMSAmethod(String mSAmethod)
	{
		MSAmethod = mSAmethod;
	}


	public Map<String, String> getSequenceclasses()
	{
		return sequenceclasses;
	}


	public void setSequenceclasses(Map<String, String> sequenceclasses)
	{
		this.sequenceclasses = sequenceclasses;
	}


	public Map<String, Double> getSeqIds()
	{
		return seqids;
	}


	public void setSeqIds(Map<String, Double> seqids)
	{
		this.seqids = seqids;
	}


	public boolean hasThreeDCoffeeId()
	{
		return (THREEDCOFFEEID!=null);
	}


	public String getThreeDCoffeeId()
	{
		return THREEDCOFFEEID;
	}

	public void setThreeDCoffeeId(String tHREEDCOFFEEID)
	{
		THREEDCOFFEEID = tHREEDCOFFEEID;
	}

	public Map<String, List<Double>> getHyperParams()
	{
		return hyperparams;
	}

	public void setHyperParams(Map<String, List<Double>> params)
	{
		this.hyperparams = params;
	}

	public String getModelType()
	{
		return modeltype;
	}

	public void setModelType(String modeltype)
	{
		this.modeltype = modeltype;
	}

	public void setScore(Integer c1, Integer c2, double[] score)
	{
		if(!scores.containsKey(c1))
		{
			scores.put(c1, new HashMap<Integer,double[]>());
		}
		scores.get(c1).put(c2, score);
	}

	public double[] getScore(Integer c1, Integer c2)
	{
		return scores.get(c1).get(c2);
	}

	public int[] getRanks(Integer c1, Integer c2)
	{
		double scores[] = getScore(c1,c2);
		List<Double> data = new ArrayList<Double>();
		for(Double d:scores)
		{
			data.add(d);
		}
		int ranks[] = new int[scores.length];
		List<Pair<Integer,Double>> sorted = Statistics.sort(data, false);

		int rnk = 1;
		for(Pair<Integer,Double> p: sorted)
		{
			ranks[p.getFirst()]=rnk;
			rnk++;
		}
		return ranks;
	}

	public String getTemplateSequence()
	{
		return templatesequence;
	}

	public void setTemplateSequence(String templatesequence)
	{
		this.templatesequence = templatesequence;
	}

	public String getOutputDir()
	{
		return outputdir;
	}

	public void setOutputDir(String outputdir)
	{
		this.outputdir = outputdir;
	}

	public int getNInnerFolds()
	{
		return nInnerFolds;
	}

	public void setNInnerFolds(int nInnerFolds)
	{
		this.nInnerFolds = nInnerFolds;
	}

	public int getNOuterFolds()
	{
		return nOuterFolds;
	}

	public void setNOuterFolds(int nOuterFolds)
	{
		this.nOuterFolds = nOuterFolds;
	}

	public List<ResidueLocator> getReslocs()
	{
		return reslocs;
	}

	public void setReslocs(List<ResidueLocator> reslocs)
	{
		this.reslocs = reslocs;
	}

	public double getPid_threshold()
	{
		return pid_threshold;
	}

	public void setPid_threshold(double pidThreshold)
	{
		pid_threshold = pidThreshold;
	}

	public String getPdbfile()
	{
		return pdbfile;
	}

	public void setPdbfile(String pdbfile)
	{
		this.pdbfile = pdbfile;
	}

	public List<String> getFilenames()
	{
		return filenames;
	}

	public void setFilenames(List<String> filenames)
	{
		this.filenames = filenames;
	}

	public double getMttsi()
	{
		return mttsi;
	}

	public void setMttsi(double mttsi)
	{
		this.mttsi = mttsi;
	}

	public double getFrac()
	{
		return frac;
	}

	public void setFrac(double frac)
	{
		this.frac = frac;
	}

	public boolean isUsesigs()
	{
		return usesigs;
	}

	public void setUsesigs(boolean usesigs)
	{
		this.usesigs = usesigs;
	}

	public String getKernelType()
	{
		return kerneltype;
	}

	public void setKernelType(String kernel)
	{
		kerneltype = kernel;
	}

	public List<Vector3> getPositions()
	{
		return positions;
	}

	public void setPositions(List<Vector3> positions)
	{
		this.positions = positions;
	}

	public List<SequenceAnchoredResidue> getASCresidues()
	{
		return ascresidues;
	}

	public void setASCresidues(List<SequenceAnchoredResidue> ascresidues)
	{
		this.ascresidues = ascresidues;
	}

	public Map<String, Integer> getOccupancies()
	{
		return occupancies;
	}

	public void setOccupancies(Map<String, Integer> occupancies)
	{
		this.occupancies = occupancies;
	}

	public String getJobid()
	{
		return jobid;
	}

	public void setJobid(String jobid)
	{
		this.jobid = jobid;
	}

	public int getNumberOfClasses()
	{
		return numberOfClasses;
	}

	public void setNumberOfClasses(int numberOfClasses)
	{
		this.numberOfClasses = numberOfClasses;
	}

	public Map<String, Double> getPrecs()
	{
		return precs;
	}

	public void setPrecs(Map<String, Double> precs)
	{
		this.precs = precs;
	}

	public Map<String, Double> getRecs()
	{
		return recs;
	}

	public void setRecs(Map<String, Double> recs)
	{
		this.recs = recs;
	}

	public List<String> getClassnames()
	{
		return classnames;
	}

	public void setClassnames(List<String> classnames)
	{
		this.classnames = classnames;
		this.setNumberOfClasses(classnames.size());
	}

	public double getDistance()
	{
		return distance;
	}

	public void setDistance(double distance)
	{
		this.distance = distance;
	}

	public String getTemplateName()
	{
		return templatename;
	}

	public void setTemplateName(String template)
	{
		this.templatename = template;
	}

	public ASAJob()
	{

	}

	public double getFmeasure()
	{
		return fmeasure;
	}

	public void setFmeasure(double fmeasure)
	{
		this.fmeasure = fmeasure;
	}

	public String getTable()
	{
		return "<table border='1'><tr><td>1</td><td>2</td></tr><tr><td>3</td><td>4</td></tr></table>";
	}

	public void save(String filename) throws Exception
	{
		FileOutputStream f_out = new FileOutputStream(filename);
		// Write object with ObjectOutputStream
		ObjectOutputStream obj_out = new ObjectOutputStream (f_out);
		// Write object out to disk
		obj_out.writeObject( this );
		f_out.close();
	}

	public static ASAJob load(String filename) throws Exception
	{
		// restore object from file ...
		// Read from disk using FileInputStream
		FileInputStream f_in = new FileInputStream(filename);

		// Read object using ObjectInputStream
		ObjectInputStream obj_in = new ObjectInputStream (f_in);
		ASAJob obj = (ASAJob) obj_in.readObject();
		f_in.close();
		return obj;
	}
}
