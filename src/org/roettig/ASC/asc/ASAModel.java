/**
 * 
 */
package org.roettig.ASC.asc;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;

import org.roettig.MLToolbox.model.Model;
import org.roettig.PDBTools.SequenceAnchoredResidue;
import org.roettig.SequenceTools.MSA;

/**
 * @author roettig
 *
 */
public class ASAModel implements Serializable
{
	private static final long serialVersionUID = -3918205273629442734L;

	private Model model;
	private MSA   msa;
	private List<SequenceAnchoredResidue> ascresidues = new Vector<SequenceAnchoredResidue>();

	public ASAModel(Model _model, MSA _msa)
	{
		model = _model;
		msa   = _msa;
	}

	public Model getModel()
	{
		return model;
	}
	public void setModel(Model model)
	{
		this.model = model;
	}
	public MSA getMSA()
	{
		return msa;
	}

	public void setMSA(MSA msa)
	{
		this.msa = msa;
	}

	public List<SequenceAnchoredResidue> getASCresidues()
	{
		return ascresidues;
	}

	public void setASCresidues(List<SequenceAnchoredResidue> ascresidues)
	{
		this.ascresidues = ascresidues;
	}

	public Set<Integer> getASCresidueIndices()
	{
		Set<Integer> ret = new TreeSet<Integer>();
		for(SequenceAnchoredResidue aa: getASCresidues())
		{
			ret.add(aa.seqIdx);
		}
		return ret;
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

	public static ASAModel load(String filename) throws Exception
	{
		// restore object from file ...
		// Read from disk using FileInputStream
		FileInputStream f_in = new FileInputStream(filename);

		// Read object using ObjectInputStream
		ObjectInputStream obj_in = new ObjectInputStream (f_in);
		ASAModel ret = (ASAModel) obj_in.readObject();
		f_in.close();
		return ret;
	}
}
