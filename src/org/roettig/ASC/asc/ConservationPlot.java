/**
 * 
 */
package org.roettig.ASC.asc;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.PDBFileReader;
import org.roettig.PDBTools.PDBTools;
import org.roettig.SequenceTools.MSA;
import org.roettig.SequenceTools.SeqTools;
import org.roettig.SequenceTools.SequenceSet;


/**
 * @author roettig
 *
 */
public class ConservationPlot
{

	private static Structure struc    = null;
	private static Sequence  template = null;
	private static String pdbfile;

	public static void main(String[] args) throws Exception
	{
		pdbfile = args[0];

		readStructure();

		SequenceSet seqs = SequenceSet.readFromFile(args[1]);

		String queryId = seqs.getIDs().get(0);

		seqs.add(template);

		MSA msa = MSA.createMuscleMSA(seqs);

		msa.store("/tmp/msa.afa");

		Sequence tmpl = msa.getById("pdb");
		Sequence qry  = msa.getById(queryId);
		
		System.out.println("load "+pdbfile);
		System.out.println("show cartoon, chain A\nhide lines, chain A\nhide (resnam HOH)\n");
		System.out.println("select cons, resi 9999999\n");
		System.out.println("select indel, resi 9999999\n");
		System.out.println("select div, resi 9999999\n");
		String t = tmpl.seqString();
		String q = qry.seqString();
		int offset = 1;
		int pos    = 1;
		for(int i=0;i<t.length();i++)
		{
			char tC = t.charAt(i);
			char qC = q.charAt(i);

			// template deletion
			if(tC=='-')
			{
				continue;
			}

			//System.out.println(tC+" - "+qC+" "+(pos+offset));
			if(tC==qC)
			{
				System.out.println("select cons, cons or resi "+(pos+offset));
			}
			else
			{
				System.out.println("select div, div or resi "+(pos+offset));
			}
			pos++;
		}
		System.out.println("color green, cons");
		System.out.println("color white, div");


	}



	protected static void readStructure()
	{
		PDBFileReader pdbreader = new PDBFileReader();

		struc = null;
		try
		{
			struc = pdbreader.getStructure(pdbfile);
		}
		catch(Exception e)
		{
			System.out.println("critical error: could not load PDB file "+pdbfile);
			System.exit(1);
		}

		Sequence seq = PDBTools.getSequence(struc.getChain(0),null,null);
		template = SeqTools.makeProteinSequence("pdb", seq.seqString());
	}

}
