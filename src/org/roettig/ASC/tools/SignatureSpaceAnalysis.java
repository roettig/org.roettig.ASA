/**
 * 
 */
package org.roettig.ASC.tools;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Locale;
import org.biojava.bio.seq.Sequence;
import org.roettig.SequenceTools.SequenceSet;
import org.roettig.SequenceTools.exception.FileParseErrorException;

/**
 * @author roettig
 *
 */
public class SignatureSpaceAnalysis
{

	public static void doDistanceAnalysis(SequenceSet x, SequenceSet y, SimilarityDelegate sim)
	{
		for(Sequence sx: x)
		{
			String ssx = sx.seqString();
			double maxSim = Double.NEGATIVE_INFINITY;
			double selfsimX = sim.getSimilarity(ssx, ssx);
			String bestHit="";

			for(Sequence sy: y)
			{
				String ssy = sy.seqString();
				double selfsimY = sim.getSimilarity(ssy, ssy);
				double similarity = sim.getSimilarity(ssx, ssy)/Math.sqrt(selfsimX*selfsimY);
				if(similarity>maxSim)
				{
					maxSim  = similarity;
					bestHit = ssy; 
				}
			}
			System.out.println(String.format(Locale.ENGLISH,"%-7.3f %s-%s %s",maxSim,sx.seqString(),bestHit,sx.getName()));
		}
	}

	public static void calcDistanceMatrix(SequenceSet seqs, BufferedWriter writer, SimilarityDelegate sim) throws IOException
	{
		for(Sequence sx: seqs)
		{
			String ssx  = sx.seqString();
			String rssx = String.format("%c%c%c%c",ssx.charAt(3),ssx.charAt(4),ssx.charAt(7),ssx.charAt(17)); 
			//System.out.println(rssx);
			double selfsimX = sim.getSimilarity(rssx, rssx);


			for(Sequence sy: seqs)
			{
				String ssy = sy.seqString();
				String rssy = String.format("%c%c%c%c",ssy.charAt(3),ssy.charAt(4),ssy.charAt(7),ssy.charAt(17));
				double selfsimY = sim.getSimilarity(rssy, rssy);
				double similarity = sim.getSimilarity(rssx, rssy)/Math.sqrt(selfsimX*selfsimY);
				writer.write(String.format(Locale.ENGLISH,"%.5f ",similarity));
			}
			writer.write("\n");

		}
	}

	/**
	 * @param args
	 * @throws FileParseErrorException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws FileParseErrorException, IOException
	{

		SequenceSet knownSeqs   = SequenceSet.readFromFile("/home/roettig/coops/nrps/nrps_unknown/sigs-known.afa");
		//SequenceSet unknownSeqs = SequenceSet.readFromFile("/home/roettig/coops/nrps/nrps_unknown/sigs-unknown.afa");
		SequenceSet unknownSeqs = SequenceSet.readFromFile("/home/roettig/coops/wohlleben2/bact/big1/sig-orig.afa");
		doDistanceAnalysis(unknownSeqs,knownSeqs, new BLOSUMSimilarityDelegate());

		/*
	System.out.println("haha");
	BufferedWriter out = new BufferedWriter(new FileWriter("/tmp/mat-chemsim-red")); 
	SequenceSet allSeqs   = SequenceSet.readFromFile("/tmp/allsigs.afa");
	//calcDistanceMatrix(allSeqs,out,new BLOSUMSimilarityDelegate());
	calcDistanceMatrix(allSeqs,out,new ChemicalSimilarityDelegate());
	out.close();
		 */
	}

}
