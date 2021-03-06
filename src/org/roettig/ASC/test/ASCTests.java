package org.roettig.ASC.test;

import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;
import org.junit.Test;
import org.roettig.ASC.asc.ASAFlow;
import org.roettig.ASC.asc.ASAJob;
import org.roettig.ASC.asc.ASC;
import org.roettig.ASC.asc.ASR;
import org.roettig.ASC.test.data.Resources;
import org.roettig.PDBTools.ResidueLocator;

public class ASCTests extends TestCase
{

	@Test
	public void testASCNRPS() throws Exception
	{
		ASAJob jd = new ASAJob();
		
		String tmpdir = System.getProperty("java.io.tmpdir");
		
		jd.setPdbfile("data/1AMU.pdb");
		jd.setOutputDir(tmpdir);
		
		ResidueLocator       resloc  = new ResidueLocator("A","PHE",566);
		List<ResidueLocator> reslocs = new ArrayList<ResidueLocator>();
		reslocs.add(resloc);
		
		jd.setReslocs(reslocs);
		jd.setDistance(8.0);
		
		List<String> filenames = new ArrayList<String>();
		
		
		filenames.add("data/ala.fa");
		filenames.add("data/tyr.fa");
		filenames.add("data/asp.fa");
		
		jd.setFilenames(filenames);
		List<String> classnames = new ArrayList<String>();
		classnames.add("ala");
		classnames.add("tyr");
		classnames.add("asp");
		
		jd.setClassnames(classnames);
		
		ASAFlow flow = new ASC(jd);
		flow.run();
		
		assertEquals(Resources.readString(Resources.class.getResourceAsStream("activesite.pdb")),Resources.readString(new FileInputStream(tmpdir+"/activesite.pdb")));
		assertEquals(Resources.readString(Resources.class.getResourceAsStream("all.asc")),Resources.readString(new FileInputStream(tmpdir+"/all.asc")));
		
		assertEquals(0.94118,jd.getPrecs().get("ala"),1e-4);
		assertEquals(0.92308,jd.getPrecs().get("tyr"),1e-4);
		assertEquals(1.00000,jd.getPrecs().get("asp"),1e-4);
		
		assertEquals(0.94118,jd.getRecs().get("ala"),1e-4);
		assertEquals(0.92308,jd.getRecs().get("tyr"),1e-4);
		assertEquals(1.00000,jd.getRecs().get("asp"),1e-4);
		
		assertEquals(0.954751131221719,jd.getQuality(),1e-4);
		
		// check position importance
		
		// 0 1: ala vs. tyr
		assertEquals(1,jd.getRanks(0, 1)[9]);
		assertEquals(2,jd.getRanks(0, 1)[5]);
		assertEquals(3,jd.getRanks(0, 1)[14]);
		
		// 0 2: ala vs. asp
		assertEquals(1,jd.getRanks(0, 2)[19]);
		assertEquals(2,jd.getRanks(0, 2)[3]);
		assertEquals(3,jd.getRanks(0, 2)[33]);
		
		// 1 2: tyr vs. asp
		assertEquals(1,jd.getRanks(1, 2)[3]);
		assertEquals(2,jd.getRanks(1, 2)[5]);
		assertEquals(3,jd.getRanks(1, 2)[9]);
		
		
		jd.save(tmpdir+"/job.ini");
	}
	
	
	@Test
	public void testASRNRPS() throws Exception
	{
		ASAJob jd = new ASAJob();
		
		String tmpdir = System.getProperty("java.io.tmpdir");
		
		jd.setPdbfile("data/1AMU.pdb");
		jd.setOutputDir(tmpdir);
		
		ResidueLocator       resloc  = new ResidueLocator("A","PHE",566);
		List<ResidueLocator> reslocs = new ArrayList<ResidueLocator>();
		reslocs.add(resloc);
		
		jd.setReslocs(reslocs);
		jd.setDistance(8.0);
		
		List<String> filenames = new ArrayList<String>();
		
		filenames.add("data/nrpsR.fa");		
		jd.setFilenames(filenames);
		
		ASAFlow flow = new ASR(jd);
		
		flow.run();
		
		assertEquals(0.8373376080011307,jd.getQuality(),1e-4);
		

		// check position importance
		
		// 0 1
		assertEquals(1,jd.getRanks(0, 1)[33]);
		assertEquals(2,jd.getRanks(0, 1)[15]);
		assertEquals(3,jd.getRanks(0, 1)[3]);

		
		jd.save(tmpdir+"/job.ini");
	}
}
