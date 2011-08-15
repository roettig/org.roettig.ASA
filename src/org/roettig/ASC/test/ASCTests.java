package org.roettig.ASC.test;

import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;
import org.junit.Test;
import org.roettig.ASC.asc.ASAFlow;
import org.roettig.ASC.asc.ASAJob;
import org.roettig.ASC.asc.ASC;
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
		
		jd.save(tmpdir+"/job.ini");
	}

}
